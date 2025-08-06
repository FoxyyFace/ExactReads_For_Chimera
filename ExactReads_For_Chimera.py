import os
import logging
import time
import psutil
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from Geneinfo_Exactreads import GeneInfo, ReadIndex, ExactReads

ALL_start_time = time.time()

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("processing.log"),
        logging.StreamHandler()
    ]
)

# ###########################
# Configuration Parameters
# ###########################

# Get available CPU cores
TOTAL_CORES = psutil.cpu_count(logical=True)
BUILD_INDEX_THREADS = min(8, TOTAL_CORES)      # Threads per index build, less than half of ALL cores,
INDEX_SEARCH_THREADS = min(16, TOTAL_CORES)   # Threads per index search
SAMPLE_PROCESS_POOL = min(8, TOTAL_CORES)     # Max concurrent samples

logging.info(f"System Configuration: {TOTAL_CORES} logical cores")
logging.info(f"Index build threads per sample: {BUILD_INDEX_THREADS}")
logging.info(f"Index search threads per sample: {INDEX_SEARCH_THREADS}")
logging.info(f"Max concurrent samples: {SAMPLE_PROCESS_POOL}")

# Path configurations
seqtk_path = "/PATH/TO/seqtk/seqtk"
RNA1 = {
    "name": "EFV83_RS15570:EFV83_RS15625",
    "type": "IGR",
    "strand": 0
}
RNA2 = {
    "name": "PA14_RS29125",
    "type": "gene",
    "strand": 0
}
UTR_LENGTH = 150
SAMPLE_GROUPS = {
    ("Mabin-1-LHF18430_L5", "MabcoPA14_IN"),
    ("Mabin-2-LHF18431_L5", "MabcoPA14_IN"),
    ("Mabin-3-LHF18432_L5", "MabcoPA14_IN")
}
GFF_FILE = "/PATH/TO/PA14_Mab.gff"
BAM_DIR = "/PATH/TO/250629-LH00209B"    # 存放bam文件和fastq文件的目录
FASTQ_BAM_PREFIX = "trimmed_"
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Output format
OUTPUT_FORMAT = "fasta"    # 或者fastq

# ###########################
# Core Functions
# ###########################

def build_index(sample_info):
    """Parallel index builder"""
    bam_file = sample_info["bam_file"]
    index_file = f"{bam_file}.ridx"
    
    if os.path.exists(index_file):
        logging.info(f"Using existing index for {bam_file}")
        return True

    try:
        ReadIndex.build_index(
            bam_file,
            index_file,
            num_threads=BUILD_INDEX_THREADS,
            bam_threads=BUILD_INDEX_THREADS,
            chunk_size=50000
        )
        return True
    except Exception as e:
        logging.error(f"Index build failed for {bam_file}: {str(e)}")
        if os.path.exists(index_file):
            os.remove(index_file)
        return False

def process_sample(sample_info, rna1_info, rna2_info):
    """Parallel sample processor"""
    try:
        sample_name = sample_info["sample_name"]
        bam_file = sample_info["bam_file"]
        
        # Load index
        index_data = ReadIndex.load_index(f"{bam_file}.ridx")
        if not index_data:
            raise RuntimeError("Index load failed")

        # Parallel search
        logging.info(f"Searching chimeric reads in {sample_name} with {INDEX_SEARCH_THREADS} threads")
        chimeric_ids = ReadIndex.find_in_index(
            index_data,
            rna1_info,
            rna2_info,
            num_threads=INDEX_SEARCH_THREADS
        )

        if not chimeric_ids:
            logging.warning(f"No chimeric reads found in {sample_name}")
            return sample_name, False, 0, "No chimeric reads", None

        # Prepare output
        output_file = os.path.join(OUTPUT_DIR, f"{sample_name}_chimeric.fasta")

        # Extract sequences
        fastq_files = (
            (sample_info["fastq_r1"], sample_info["fastq_r2"]) 
            if sample_info["is_paired_end"] 
            else sample_info["fastq"]
        )

        success = ExactReads.extract_reads_with_seqtk(
            fastq_files=fastq_files,
            read_ids=chimeric_ids,
            output_file=output_file,
            seqtk_path=seqtk_path,
            output_format=OUTPUT_FORMAT
        )

        if success and os.path.getsize(output_file) > 0:
            return sample_name, True, len(chimeric_ids), "Success", output_file
        else:
            raise RuntimeError("Extraction failed")
    except Exception as e:
        logging.error(f"Sample {sample_name} failed: {str(e)}")
        return sample_name, False, 0, str(e), None

# ###########################
# Main Workflow
# ###########################

def main():
    # Initialize
    sample_info_list = []
    group_files = defaultdict(list)

    # Prepare samples
    logging.info("Preparing sample information...")
    for sample_name, group_name in SAMPLE_GROUPS:
        bam_file = os.path.join(BAM_DIR, f"{FASTQ_BAM_PREFIX}{sample_name}.bam")
        r1 = os.path.join(BAM_DIR, f"{FASTQ_BAM_PREFIX}{sample_name}_1.fq.gz")
        r2 = os.path.join(BAM_DIR, f"{FASTQ_BAM_PREFIX}{sample_name}_2.fq.gz")
        
        if all(os.path.exists(f) for f in [bam_file, r1, r2]):
            sample_info_list.append({
                "sample_name": sample_name,
                "group_name": group_name,
                "bam_file": bam_file,
                "is_paired_end": True,
                "fastq_r1": r1,
                "fastq_r2": r2
            })
            logging.info(f"Sample '{sample_name}' files verified")
        else:
            logging.warning(f"Missing files for sample '{sample_name}'")

    if not sample_info_list:
        logging.error("No valid samples found!")
        return

    # Get RNA positions
    logging.info("Getting RNA positions...")
    rna1_info = {
        "gene": GeneInfo.find_gene_in_gff,
        "IGR": GeneInfo.find_igr_in_gff,
        "5UTR": lambda f,n: GeneInfo.find_utr_in_gff(f,n,"5UTR",UTR_LENGTH),
        "3UTR": lambda f,n: GeneInfo.find_utr_in_gff(f,n,"3UTR",UTR_LENGTH)
    }[RNA1["type"]](GFF_FILE, RNA1["name"])
    
    rna2_info = {
        "gene": GeneInfo.find_gene_in_gff,
        "IGR": GeneInfo.find_igr_in_gff,
        "5UTR": lambda f,n: GeneInfo.find_utr_in_gff(f,n,"5UTR",UTR_LENGTH),
        "3UTR": lambda f,n: GeneInfo.find_utr_in_gff(f,n,"3UTR",UTR_LENGTH)
    }[RNA2["type"]](GFF_FILE, RNA2["name"])

    if not all([rna1_info, rna2_info]):
        logging.error("Failed to get RNA positions!")
        return

    # Build indexes in parallel
    logging.info(f"Building indexes for {len(sample_info_list)} samples...")
    with ProcessPoolExecutor(max_workers=SAMPLE_PROCESS_POOL) as executor:
        build_results = list(executor.map(build_index, sample_info_list))
    
    if not all(build_results):
        logging.error("Index building failed for some samples!")
        return

    # Process samples in parallel
    logging.info(f"Processing {len(sample_info_list)} samples with {SAMPLE_PROCESS_POOL} processes...")
    results = []
    with ProcessPoolExecutor(max_workers=SAMPLE_PROCESS_POOL) as executor:
        futures = {
            executor.submit(process_sample, sample, rna1_info, rna2_info): sample
            for sample in sample_info_list
        }
        
        for future in as_completed(futures):
            results.append(future.result())
            sample_name, success, count, status, _ = results[-1]
            logging.info(f"Processed {sample_name}: {status} ({count} reads)")

    # Generate report
    logging.info("\nProcessing Summary:")
    logging.info("="*80)
    logging.info(f"{'Sample':<25} | {'Status':<10} | {'Reads':<10} | {'Details'}")
    logging.info("-"*80)
    
    success_count = 0
    total_chimeric = 0
    
    for sample_name, success, count, status, output_file in results:
        if success:
            group = next(g for s,g in SAMPLE_GROUPS if s == sample_name)
            group_files[group].append(output_file)
            success_count += 1
            total_chimeric += count
        logging.info(f"{sample_name:<25} | {'Success' if success else 'Failed':<10} | {count:<10} | {status}")
    
    logging.info("-"*80)
    logging.info(f"Total samples: {len(sample_info_list)}")
    logging.info(f"Successfully processed: {success_count}")
    logging.info(f"Failed: {len(sample_info_list) - success_count}")
    logging.info(f"Total chimeric sequences: {total_chimeric}")
    logging.info("="*80)

    # Merge group files
    logging.info("\nMerging group files...")
    for group, files in group_files.items():
        merged_file = os.path.join(OUTPUT_DIR, f"{group}_merged.fasta")
        with open(merged_file, 'w') as out:
            for f in files:
                with open(f) as infile:
                    out.write(infile.read())
        logging.info(f"Merged {len(files)} files -> {merged_file}")

    logging.info("Processing completed!")

    print(f"spend {time.time() - ALL_start_time}")

if __name__ == '__main__':
    main()
