import os
import concurrent.futures
import logging
from tqdm import tqdm
from Chimera_ExactReads import ExactReads, GeneInfo

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("processing.log"),
        logging.StreamHandler()
    ]
)

# Configuration parameters
seqtk_path = "/path/to/seqtk"
RNA1 = {
    "name" : "EFV83_RS15570:EFV83_RS15625",
    "type" : "IGR",     # gene, IGR , 5UTR, 3UTR
    "strand" : 0       # 0 is + , 1 is -
}
RNA2 = {
    "name" : "PA14_RS29125",
    "type" : "gene",     # gene, IGR , 5UTR, 3UTR
    "strand" : 0       # 0 is + , 1 is -
}
UTR_LENGTH = 150    # bp
SAMPLE_GROUP_NAME = {
    ("Mabin-1-LHF18430_L5", "MabcoPA14_IN"),
    ("Mabin-2-LHF18431_L5", "MabcoPA14_IN"),
    ("Mabin-3-LHF18432_L5", "MabcoPA14_IN"),
}
FASTA_FILE = "/path/to/reference.fasta"
GFF_FILE = "/path/to/annotation.gff"
BAM_DIR = "/path/to/bam_files"
FASTQ_BAM_PREFIX = "trimmed_"
OUTPUT_DIR = "output"
is_paired_end = True
suffix_read1 = "_1"
suffix_read2 = "_2"
file_type = ".fq.gz"
MAX_WORKERS = 3  # Adjust based on CPU cores

# Output format (fasta or fastq)
OUTPUT_FORMAT = "fasta" 

################
#              #
#   RUNNING    #
#              #
################

os.makedirs(OUTPUT_DIR, exist_ok=True)


sample_info_list = []
group_files = {} 

# Process all samples
for sample_name, group_name in SAMPLE_GROUP_NAME:
    # Build BAM file path
    bam_file = os.path.join(BAM_DIR, f"{FASTQ_BAM_PREFIX}{sample_name}.bam")
    
    # Build FASTQ file paths
    if is_paired_end:
        fastq_r1 = os.path.join(BAM_DIR, f"{FASTQ_BAM_PREFIX}{sample_name}{suffix_read1}{file_type}")
        fastq_r2 = os.path.join(BAM_DIR, f"{FASTQ_BAM_PREFIX}{sample_name}{suffix_read2}{file_type}")
        
        if not os.path.exists(fastq_r1) or not os.path.exists(fastq_r2):
            logging.warning(f"FASTQ files missing: {fastq_r1} or {fastq_r2}")
            continue
    else:
        fastq_file = os.path.join(BAM_DIR, f"{FASTQ_BAM_PREFIX}{sample_name}{file_type}")
        
        if not os.path.exists(fastq_file):
            logging.warning(f"FASTQ file missing: {fastq_file}")
            continue
    
    if not os.path.exists(bam_file):
        logging.warning(f"BAM file missing: {bam_file}")
        continue
    
    # Store sample info
    sample_info = {
        "sample_name": sample_name,
        "group_name": group_name,
        "bam_file": bam_file,
        "is_paired_end": is_paired_end
    }
    
    if is_paired_end:
        sample_info["fastq_r1"] = fastq_r1
        sample_info["fastq_r2"] = fastq_r2
    else:
        sample_info["fastq"] = fastq_file
    
    sample_info_list.append(sample_info)
    
    if group_name not in group_files:
        group_files[group_name] = []

# Get RNA1 and RNA2 position info
rna1_info = None
rna2_info = None

# Get position info based on RNA type
if RNA1["type"] == "gene":
    logging.info("Getting RNA1 position (gene)")
    rna1_info = GeneInfo.find_gene_in_gff(GFF_FILE, RNA1["name"])
elif RNA1["type"] == "IGR":
    logging.info("Getting RNA1 position (IGR)")
    rna1_info = GeneInfo.find_igr_in_gff(GFF_FILE, RNA1["name"])
elif RNA1["type"] in ["5UTR", "3UTR"]:
    logging.info(f"Getting RNA1 position ({RNA1['type']})")
    rna1_info = GeneInfo.find_utr_in_gff(GFF_FILE, RNA1["name"], RNA1["type"], UTR_LENGTH)

if RNA2["type"] == "gene":
    logging.info("Getting RNA2 position (gene)")
    rna2_info = GeneInfo.find_gene_in_gff(GFF_FILE, RNA2["name"])
elif RNA2["type"] == "IGR":
    logging.info("Getting RNA2 position (IGR)")
    rna2_info = GeneInfo.find_igr_in_gff(GFF_FILE, RNA2["name"])
elif RNA2["type"] in ["5UTR", "3UTR"]:
    logging.info(f"Getting RNA2 position ({RNA2['type']})")
    rna2_info = GeneInfo.find_utr_in_gff(GFF_FILE, RNA2["name"], RNA2["type"], UTR_LENGTH)

if not rna1_info or not rna2_info:
    logging.error("Failed to get RNA1 or RNA2 position info")
    exit(1)

# Function to process a single sample
def process_sample(sample_info):
    try:
        logging.info(f"\nProcessing sample: {sample_info['sample_name']} (Group: {sample_info['group_name']})")
        
        # Identify chimeric reads
        logging.info(f"Step 1: Scanning BAM for chimeric reads")
        chimeric_ids = ExactReads.find_chimeric_read_ids(
            sample_info["bam_file"], 
            rna1_info, 
            rna2_info
        )
        
        if not chimeric_ids:
            logging.warning(f"No chimeric reads found")
            return sample_info["sample_name"], False, 0, "No chimeric reads", None
        
        logging.info(f"Found {len(chimeric_ids)} chimeric reads")
        
        # Extract sequences
        logging.info(f"Step 2: Extracting chimeric sequences")
        
        if OUTPUT_FORMAT == "fasta":
            output_file = os.path.join(OUTPUT_DIR, f"{sample_info['sample_name']}_chimeric.fasta")
        else:
            output_file = os.path.join(OUTPUT_DIR, f"{sample_info['sample_name']}_chimeric.fastq")
        
        if sample_info["is_paired_end"]:
            success = ExactReads.extract_reads_with_seqtk(
                fastq_files=(sample_info["fastq_r1"], sample_info["fastq_r2"]),
                read_ids=chimeric_ids,
                output_file=output_file,
                seqtk_path=seqtk_path,
                output_format=OUTPUT_FORMAT
            )
        else:
            success = ExactReads.extract_reads_with_seqtk(
                fastq_files=sample_info["fastq"],
                read_ids=chimeric_ids,
                output_file=output_file,
                seqtk_path=seqtk_path,
                output_format=OUTPUT_FORMAT
            )
        
        if success:
            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                logging.info(f"Successfully extracted to {output_file} (Size: {os.path.getsize(output_file)} bytes)")
                group_files[sample_info["group_name"]].append(output_file)
                return sample_info["sample_name"], True, len(chimeric_ids), "Success", output_file
            else:
                logging.error(f"Empty output file: {output_file}")
                return sample_info["sample_name"], False, len(chimeric_ids), "Empty output", None
        else:
            logging.error(f"Failed to extract sequences")
            return sample_info["sample_name"], False, len(chimeric_ids), "Extraction failed", None
    
    except Exception as e:
        logging.error(f"Error processing sample {sample_info['sample_name']}: {str(e)}", exc_info=True)
        return sample_info["sample_name"], False, 0, f"Error: {str(e)}", None

# Function to merge group files
def merge_group_files(group_name, file_list):
    """Merge all sample files for a group"""
    if not file_list:
        logging.warning(f"No valid output files to merge for group {group_name}")
        return
    
    if OUTPUT_FORMAT == "fasta":
        merged_file = os.path.join(OUTPUT_DIR, f"{group_name}_merged_chimeric.fasta")
    else:
        merged_file = os.path.join(OUTPUT_DIR, f"{group_name}_merged_chimeric.fastq")
    
    logging.info(f"Merging {len(file_list)} files for group {group_name} into {merged_file}")
    
    try:
        with open(merged_file, "w") as out_f:
            for file_path in file_list:
                if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
                    with open(file_path, "r") as in_f:
                        out_f.write(in_f.read())
                    logging.debug(f"Merged file: {file_path}")
                else:
                    logging.warning(f"Skipping invalid file: {file_path}")
        
        if os.path.exists(merged_file) and os.path.getsize(merged_file) > 0:
            logging.info(f"Successfully merged group {group_name} (Size: {os.path.getsize(merged_file)} bytes)")
            return merged_file
        else:
            logging.error(f"Empty merged file: {merged_file}")
            return None
    
    except Exception as e:
        logging.error(f"Error merging group {group_name}: {str(e)}", exc_info=True)
        return None

# Process samples in parallel
total_samples = len(sample_info_list)

logging.info(f"Starting parallel processing of {total_samples} samples (using {MAX_WORKERS} threads)")
logging.info(f"Output format: {OUTPUT_FORMAT.upper()}")

# Initialize progress bar
pbar = tqdm(total=total_samples, desc="Processing samples", unit="sample", position=0)
results = []

# Use ThreadPoolExecutor
with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
    futures = []
    for sample_info in sample_info_list:
        future = executor.submit(process_sample, sample_info)
        futures.append(future)
    
    for i, future in enumerate(concurrent.futures.as_completed(futures)):
        try:
            result = future.result()
            results.append(result)
            sample_name, success, count, status, output_file = result
            
            pbar.update(1)
            pbar.set_description(f"Processing: {sample_name}")
            pbar.set_postfix_str(f"Status: {status}")
            
            logging.info(f"Progress: {i+1}/{total_samples} - {sample_name}: {status}")
        except Exception as e:
            logging.error(f"Task failed: {str(e)}", exc_info=True)
            pbar.update(1)
            pbar.set_postfix_str(f"Status: Error")
            results.append((None, False, 0, f"Error: {str(e)}", None))

pbar.close()

logging.info("\nAll samples processed! Summary:")
logging.info("-" * 80)
logging.info(f"{'Sample Name':<25} | {'Status':<15} | {'Chimeric Count':<10} | {'Details'}")
logging.info("-" * 80)

success_count = 0
total_chimeric = 0

for result in results:
    if result is None:
        continue
    sample_name, success, count, status, _ = result
    status_text = "Success" if success else "Failed"
    logging.info(f"{sample_name:<25} | {status_text:<15} | {count:<10} | {status}")
    if success:
        success_count += 1
    total_chimeric += count

logging.info("-" * 80)
logging.info(f"Total samples: {total_samples}")
logging.info(f"Successfully processed: {success_count}")
logging.info(f"Failed: {total_samples - success_count}")
logging.info(f"Total chimeric sequences: {total_chimeric}")
logging.info(f"Output format: {OUTPUT_FORMAT.upper()}")
logging.info("-" * 80)

# Merge group files
logging.info("\nMerging group files...")
merged_files = {}
for group_name, file_list in group_files.items():
    merged_file = merge_group_files(group_name, file_list)
    if merged_file:
        merged_files[group_name] = merged_file

if merged_files:
    logging.info("\nGroup files merged:")
    for group_name, merged_file in merged_files.items():
        logging.info(f"Group {group_name}: {merged_file} (Size: {os.path.getsize(merged_file)} bytes)")
else:
    logging.warning("No group files were merged")

# Check for empty output files
for sample_name, _, _, _, output_file in results:
    if output_file and os.path.exists(output_file):
        if os.path.getsize(output_file) == 0:
            logging.warning(f"Empty output file: {output_file}")
    elif output_file:
        logging.warning(f"Output file missing: {output_file}")

# Check for empty merged files
for group_name, merged_file in merged_files.items():
    if os.path.exists(merged_file):
        if os.path.getsize(merged_file) == 0:
            logging.warning(f"Empty merged file: {merged_file}")
    else:
        logging.warning(f"Merged file missing: {merged_file}")