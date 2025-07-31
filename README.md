# Chimera Exact Reads Tool

## Overview
This tool identifies and extracts chimeric RNA reads that span two specific genomic regions (RNA1 and RNA2) from RNA-seq data. It processes BAM files to identify chimeric reads and extracts their sequences from FASTQ files using seqtk.

## Key Features
- Identifies chimeric reads spanning two specified genomic regions
- Supports gene, IGR, 5'UTR, and 3'UTR regions
- Parallel processing for efficient analysis
- Generates both sample-specific and group-merged output files
- Comprehensive logging and progress tracking

## Dependencies
### Python Packages
| Package | Version |
|---------|---------|
| Python  | 3.13.5  |
| pysam   | 0.23.3  |
| tqdm    | 4.67.1  |

### External Tools
- **seqtk**: Required for sequence extraction
  - Download: https://github.com/lh3/seqtk
  - Must be installed and path specified in configuration

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/FoxyyFace/ExactReads_For_Chimeras.git
   ```

2. Install required Python packages:
   ```bash
   pip install pysam==0.23.3 tqdm==4.67.1
   ```

3. Install seqtk:
   ```bash
   git clone https://github.com/lh3/seqtk.git
   cd seqtk
   make
   ```

4. Verify installation:
   ```bash
   ./seqtk
   ```

## Configuration
Modify the configuration section in `ExactReads_For_Chimera.py`:

```python
# Configuration parameters
seqtk_path = "/path/to/seqtk"  # Path to seqtk executable
RNA1 = {
    "name" : "EFV83_RS15570:EFV83_RS15625",  # Target region 1
    "type" : "IGR",     # Options: gene, IGR, 5UTR, 3UTR
    "strand" : 0        # 0 = + strand, 1 = - strand
}
RNA2 = {
    "name" : "PA14_RS29125",    # Target region 2
    "type" : "gene",     # Options: gene, IGR, 5UTR, 3UTR
    "strand" : 0         # 0 = + strand, 1 = - strand
}
UTR_LENGTH = 150         # Length of UTR regions (bp)
SAMPLE_GROUP_NAME = {    # Sample names and group associations
    ("Sample1", "GroupA"),
    ("Sample2", "GroupA"),
    ("Sample3", "GroupB"),
}
FASTA_FILE = "/path/to/reference.fasta"  # Reference genome. It is temporarily useless, so you can leave it out
GFF_FILE = "/path/to/annotation.gff"     # Genome annotation
BAM_DIR = "/path/to/bam_files"           # Directory with BAM/FASTQ files
FASTQ_BAM_PREFIX = "trimmed_"            # File prefix
OUTPUT_DIR = "output"                    # Output directory
is_paired_end = True                     # Paired-end data
suffix_read1 = "_1"
suffix_read2 = "_2"
file_type = ".fq.gz"
MAX_WORKERS = 4                          # Number of parallel processes
OUTPUT_FORMAT = "fasta"                  # Output format: fasta or fastq
```

## Input Files
1. **Reference Genome**: FASTA format
2. **Genome Annotation**: GFF format
3. **Alignment Files**:
   - BAM files (sorted by read names)
   - Corresponding FASTQ files (gzipped or uncompressed)

## File Naming Convention
Files must follow this pattern:
```
{BAM_DIR}/{FASTQ_BAM_PREFIX}{sample_name}.bam
{BAM_DIR}/{FASTQ_BAM_PREFIX}{sample_name}_{suffix_read1}.{file_type}  # Read 1
{BAM_DIR}/{FASTQ_BAM_PREFIX}{sample_name}_{suffix_read2}.{file_type}  # Read 2 (if paired-end)
```

## Usage
Run the main script:
```bash
python ExactReads_For_Chimera.py
```

## Output Files
- `processing.log`: Detailed processing log
- Sample-specific outputs:
  - `output/{sample_name}_chimeric.fasta`: Chimeric reads for each sample
- Group-merged outputs:
  - `output/{group_name}_merged_chimeric.fasta`: Combined chimeric reads for each group

## Workflow
1. **Configuration**: Set parameters in the script
2. **Gene Position Identification**: 
   - Extract genomic coordinates for RNA1 and RNA2 from GFF
3. **Chimeric Read Identification**:
   - Scan BAM files for reads spanning both regions
4. **Sequence Extraction**:
   - Extract identified reads from FASTQ files using seqtk
5. **Output Generation**:
   - Create sample-specific FASTA/FASTQ files
   - Merge files by experimental group

## Troubleshooting
1. **Missing Files Warning**:
   - Verify all input files exist at specified paths
   - Check file naming conventions

2. **No Chimeric Reads Found**:
   - Verify RNA regions are correctly specified
   - Check strand orientation matches biological reality
   - Ensure sufficient sequencing depth
