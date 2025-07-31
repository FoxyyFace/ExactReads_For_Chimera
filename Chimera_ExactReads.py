import os
import subprocess
import pysam
import logging
import uuid

class GeneInfo:
    """Utility class for extracting gene information from GFF files"""
    
    @staticmethod
    def find_gene_in_gff(gff_file, gene_name):
        """
        Finds gene position information in GFF file
        
        Args:
            gff_file: Path to GFF file
            gene_name: Gene name to search for (e.g. "PA14_RS00005")
            
        Returns:
            dict with chromosome, start/end positions, and strand information
        """
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
        
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue

                attr_dict = dict(attr.split('=', 1) for attr in parts[8].split(';') if '=' in attr)

                if ('locus_tag' in attr_dict and attr_dict['locus_tag'] == gene_name) or \
                   ('locus' in attr_dict and attr_dict['locus'] == gene_name):
                    
                    strand_flag = 0 if parts[6] == '+' else 1
                    return {
                        "chrid": parts[0],
                        "start": int(parts[3]),
                        "end": int(parts[4]),
                        "strand": strand_flag
                    }
                
        logging.error(f"Gene '{gene_name}' not found in GFF file")
        return None

    @staticmethod
    def find_igr_in_gff(gff_file, igr_name):
        """
        Finds position information for an Intergenic Region (IGR)
        
        Args:
            gff_file: Path to GFF file
            igr_name: IGR name in "geneA:geneB" format
            
        Returns:
            dict with IGR chromosome, start/end positions, and strand information
        """
        if ":" not in igr_name:
            logging.error(f"Invalid IGR name format '{igr_name}', should be 'geneA:geneB'")
            return None
        
        gene1, gene2 = igr_name.split(":", 1)
        gene1_info = GeneInfo.find_gene_in_gff(gff_file, gene1)
        gene2_info = GeneInfo.find_gene_in_gff(gff_file, gene2)
        
        if not gene1_info or not gene2_info:
            return None
        
        if gene1_info["chrid"] != gene2_info["chrid"]:
            logging.error(f"Genes '{gene1}' and '{gene2}' on different chromosomes")
            return None
        
        if gene1_info["strand"] != gene2_info["strand"]:
            logging.warning(f"Genes '{gene1}' and '{gene2}' on different strands")
        
        # Determine IGR boundaries
        if gene1_info["end"] < gene2_info["start"]:
            start = gene1_info["end"] + 1
            end = gene2_info["start"] - 1
            strand = gene1_info["strand"]
        else:
            start = gene2_info["end"] + 1
            end = gene1_info["start"] - 1
            strand = gene2_info["strand"]
        
        if start > end:
            logging.error(f"Invalid IGR position {start}-{end}")
            return None
        
        return {"chrid": gene1_info["chrid"], "start": start, "end": end, "strand": strand}

    @staticmethod
    def find_utr_in_gff(gff_file, gene_name, utr_type, utr_length=150):
        """
        Finds UTR position information for a gene
        
        Args:
            gff_file: Path to GFF file
            gene_name: Target gene name
            utr_type: "5UTR" or "3UTR"
            utr_length: Length of UTR region in bp
            
        Returns:
            dict with UTR chromosome, start/end positions, and strand information
        """
        gene_info = GeneInfo.find_gene_in_gff(gff_file, gene_name)
        if not gene_info:
            return None
            
        # Calculate UTR coordinates based on strand orientation
        if utr_type == "5UTR":
            if gene_info["strand"] == 0:  # Forward strand
                start = max(1, gene_info["start"] - utr_length)
                end = gene_info["start"] - 1
            else:  # Reverse strand
                start = gene_info["end"] + 1
                end = gene_info["end"] + utr_length
        elif utr_type == "3UTR":
            if gene_info["strand"] == 0:  # Forward strand
                start = gene_info["end"] + 1
                end = gene_info["end"] + utr_length
            else:  # Reverse strand
                start = max(1, gene_info["start"] - utr_length)
                end = gene_info["start"] - 1
        else:
            logging.error(f"Invalid UTR type '{utr_type}', must be '5UTR' or '3UTR'")
            return None
        
        if start > end:
            logging.error(f"Invalid UTR position {start}-{end}")
            return None
        
        return {
            "chrid": gene_info["chrid"],
            "start": start,
            "end": end,
            "strand": gene_info["strand"]
        }

class ExactReads:
    """Utility class for identifying and extracting chimeric reads"""

    @staticmethod
    def is_chimera(segments, rna1_info, rna2_info):
        """
        Determines if a read is a chimeric RNA1-RNA2 fragment
        
        Args:
            segments: Alignment segments of a read
            rna1_info: Position info for first RNA region
            rna2_info: Position info for second RNA region
            
        Returns:
            True if segments show RNA1->RNA2 chimera, False otherwise
        """
        for i in range(len(segments) - 1):
            seg1, seg2 = segments[i], segments[i+1]
            
            # Check RNA1 overlap
            if not (seg1["chrid"] == rna1_info["chrid"] and
                    seg1["strand"] == rna1_info["strand"] and
                    seg1["start"] <= rna1_info["end"] and
                    seg1["end"] >= rna1_info["start"]):
                continue
            
            # Check RNA2 overlap
            if (seg2["chrid"] == rna2_info["chrid"] and
                seg2["strand"] == rna2_info["strand"] and
                seg2["start"] <= rna2_info["end"] and
                seg2["end"] >= rna2_info["start"]):
                return True
        return False

    @staticmethod
    def find_chimeric_read_ids(bam_file, rna1_info, rna2_info):
        """
        Identifies chimeric RNA1->RNA2 reads in BAM file
        
        Args:
            bam_file: Path to BAM file
            rna1_info: Position info for first RNA region
            rna2_info: Position info for second RNA region
            
        Returns:
            Set of chimeric read IDs
        """
        chimeric_ids = set()
        current_read_segments = []
        current_read_id = None

        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for seg in bam:
                if seg.is_unmapped or not seg.query_name:
                    continue

                if seg.query_name != current_read_id:
                    if current_read_id and len(current_read_segments) > 1:
                        if ExactReads.is_chimera(current_read_segments, rna1_info, rna2_info):
                            chimeric_ids.add(current_read_id)
                    current_read_id = seg.query_name
                    current_read_segments = []

                current_read_segments.append({
                    "chrid": seg.reference_name,
                    "start": seg.reference_start + 1,
                    "end": seg.reference_end,
                    "strand": 0 if not seg.is_reverse else 1
                })

            if current_read_id and len(current_read_segments) > 1:
                if ExactReads.is_chimera(current_read_segments, rna1_info, rna2_info):
                    chimeric_ids.add(current_read_id)

        logging.info(f"Found {len(chimeric_ids)} chimeric reads")
        return chimeric_ids

    @staticmethod
    def extract_reads_with_seqtk(fastq_files, read_ids, output_file, seqtk_path, output_format="fasta"):
        """
        Extracts specific reads using seqtk
        
        Args:
            fastq_files: FASTQ file(s) (single or pair)
            read_ids: Collection of read IDs to extract
            output_file: Destination file path
            seqtk_path: Path to seqtk executable
            output_format: "fasta" or "fastq"
            
        Returns:
            True on success, False on failure
        """
        if not read_ids:
            logging.warning("No read IDs to extract")
            return False

        temp_id_file = f"temp_ids_{uuid.uuid4().hex}.txt"
        with open(temp_id_file, "w") as f:
            f.write("\n".join(read_ids) + "\n")
        
        logging.info(f"Extracting {len(read_ids)} reads using seqtk")

        try:
            if isinstance(fastq_files, tuple):  # Paired-end
                temp1 = f"temp_r1_{uuid.uuid4().hex}.fq"
                temp2 = f"temp_r2_{uuid.uuid4().hex}.fq"
                
                # Process R1
                cmd = f"{seqtk_path} subseq {fastq_files[0]} {temp_id_file} > {temp1}"
                if subprocess.run(cmd, shell=True, timeout=7200).returncode != 0:
                    return False
                
                # Process R2
                cmd = f"{seqtk_path} subseq {fastq_files[1]} {temp_id_file} > {temp2}"
                if subprocess.run(cmd, shell=True, timeout=7200).returncode != 0:
                    return False
                
                # Combine and convert if needed
                combine_cmd = f"cat {temp1} {temp2}"
                if output_format == "fasta":
                    combine_cmd += f" | {seqtk_path} seq -A"
                
                with open(output_file, "w") as f:
                    subprocess.run(combine_cmd, shell=True, stdout=f)
                
                # Cleanup
                os.remove(temp1)
                os.remove(temp2)
            else:  # Single-end
                cmd = f"{seqtk_path} subseq {fastq_files} {temp_id_file}"
                if output_format == "fasta":
                    cmd += f" | {seqtk_path} seq -A"
                
                with open(output_file, "w") as f:
                    result = subprocess.run(cmd, shell=True, stdout=f, timeout=7200)
                    if result.returncode != 0:
                        return False
        
        except Exception as e:
            logging.error(f"seqtk extraction failed: {str(e)}")
            return False
        finally:
            if os.path.exists(temp_id_file):
                os.remove(temp_id_file)
        
        return True