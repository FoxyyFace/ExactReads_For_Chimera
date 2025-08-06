import pysam
import argparse
import os
import re
import gzip
from collections import defaultdict

def extract_read_ids(fasta_path):
    """从FASTA文件中提取所有read ID"""
    read_ids = set()
    
    open_func = gzip.open if fasta_path.endswith('.gz') else open
    mode = 'rt' if fasta_path.endswith('.gz') else 'r'
    
    with open_func(fasta_path, mode) as f:
        for line in f:
            if line.startswith('>'):
                # 提取 > 后的第一个非空白字符串作为read ID
                read_id = re.split(r'\s+', line[1:])[0]
                read_ids.add(read_id)
    return read_ids

def extract_segments_for_reads(bam_path, read_ids, output_path):
    """提取指定read IDs的所有比对segment"""
    segments_by_read = defaultdict(list)
    segment_counter = defaultdict(int)
    
    with pysam.AlignmentFile(bam_path, "rb", threads= 8) as bam:
        for read in bam:
            if read.query_name in read_ids and not read.is_unmapped:
                chrom = read.reference_name
                start = read.reference_start + 1
                end = read.reference_end
                strand = '-' if read.is_reverse else '+'
                
                segment_info = {
                    'segment_id': f"{chrom}:{start}-{end}({strand})",
                    'sequence': read.query_sequence,
                    'cigar': read.cigarstring
                }
                
                segments_by_read[read.query_name].append(segment_info)
                segment_counter[read.query_name] += 1
    
    # 写入输出文件
    with open(output_path, 'w') as out:
        for read_id, segments in segments_by_read.items():
            for i, seg in enumerate(segments):
                header = f">{read_id}|segment_{i+1}|{seg['segment_id']}"
                out.write(f"{header}\n{seg['sequence']}\n")
    
    # 打印统计信息
    total_reads = len(segments_by_read)
    total_segments = sum(segment_counter.values())
    print(f"提取完成: 共处理 {total_reads} 个reads, 提取 {total_segments} 个segments")
    print(f"平均每个read的segment数量: {total_segments/total_reads:.2f}")
    print(f"结果已保存至: {output_path}")

def main():
    parser = argparse.ArgumentParser(description='从BAM中提取嵌合reads的所有segment序列')
    parser.add_argument('-f', '--fasta', required=True, help='输入的FASTA文件路径 (嵌合reads提取结果)')
    parser.add_argument('-b', '--bam', required=True, help='样本对应的BAM文件路径')
    parser.add_argument('-o', '--output', default='extracted_segments.fa', help='输出文件名')
    args = parser.parse_args()

    print(f"从 {args.fasta} 提取read IDs...")
    read_ids = extract_read_ids(args.fasta)
    print(f"共找到 {len(read_ids)} 个reads")
    
    print(f"从 {args.bam} 提取segment信息...")
    extract_segments_for_reads(args.bam, read_ids, args.output)

if __name__ == '__main__':
    main()