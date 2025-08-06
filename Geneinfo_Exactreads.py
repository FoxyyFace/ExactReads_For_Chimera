import os
import subprocess
import pysam
import logging
import uuid
import struct
import mmap
import threading
import time
import shutil
import queue
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

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

class ReadIndex:
    """高效二进制索引系统，用于快速查询read比对信息"""
    """
    # 这里构建的是组合索引，https://developer.aliyun.com/article/841106
    # 这里的Struct和C语言的类似，用于定义数据类型，大小和字节序，和C语言不同的是，它是用于操作二进制数据的工具，而非创造一个数据类型变量
    # python 通过struct.pack生成一个bytes对象，unpack解包后得到元组tuple，通过索引访问
    # 是用于处理特定格式的二进制文件
    """
    """  
    # 索引文件头结构
    # 定义文件开头（head）数据块的格式  # 4SBIIQ请参考struct的格式字符定义，说明了head的C 数据块格式
    # < 代表小端序 (Little-endian)， 表示后面这些字节数据采用小端字节序
    # 小端代表低位字节存放在小地址，符合计算机读取内存的方式（效率更高），但是和人类阅读习惯相反
    # https://www.cnblogs.com/gremount/p/8830707.html 这篇文章很形象地说明了
    # https://docs.python.org/zh-cn/3/library/struct.html 这是python对于Struct的中文文档
    # 用C语言结构体表示如下：
    # struct {
    #     char magic[4];     // 'R' 'I' 'D' 'X'
    #     uint8_t version;   // 版本号 (1)
    #     uint32_t chrom_count; // 染色体数量
    #     uint32_t read_count;  // 总reads数
    #     uint64_t chrom_offset; // 染色体表偏移量
    # } header;
    """
    
    # 文件格式版本
    INDEX_VERSION = 1
    
    # 索引文件头结构
    HEADER_FORMAT = struct.Struct('<4sBIIQ')
    HEADER_MAGIC = b'RIDX'  # 索引文件标识
    
    # 染色体映射结构
    CHROM_MAP_FORMAT = struct.Struct('<I')
    
    # Read条目结构
    READ_HEADER_FORMAT = struct.Struct('<H')  # read ID长度 (2字节)
    SEGMENT_FORMAT = struct.Struct('<I B II')   # chrom_id, strand, start, end
    
    @staticmethod
    def build_index(bam_file, index_file, num_threads=4, bam_threads=4, chunk_size=50000):
        """
        构建二进制索引文件（多线程优化版）
        
        Args:
            bam_file: 输入BAM文件路径
            index_file: 输出索引文件路径
            num_threads: 使用的线程数
            bam_threads: BAM文件读取线程数
            chunk_size: 每个任务处理的read数量
            
        Returns:
            成功时返回索引文件路径，失败时抛出RuntimeError
        """
        logging.info(f"Building read index for {bam_file} with {num_threads} threads")
        start_time = time.time()
        
        # 创建染色体映射
        chrom_map = {}
        chrom_reverse_map = {}
        
        # 从BAM头部获取染色体信息
        with pysam.AlignmentFile(bam_file, "rb", threads=bam_threads) as bam:
            for i, chrom in enumerate(bam.references):
                chrom_map[chrom] = i
                chrom_reverse_map[i] = chrom
        
        # 创建临时目录用于存储中间结果
        temp_dir = f"{index_file}.tmp"
        os.makedirs(temp_dir, exist_ok=True)
        
        # 第一步：主线程读取BAM文件并生成任务队列
        task_queue = queue.Queue(maxsize=num_threads * 8)
        error_flag = threading.Event()  # 错误标志
        lock = threading.Lock()  # 用于保护共享变量
        chunk_files = []
        
        # 启动生产者线程
        def producer():
            try:
                with pysam.AlignmentFile(bam_file, "rb", threads=bam_threads) as bam:
                    current_chunk = []
                    current_count = 0
                    
                    for seg in bam:
                        if seg.is_unmapped or not seg.query_name:
                            continue
                        
                        # 收集当前read的所有segment
                        if not current_chunk or seg.query_name != current_chunk[-1][0]: 
                            current_chunk.append((seg.query_name, []))
                            current_count += 1
                        
                        chrom_id = chrom_map.get(seg.reference_name, 0xFFFFFFFF)
                        # 添加seg信息
                        current_chunk[-1][1].append((
                            chrom_id,
                            0 if not seg.is_reverse else 1,
                            seg.reference_start + 1,
                            seg.reference_end
                        ))
                        
                        # 当达到处理阈值时提交任务
                        if current_count >= chunk_size:
                            task_queue.put(current_chunk)
                            current_chunk = []
                            current_count = 0
                    
                    # 处理最后一批
                    if current_chunk:
                        task_queue.put(current_chunk)
            except Exception as e:
                logging.error(f"Producer thread failed: {str(e)}")
                error_flag.set()  # 设置错误标志
            finally:
                # 放入哨兵值，每个消费者线程一个
                for _ in range(num_threads):
                    task_queue.put(None)
        
        producer_thread = threading.Thread(target=producer)
        producer_thread.start()
        
        # 第二步：消费者线程处理任务队列
        chunk_index = 0

        def consumer():
            nonlocal chunk_index
            while not error_flag.is_set():
                try:
                    # 使用超时避免永久阻塞
                    current_chunk = task_queue.get(timeout=1.0)
                    
                    # 检查哨兵值
                    if current_chunk is None:
                        break

                    # 处理数据块
                    with lock:
                        chunk_file = os.path.join(temp_dir, f"chunk_{chunk_index}.dat")
                        chunk_index += 1
                    
                    ReadIndex._process_chunk(current_chunk, chrom_map, chunk_file)
                    
                    with lock:
                        chunk_files.append(chunk_file)
                except queue.Empty:
                    # 检查错误标志
                    if error_flag.is_set():
                        break
                    # 继续等待
                    continue
                except Exception as e:
                    logging.error(f"Consumer error: {str(e)}")
                    error_flag.set()  # 设置错误标志
                    break
        
        # 创建消费者线程池
        consumer_threads = []
        for _ in range(num_threads):
            t = threading.Thread(target=consumer)
            t.start()
            consumer_threads.append(t)
        
        # 等待生产者完成
        producer_thread.join()
        
        # 等待消费者完成
        for t in consumer_threads:
            t.join()
        
        # 检查错误
        if error_flag.is_set():
            # 清理临时文件
            try:
                shutil.rmtree(temp_dir)
            except Exception as e:
                logging.error(f"Failed to clean temp directory: {str(e)}")
            
            # 抛出异常
            raise RuntimeError("Error occurred during index building")
        
        # 第三步：合并所有临时分片到最终索引文件
        logging.info(f"Merging {len(chunk_files)} chunks into final index")
        total_reads = 0
        try:
            with open(index_file, "wb") as idx_file:
                # 写入文件头（占位）
                header = ReadIndex.HEADER_FORMAT.pack(
                    ReadIndex.HEADER_MAGIC,
                    ReadIndex.INDEX_VERSION,
                    len(chrom_map),
                    0,  # 占位符 (后续更新read计数)
                    0   # 占位符 (后续更新文件偏移)
                )
                idx_file.write(header)
                
                # 写入染色体映射表
                chrom_offset = idx_file.tell()
                for chrom, id in chrom_map.items():
                    chrom_bytes = chrom.encode('utf-8')
                    idx_file.write(ReadIndex.CHROM_MAP_FORMAT.pack(len(chrom_bytes)))
                    idx_file.write(chrom_bytes)
                
                # 合并所有分片
                data_offset = idx_file.tell()
                
                # 写入数据块头：总read数
                idx_file.write(struct.pack('Q', 0))  # 占位，稍后更新
                
                # 合并所有chunk文件
                for chunk_file in chunk_files:
                    with open(chunk_file, "rb") as cf:
                        # 读取chunk文件头（read数量）
                        chunk_reads = struct.unpack('Q', cf.read(8))[0]
                        total_reads += chunk_reads
                        # 复制read数据
                        shutil.copyfileobj(cf, idx_file)
                
                # 更新数据块头的总read数
                end_pos = idx_file.tell()
                idx_file.seek(data_offset)
                idx_file.write(struct.pack('Q', total_reads))
                idx_file.seek(end_pos)
                
                # 更新文件头
                idx_file.seek(0)
                header = ReadIndex.HEADER_FORMAT.pack(
                    ReadIndex.HEADER_MAGIC,
                    ReadIndex.INDEX_VERSION,
                    len(chrom_map),
                    total_reads,
                    chrom_offset
                )
                idx_file.write(header)
        except Exception as e:
            # 清理临时文件
            try:
                shutil.rmtree(temp_dir)
            except Exception as e:
                logging.error(f"Failed to clean temp directory: {str(e)}")
            
            # 抛出异常
            raise RuntimeError(f"Failed to merge index files: {str(e)}")
        
        # 清理临时文件
        try:
            shutil.rmtree(temp_dir)
        except Exception as e:
            logging.error(f"Failed to clean temp directory: {str(e)}")
        
        elapsed = time.time() - start_time
        logging.info(f"Built index with {total_reads} reads in {elapsed:.1f} seconds")
        return index_file
    
    @staticmethod
    def _process_chunk(chunk_data, chrom_map, output_file):
        """处理一个数据块并写入临时文件"""
        try:
            with open(output_file, "wb") as f:
                # 写入块头：块中的read数量
                f.write(struct.pack('Q', len(chunk_data)))
                
                for read_id, segments in chunk_data:
                    # 写入read ID
                    id_bytes = read_id.encode('utf-8')
                    f.write(ReadIndex.READ_HEADER_FORMAT.pack(len(id_bytes)))
                    f.write(id_bytes)
                    
                    # 写入segment数量
                    f.write(struct.pack('B', len(segments)))
                    
                    # 写入每个segment
                    for seg in segments:
                        f.write(ReadIndex.SEGMENT_FORMAT.pack(*seg))
        except Exception as e:
            # 删除可能损坏的文件
            if os.path.exists(output_file):
                try:
                    os.remove(output_file)
                except:
                    pass
            raise
    
    @staticmethod
    def load_index(index_file):
        """加载索引文件，返回内存映射对象和元数据"""
        if not os.path.exists(index_file):
            return None
        
        # 打开文件并内存映射
        with open(index_file, "rb") as f:
            mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        
        # 解析文件头
        magic, version, chrom_count, read_count, chrom_offset = \
            ReadIndex.HEADER_FORMAT.unpack_from(mm, 0)
        
        if magic != ReadIndex.HEADER_MAGIC:
            raise ValueError("Invalid index file format")
        if version != ReadIndex.INDEX_VERSION:
            raise ValueError(f"Unsupported index version: {version}")
        
        # 解析染色体映射
        chrom_map = {}
        chrom_reverse_map = {}
        pos = chrom_offset
        
        for _ in range(chrom_count):
            name_len = ReadIndex.CHROM_MAP_FORMAT.unpack_from(mm, pos)[0]
            pos += ReadIndex.CHROM_MAP_FORMAT.size
            
            chrom_name = mm[pos:pos+name_len].decode('utf-8')
            pos += name_len
            
            chrom_id = len(chrom_map)
            chrom_map[chrom_name] = chrom_id
            chrom_reverse_map[chrom_id] = chrom_name
        
        # 数据块起始位置
        data_offset = pos
        # 读取总read数
        total_reads = struct.unpack_from('Q', mm, data_offset)[0]
        data_offset += 8  # 跳过8字节的read计数头
        
        return {
            'mmap': mm,
            'read_count': total_reads,
            'chrom_map': chrom_map,
            'chrom_reverse_map': chrom_reverse_map,
            'data_offset': data_offset
        }
    
    @staticmethod
    def find_in_index(index_data, rna1_info, rna2_info, num_threads=4):
        """
        在索引中查找嵌合read
        
        Args:
            index_data: load_index返回的索引数据
            rna1_info: 第一个RNA区域信息
            rna2_info: 第二个RNA区域信息
            num_threads: 使用的线程数
            
        Returns:
            嵌合read ID集合
        """
        if not index_data:
            return set()
        
        mm = index_data['mmap']
        chrom_map = index_data['chrom_map']
        data_offset = index_data['data_offset']
        read_count = index_data['read_count']
        
        # 转换RNA信息为数字ID
        rna1_chrom = chrom_map.get(rna1_info["chrid"], 0xFFFFFFFF)
        rna2_chrom = chrom_map.get(rna2_info["chrid"], 0xFFFFFFFF)
        rna1_strand = rna1_info["strand"]
        rna2_strand = rna2_info["strand"]
        
        # 准备区域边界检查
        rna1_start, rna1_end = rna1_info["start"], rna1_info["end"]
        rna2_start, rna2_end = rna2_info["start"], rna2_info["end"]
        
        # 计算每个线程处理的范围
        chunk_size = max(1, read_count // num_threads)
        chunks = []
        pos = data_offset
        
        for i in range(num_threads):
            start = pos
            # 找到当前块的结束位置
            if i == num_threads - 1:
                end = len(mm)  # 最后一个块到文件末尾
            else:
                # 跳过指定数量的read
                for _ in range(chunk_size):
                    # 检查是否超出文件范围
                    if pos + ReadIndex.READ_HEADER_FORMAT.size > len(mm):
                        break
                    
                    # 读取read ID长度
                    id_len = ReadIndex.READ_HEADER_FORMAT.unpack_from(mm, pos)[0]
                    pos += ReadIndex.READ_HEADER_FORMAT.size
                    
                    # 检查是否超出文件范围
                    if pos + id_len > len(mm):
                        break
                    
                    # 跳过read ID
                    pos += id_len
                    
                    # 检查是否超出文件范围
                    if pos + 1 > len(mm):
                        break
                    
                    # 读取segment数量
                    seg_count = mm[pos]
                    pos += 1
                    
                    # 检查是否超出文件范围
                    if pos + seg_count * ReadIndex.SEGMENT_FORMAT.size > len(mm):
                        break
                    
                    # 跳过所有segment
                    pos += seg_count * ReadIndex.SEGMENT_FORMAT.size
                
                end = pos
            chunks.append((start, end))
        
        # 并行处理每个块
        chimeric_ids = set()
        lock = threading.Lock()
        
        def process_chunk(start, end):
            local_ids = set()
            pos = start
            
            while pos < end:
                # 检查是否超出文件范围
                if pos + ReadIndex.READ_HEADER_FORMAT.size > end:
                    break
                
                # 读取read ID长度
                id_len = ReadIndex.READ_HEADER_FORMAT.unpack_from(mm, pos)[0]
                pos += ReadIndex.READ_HEADER_FORMAT.size
                
                # 检查是否超出文件范围
                if pos + id_len > end:
                    break
                
                # 读取read ID
                try:
                    read_id = mm[pos:pos+id_len].decode('utf-8')
                except UnicodeDecodeError:
                    # 使用错误处理机制替换无效字节
                    read_id = mm[pos:pos+id_len].decode('utf-8', errors='replace')
                    logging.warning(f"Invalid UTF-8 sequence in read ID at position {pos}, using replacement characters")
                
                pos += id_len
                
                # 检查是否超出文件范围
                if pos + 1 > end:
                    break
                
                # 读取segment数量
                seg_count = mm[pos]
                pos += 1
                
                # 检查是否超出文件范围
                if pos + seg_count * ReadIndex.SEGMENT_FORMAT.size > end:
                    break
                
                # 读取所有segment
                segments = []
                for _ in range(seg_count):
                    seg_data = ReadIndex.SEGMENT_FORMAT.unpack_from(mm, pos)
                    pos += ReadIndex.SEGMENT_FORMAT.size
                    segments.append(seg_data)
                
                # 检查是否为嵌合read
                if ReadIndex._is_chimera(segments, rna1_chrom, rna1_strand, 
                                        rna1_start, rna1_end,
                                        rna2_chrom, rna2_strand,
                                        rna2_start, rna2_end):
                    local_ids.add(read_id)
            
            with lock:
                chimeric_ids.update(local_ids)
        
        # 使用线程池处理
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            futures = []
            for start, end in chunks:
                futures.append(executor.submit(process_chunk, start, end))
            
            for future in futures:
                future.result()
        
        return chimeric_ids
    
    @staticmethod
    def _is_chimera(segments, rna1_chrom, rna1_strand, rna1_start, rna1_end,
                   rna2_chrom, rna2_strand, rna2_start, rna2_end):
        """快速检查嵌合模式，使用预计算的数字ID"""
        for i in range(len(segments) - 1):
            seg1 = segments[i]
            seg2 = segments[i+1]
            
            # 检查染色体和链方向
            if not (seg1[0] == rna1_chrom and seg1[1] == rna1_strand):
                continue
            
            # 检查位置重叠
            if not (seg1[2] <= rna1_end and seg1[3] >= rna1_start):
                continue
            
            # 检查第二个segment
            if (seg2[0] == rna2_chrom and 
                seg2[1] == rna2_strand and
                seg2[2] <= rna2_end and 
                seg2[3] >= rna2_start):
                return True
        
        return False

class ExactReads:
    """Utility class for identifying and extracting chimeric reads"""
    
    # 索引缓存
    INDEX_CACHE = {}
    INDEX_LOCK = threading.Lock()
    
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
    def find_chimeric_read_ids(bam_file, rna1_info, rna2_info, use_index=True, 
                              num_threads=4, bam_threads=4, chunk_size=50000):
        """
        Identifies chimeric RNA1->RNA2 reads in BAM file
        
        Args:
            bam_file: Path to BAM file
            rna1_info: Position info for first RNA region
            rna2_info: Position info for second RNA region
            use_index: 是否使用索引加速
            num_threads: 使用的线程数
            bam_threads: BAM文件读取线程数
            chunk_size: 每个任务处理的read数量
            
        Returns:
            Set of chimeric read IDs
        """
        if use_index:
            return ExactReads._find_with_index(bam_file, rna1_info, rna2_info, 
                                             num_threads, bam_threads, chunk_size)
        else:
            return ExactReads._find_without_index(bam_file, rna1_info, rna2_info, 
                                                 bam_threads)
    
    @staticmethod
    def _find_with_index(bam_file, rna1_info, rna2_info, num_threads, bam_threads, chunk_size):
        """使用索引加速查找"""
        # 获取或创建索引
        index_file = f"{bam_file}.ridx"
        index_data = None
        
        with ExactReads.INDEX_LOCK:
            if index_file in ExactReads.INDEX_CACHE:
                index_data = ExactReads.INDEX_CACHE[index_file]
            else:
                # 检查索引文件是否存在
                if not os.path.exists(index_file):
                    logging.info(f"Building index for {bam_file}")
                    ReadIndex.build_index(bam_file, index_file, 
                                        num_threads=num_threads, 
                                        bam_threads=bam_threads,
                                        chunk_size=chunk_size)
                
                # 加载索引
                index_data = ReadIndex.load_index(index_file)
                if index_data:
                    ExactReads.INDEX_CACHE[index_file] = index_data
        
        if not index_data:
            logging.warning("Index not available, falling back to standard method")
            return ExactReads._find_without_index(bam_file, rna1_info, rna2_info, 
                                                 bam_threads)
        
        # 使用索引查找
        logging.info(f"Searching for chimeric reads using index")
        chimeric_ids = ReadIndex.find_in_index(index_data, rna1_info, rna2_info, num_threads)
        logging.info(f"Found {len(chimeric_ids)} chimeric reads")
        return chimeric_ids
    
    @staticmethod
    def _find_without_index(bam_file, rna1_info, rna2_info, bam_threads):
        """不使用索引的标准查找方法（使用多线程读取BAM）"""
        chimeric_ids = set()
        current_read_segments = []
        current_read_id = None

        # 使用多线程读取BAM文件
        with pysam.AlignmentFile(bam_file, "rb", threads=bam_threads) as bam:
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
