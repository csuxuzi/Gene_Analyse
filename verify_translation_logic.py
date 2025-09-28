import vcfpy
import requests
from datetime import datetime
import json
import math

ENSEMBL_SERVER = "https://rest.ensembl.org"
JSON_HEADERS = {"Content-Type": "application/json"}

CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def translate_dna_std(dna):
    # 你自己的翻译函数，确保用核基因的标准遗传密码表，最后把末端终止符去掉
    # 这里给个极简示意，替换成你现有的 translate_dna 即可
    from Bio.Seq import Seq
    aa = str(Seq(dna).translate(table=1))  # 核基因: Standard code
    return aa.rstrip('*')

def reverse_complement(seq):
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]

def translate_dna_std(dna):
    # 用你现有的 translate_dna 也行；注意最后去掉末尾'*'
    from Bio.Seq import Seq
    return str(Seq(dna).translate(table=1)).rstrip('*')

def verify_translation_logic_fixed(transcript_id: str):
    server = "https://rest.ensembl.org"
    headers = {"Content-Type": "application/json"}

    print(f"\n===== 开始对转录本 {transcript_id} 进行双向翻译验证 (Fixed) =====")

    # 1) 官方蛋白
    r_prot = requests.get(f"{server}/sequence/id/{transcript_id}?type=protein", headers=headers)
    r_prot.raise_for_status()
    official_protein_seq = r_prot.json()["seq"]
    print(f"  > 官方蛋白长度: {len(official_protein_seq)} aa")

    # 2) 转录本信息（含外显子、Translation 的基因组坐标）
    info = requests.get(f"{server}/lookup/id/{transcript_id}?expand=1", headers=headers).json()
    strand = info["strand"]
    all_exons = info["Exon"]
    t = info["Translation"]
    cds_g_start, cds_g_end = t["start"], t["end"]  # 基因组坐标
    cds_lo, cds_hi = min(cds_g_start, cds_g_end), max(cds_g_start, cds_g_end)
    print(f"  > Transcript on strand {strand}; CDS genomic span: {cds_lo}-{cds_hi}")

    # 3) 选出与 CDS 有重叠的外显子
    coding_exons = [e for e in all_exons if max(e['start'], cds_lo) <= min(e['end'], cds_hi)]
    # 按转录本阅读方向排序
    if strand == 1:
        coding_exons.sort(key=lambda x: x['start'])
    else:
        coding_exons.sort(key=lambda x: x['start'], reverse=True)
    print(f"  > 编码外显子数: {len(coding_exons)}")

    # 4) 逐外显子取序列并按基因组重叠区切片
    cds_parts = []
    total_bp = 0
    for exon in coding_exons:
        e_start, e_end, e_strand = exon['start'], exon['end'], exon['strand']
        ovl_start = max(e_start, cds_lo)
        ovl_end   = min(e_end, cds_hi)
        if ovl_start > ovl_end:
            continue

        # 取外显子序列（已按外显子自身链方向给出）
        seq_json = requests.get(f"{server}/sequence/id/{exon['id']}", headers=headers).json()
        exon_seq = seq_json["seq"]

        # 索引换算：正链用 exon.start；负链用 exon.end
        if e_strand == 1:
            i0 = ovl_start - e_start
            i1 = ovl_end   - e_start + 1
        else:
            i0 = e_end - ovl_end
            i1 = e_end - ovl_start + 1

        part = exon_seq[i0:i1]
        cds_parts.append(part)
        total_bp += len(part)

    concatenated_cds = "".join(cds_parts)

    # 重点：此处不要再对负链做整体反向互补！
    final_manual_cds_dna = concatenated_cds

    if len(final_manual_cds_dna) % 3 != 0:
        print(f"  ! 警告：拼出的 CDS 长度不是 3 的整数倍：{len(final_manual_cds_dna)} bp")

    # 5) 翻译
    manual_protein_seq = translate_dna_std(final_manual_cds_dna)
    print(f"  > 自翻蛋白长度: {len(manual_protein_seq)} aa")

    # 6) 对比
    print("\n--- 序列片段预览 ---")
    print(f"  - 官方序列 (前120): {official_protein_seq[:]}")
    print(f"  - 自翻序列 (前120): {manual_protein_seq[:]}")
    print("\n匹配结果:", "一致 ✅" if manual_protein_seq == official_protein_seq else "不一致 ❌")

def verify_translation_logic(transcript_id: str):
    """
    接收一个转录本ID，通过两种方式获取其蛋白质序列并进行比对。
    其中“手动方式”将模拟真实的mRNA剪接和拼接过程来构建CDS。
    """
    print(f"\n===== 开始对转录本 {transcript_id} 进行双向翻译验证 =====")
    ensembl_server = "https://rest.ensembl.org"
    official_protein_seq = None
    manual_protein_seq = None

    # 1. “官方答案”：直接获取API翻译好的蛋白质序列
    print("\n[第1步] 正在获取 Ensembl 官方翻译的蛋白质序列...")
    ext_protein = f"/sequence/id/{transcript_id}?type=protein;content-type=application/json"
    try:
        r_prot = requests.get(ensembl_server + ext_protein)
        r_prot.raise_for_status()
        official_protein_seq = r_prot.json().get('seq')
        if official_protein_seq:
            print(f"  > 成功获取官方蛋白质序列，长度为 {len(official_protein_seq)} aa。")
        else:
            print("  > 此转录本可能不编码蛋白质。")
            return
    except requests.exceptions.RequestException as err:
        print(f"  > 获取官方蛋白质序列失败: {err}")
        return

    # --- “我们自己动手”的详细流程 ---
    print("\n[第2步] 开始手动模拟CDS序列的构建过程...")
    
    # 2a. 获取转录本的“说明书” (包含链向, CDS范围, 外显子列表)
    print("  - (2a) 获取转录本的完整信息...")
    ext_lookup = f"/lookup/id/{transcript_id}?expand=1;content-type=application/json"
    try:
        r_lookup = requests.get(ensembl_server + ext_lookup)
        r_lookup.raise_for_status()
        transcript_info = r_lookup.json()
        
        strand = transcript_info.get('strand')
        all_exons = transcript_info.get('Exon', [])
        cds_start = transcript_info.get('Translation', {}).get('start')
        cds_end = transcript_info.get('Translation', {}).get('end')

        if not all([strand, all_exons, cds_start, cds_end]):
            print("  > 错误：转录本信息不完整，无法继续手动构建。")
            return
        print("    > 信息获取成功。")

    except requests.exceptions.RequestException as err:
        print(f"  > 获取转录本信息失败: {err}")
        return

    # 2b. 筛选出所有与CDS区域有重叠的“编码外显子”
    print("  - (2b) 筛选编码外显子...")
    coding_exons = []
    for exon in all_exons:
        if max(exon['start'], cds_start) <= min(exon['end'], cds_end):
            coding_exons.append(exon)
    print(f"    > 找到 {len(coding_exons)} 个编码外显子。")

    # 2c. 按照生物学阅读顺序对外显子进行排序
    print("  - (2c) 按生物学顺序排序...")
    if strand == 1:
        coding_exons.sort(key=lambda x: x['start']) # 正向链，从小到大
    else:
        coding_exons.sort(key=lambda x: x['start'], reverse=True) # 反向链，从大到小
    print("    > 排序完成。")

    # 2d. 遍历排序后的外显子，获取序列、裁剪并拼接
    print("  - (2d) 获取、裁剪并拼接外显子序列...")
    cds_dna_fragments = []
    for exon in coding_exons:
        # 获取当前外显子的完整序列
        if strand == 1:
            ext_exon_seq = f"/sequence/id/{exon['id']}?content-type=application/json"
        else:
            # 反向链需要指定 strand=-1
            ext_exon_seq = f"/sequence/id/{exon['id']}?strand=-1;content-type=application/json"
        r_exon_seq = requests.get(ensembl_server + ext_exon_seq)
        full_exon_seq = r_exon_seq.json().get('seq')
        
        # 计算该外显子需要被裁剪的边界
        overlap_start = max(exon['start'], cds_start)
        overlap_end = min(exon['end'], cds_end)
        
        # 将基因组坐标转换为序列中的索引
        slice_start = overlap_start - exon['start']
        slice_end = overlap_end - exon['start'] + 1
        
        # 裁剪出属于CDS的部分
        cds_part_of_exon = full_exon_seq[slice_start:slice_end]
        cds_dna_fragments.append(cds_part_of_exon)
    
    # 拼接成一个“半成品”CDS
    concatenated_cds = "".join(cds_dna_fragments)
    print(f"    > 拼接完成，得到CDS半成品，长度 {len(concatenated_cds)} bp。")

    # 2e. 根据链向，对“半成品”进行最终的方向调整
    print("  - (2e) 根据链向调整序列方向...")
    final_manual_cds_dna = ""
    if strand == 1:
        final_manual_cds_dna = concatenated_cds
        print("    > 正向链，无需调整。")
    else:
        final_manual_cds_dna = reverse_complement(concatenated_cds)
        print("    > 反向链，已进行反向互补。")

    # 3. “我们自己翻译”：使用代码翻译手动构建的CDS
    print("\n[第3步] 正在使用本地代码翻译手动构建的CDS序列...")
    manual_protein_seq = translate_dna(final_manual_cds_dna)
    print(f"  > 本地代码翻译完成，得到蛋白质序列，长度为 {len(manual_protein_seq)} aa。")

    # 4. “最终比对”
    print("\n[第4步] 正在比对两种方式得到的蛋白质序列...")
    
    print("\n--- 序列片段预览 ---")
    print(f"  - 官方序列 (前120): {official_protein_seq[:120]}...")
    print(f"  - 自翻序列 (前120): {manual_protein_seq[:120]}...")

    if official_protein_seq == manual_protein_seq:
        print("\n✅ [验证成功] 您手动拼接并翻译的序列与官方蛋白质序列完全一致！")
        print("这证明了我们对整个流程的理解是完全正确的。")
    else:
        print("\n❌ [验证失败] 序列不一致！")
        print(f"  - 长度是否相同: 官方({len(official_protein_seq)}) vs 自翻({len(manual_protein_seq)})")

    print("=" * 60)

def analyze_exon_location(target_exon: dict, transcript_info: dict):
    """
    分析目标外显子在转录本中的位置 (UTR或CDS)。
    (修正版：以处理嵌套的Translation键名和特殊生物类型)
    """
    print("\n--- 外显子在转录本中的位置分析 ---")
    
    # 首先检查嵌套的 "Translation" 对象
    if 'Translation' in transcript_info and transcript_info['Translation']:
        translation_info = transcript_info['Translation']
        cds_start = translation_info.get('start')
        cds_end = translation_info.get('end')
        biotype = transcript_info.get('biotype', 'N/A')
    else:
        # 为其他可能使用顶级键的转录本保留备用逻辑
        cds_start = transcript_info.get('CDS_start')
        cds_end = transcript_info.get('CDS_end')
        biotype = transcript_info.get('biotype', 'N/A')

    if not cds_start or not cds_end:
        print(f"  - 生物类型: {biotype}")
        print("  - 结论: 该转录本为非编码转录本，其所有外显子均不翻译。")
        return

    exon_start = target_exon['start']
    exon_end = target_exon['end']
    
    print(f"  - 转录本生物类型: {biotype}")
    print(f"  - 转录本翻译范围 (CDS): {cds_start} - {cds_end}")
    print(f"  - 目标外显子范围: {exon_start} - {exon_end}")

    # 判断位置
    location = "位置关系未知，请检查坐标"
    if exon_end < cds_start:
        location = "完全位于 5' UTR (非翻译区)"
    elif exon_start > cds_end:
        location = "完全位于 3' UTR (非翻译区)"
    elif exon_start < cds_start and exon_end >= cds_start:
        location = "部分位于 5' UTR, 部分位于 CDS (包含起始密码子)"
    elif exon_start <= cds_end and exon_end > cds_end:
        location = "部分位于 CDS, 部分位于 3' UTR (包含终止密码子)"
    elif exon_start >= cds_start and exon_end <= cds_end:
        location = "完全位于 CDS (编码区)"
    
    print(f"  - 结论: 该外显子【{location}】。")

    # 为 "non_stop_decay" 添加特别注释
    if biotype == "non_stop_decay":
        print("  - 特别提示: 由于这是一个 'non_stop_decay' 类型的转录本，它缺少一个常规的终止密码子。从起始密码子之后的所有序列都会被翻译，直到转录本物理上结束。")

def get_uniprot_id_for_gene(gene_id: str) -> str or None:
    """第一步: 根据Ensembl基因ID获取其对应的UniProt/Swiss-Prot ID"""
    print(f"\n[步骤1] 正在为基因 {gene_id} 查询UniProt ID...")
    ensembl_server = "https://rest.ensembl.org"
    
    # 修正：URL参数应该用查询字符串格式，而不是分号
    ext_xrefs = f"/xrefs/id/{gene_id}?external_db=UniProt/SWISSPROT&content-type=application/json"
    
    try:
        r = requests.get(ensembl_server + ext_xrefs)
        r.raise_for_status()
        data = r.json()
        
        if data:
            # 可能返回多个结果，优先选择 Swiss-Prot 条目
            for entry in data:
                if entry.get('dbname') == 'UniProt/SWISSPROT':
                    uniprot_id = entry['primary_id']
                    print(f"  > 成功找到UniProt ID: {uniprot_id}")
                    return uniprot_id
            
            # 如果没有 Swiss-Prot，尝试 TrEMBL
            for entry in data:
                if entry.get('dbname') == 'UniProt/TrEMBL':
                    uniprot_id = entry['primary_id']
                    print(f"  > 找到UniProt/TrEMBL ID: {uniprot_id}")
                    return uniprot_id
                    
        print("  > 未找到对应的UniProt ID。")
        return None
        
    except requests.exceptions.RequestException as err:
        print(f"  > API请求失败: {err}")
        return None

def get_uniprot_id_for_gene_debug(gene_id: str) -> str or None:
    """调试版本：打印更多信息"""
    print(f"\n[步骤1] 正在为基因 {gene_id} 查询UniProt ID...")
    
    # 尝试不同的查询方式
    endpoints = [
        f"/xrefs/id/{gene_id}",  # 获取所有外部引用
        f"/xrefs/id/{gene_id}?external_db=UniProt/SWISSPROT",
        f"/xrefs/id/{gene_id}?external_db=UniProt"
    ]
    
    for ext in endpoints:
        try:
            print(f"  尝试端点: {ext}")
            r = requests.get(ENSEMBL_SERVER + ext, headers=JSON_HEADERS)
            r.raise_for_status()
            data = r.json()
            
            print(f"  返回数据条数: {len(data)}")
            if data:
                # 打印前几条数据以供调试
                for i, entry in enumerate(data[:3]):
                    print(f"  条目{i+1}: dbname={entry.get('dbname')}, primary_id={entry.get('primary_id')}")
                
                # 查找 UniProt 条目
                for entry in data:
                    if 'uniprot' in entry.get('dbname', '').lower():
                        return entry['primary_id']
                        
        except Exception as e:
            print(f"  该端点失败: {e}")
            
    return None

def get_protein_sequence_from_uniprot(uniprot_id: str) -> str or None:
    """第二步: 根据UniProt ID获取FASTA格式的蛋白质序列"""
    print(f"[步骤2] 正在从UniProt获取 {uniprot_id} 的蛋白质序列...")
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    try:
        r = requests.get(url)
        r.raise_for_status()
        # 解析FASTA格式，去掉第一行的头部信息
        fasta_lines = r.text.strip().split('\n')
        protein_seq = "".join(fasta_lines[1:])
        print(f"  > 成功获取序列，长度为 {len(protein_seq)} aa。")
        return protein_seq
    except requests.exceptions.RequestException as err:
        print(f"  > 请求UniProt失败: {err}")
        return None

def find_transcript_by_protein_sequence(gene_id: str, target_protein_seq: str) -> dict or None:
    """第三步: 在基因的所有转录本中，找到能翻译出目标蛋白序列的那个"""
    print(f"[步骤3] 正在基因 {gene_id} 中反向匹配转录本...")
    
    # 获取基因下所有转录本
    ext_lookup = f"/lookup/id/{gene_id}?expand=1"
    try:
        r_gene = requests.get(ENSEMBL_SERVER + ext_lookup, headers=JSON_HEADERS)
        r_gene.raise_for_status()
        transcripts = r_gene.json().get('Transcript', [])
    except requests.exceptions.RequestException as err:
        print(f"  > 获取基因信息失败: {err}")
        return None

    for transcript in transcripts:
        t_id = transcript['id']
        # 获取该转录本翻译的蛋白质序列
        ext_prot = f"/sequence/id/{t_id}?type=protein"
        try:
            r_prot = requests.get(ENSEMBL_SERVER + ext_prot, headers=JSON_HEADERS)
            r_prot.raise_for_status()
            current_prot_seq = r_prot.json().get('seq')
            
            # 序列比对
            if current_prot_seq == target_protein_seq:
                print(f"  > 匹配成功！转录本 {t_id} 翻译的蛋白质与目标序列一致。")
                return transcript
        except requests.exceptions.RequestException:
            # 某些转录本可能不编码蛋白质，请求会失败，直接跳过
            continue
            
    print("  > 未能找到与UniProt序列完全匹配的转录本。")
    return None

def reverse_complement(dna_seq: str) -> str:
    """计算DNA序列的反向互补序列"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(dna_seq))

def translate_cds(cds_seq: str) -> str:
    """将CDS序列翻译成蛋白质序列"""
    protein_seq = []
    for i in range(0, len(cds_seq) - (len(cds_seq) % 3), 3):
        codon = cds_seq[i:i+3]
        amino_acid = CODON_TABLE.get(codon, 'X') # 'X' for unknown codon
        protein_seq.append(amino_acid)
    return "".join(protein_seq)

def translate_dna(dna_seq: str) -> str:
    """将DNA序列翻译成蛋白质序列。"""
    protein = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3].upper()
        amino_acid = CODON_TABLE.get(codon, 'X') # 'X' for unknown codons
        protein.append(amino_acid)
    return "".join(protein)

def translate_single_exon(target_exon: dict, transcript_info: dict, ref_allele: str, alt_allele: str,variant_position, relative_pos,desc_id, gene_name, gene_id) -> None:
    """第四步: 使用匹配到的转录本信息，翻译单个外显子"""
    print(f"[步骤4] 正在翻译目标外显子 (ID: {target_exon.get('id')})...")
    
    strand = transcript_info.get('strand')
    
    # 1. 从转录本信息中获取准确的CDS起始和终止坐标
    cds_start, cds_end = 0, 0
    if 'Translation' in transcript_info and transcript_info['Translation']:
        translation_info = transcript_info['Translation']
        cds_start = translation_info.get('start')
        cds_end = translation_info.get('end')
    else:
        print("  > 错误：转录本信息中缺少'Translation'对象，无法计算准确相位。")
        return

    # 2. 将所有外显子按照生物学阅读顺序排序
    all_exons = sorted(transcript_info.get('Exon', []), key=lambda x: x['start'])
    if strand == -1:
        all_exons.reverse() # 反向链从高坐标到低坐标，所以反转列表

    # 3. 遍历排序后的外显子，只累加目标外显子之前的、位于CDS区域内的碱基长度
    cds_len_before_exon = 0
    for current_exon in all_exons:
        if current_exon['id'] == target_exon['id']:
            break

        overlap_start = max(current_exon['start'], cds_start)
        overlap_end = min(current_exon['end'], cds_end)

        if overlap_start <= overlap_end:
            coding_length_of_this_exon = overlap_end - overlap_start + 1
            cds_len_before_exon += coding_length_of_this_exon
            
    phase = cds_len_before_exon % 3
    print(f"  > [修正后] 根据CDS区域计算得到外显子起始相位(Phase): {phase}")

    # 获取外显子序列
    ensembl_server = "https://rest.ensembl.org"
    ext_seq = f"/sequence/id/{target_exon['id']}?content-type=application/json"
    try:
        r_seq = requests.get(ensembl_server + ext_seq)
        r_seq.raise_for_status()
        exon_seq = r_seq.json().get('seq')
    except requests.exceptions.RequestException as err:
        print(f"  > 获取外显子序列失败: {err}")
        return

    # 根据链向处理序列和变异
    local_exon_seq = exon_seq
    local_ref = ref_allele
    local_alt = alt_allele
    if strand == -1:
        local_exon_seq = reverse_complement(exon_seq)
        local_ref = reverse_complement(ref_allele)
        local_alt = reverse_complement(alt_allele)

    # 引入变异
    variant_pos_in_genome = int(target_exon['variant_pos'])
    
    local_relative_pos = 0
    if strand == 1:
        local_relative_pos = variant_pos_in_genome - target_exon['start']
    else: # 反向链
        local_relative_pos = target_exon['end'] - variant_pos_in_genome
        
    mutated_exon_seq_for_translation = local_exon_seq[:local_relative_pos] + local_alt + local_exon_seq[local_relative_pos + len(local_ref):]

    # 根据相位进行翻译
    ref_to_translate = local_exon_seq[phase:]
    mut_to_translate = mutated_exon_seq_for_translation[phase:]
    ref_protein_fragment = translate_dna(ref_to_translate)
    mut_protein_fragment = translate_dna(mut_to_translate)
    
    # --- 新增：完整打印外显子翻译的蛋白质片段 ---
    print("\n  --- 单外显子翻译完整蛋白序列 ---")
    print(f"  > 原始外显子蛋白片段 ({len(ref_protein_fragment)} aa):")
    for i in range(0, len(ref_protein_fragment), 60):
        print(f"    {ref_protein_fragment[i:i+60]}")
    
    print(f"  > 突变后外显子蛋白片段 ({len(mut_protein_fragment)} aa):")
    for i in range(0, len(mut_protein_fragment), 60):
        print(f"    {mut_protein_fragment[i:i+60]}")
    # ------------------------------------

    print("\n  --- 蛋白质变异分析 ---")
    if ref_protein_fragment == mut_protein_fragment:
        print("  - 变异类型: 同义突变 (Synonymous)")
        print("  - 变异结论: 蛋白质序列未发生改变。")
    else:
        # 计算氨基酸位置
        aa_pos_in_fragment = (local_relative_pos - phase) // 3

        # 安全检查
        if aa_pos_in_fragment < len(ref_protein_fragment) and aa_pos_in_fragment < len(mut_protein_fragment):
            ref_aa = ref_protein_fragment[aa_pos_in_fragment]
            mut_aa = mut_protein_fragment[aa_pos_in_fragment]

            mutation_type = "未知改变"
            if len(ref_protein_fragment) != len(mut_protein_fragment):
                mutation_type = "移码突变 (Frameshift)"
            elif ref_aa == '_' and mut_aa != '_':
                mutation_type = "终止密码子丢失 (Stop-loss)"
            elif ref_aa != '_' and mut_aa == '_':
                mutation_type = "无义突变 (Nonsense)"
            elif ref_aa != mut_aa:
                mutation_type = "错义突变 (Missense)"

            print(f"  - 变异类型: {mutation_type}")
            print(f"  - 影响位置: 该外显子翻译片段的第 {aa_pos_in_fragment + 1} 个氨基酸")
            print(f"  - 具体改变: {ref_aa} -> {mut_aa}")
        
            print(f"  - 原始片段: ...{ref_protein_fragment[max(0, aa_pos_in_fragment-10) : aa_pos_in_fragment]}[{ref_aa}]{ref_protein_fragment[aa_pos_in_fragment+1: aa_pos_in_fragment+11]}...")
            print(f"  - 突变后片段: ...{mut_protein_fragment[max(0, aa_pos_in_fragment-10) : aa_pos_in_fragment]}[{mut_aa}]{mut_protein_fragment[aa_pos_in_fragment+1: aa_pos_in_fragment+11]}...")
        else:
            print("  - 变异分析错误: 计算出的氨基酸位置超出了翻译片段的范围。")

    # 将原始（正向链）序列传递给报告函数
    mutated_exon_seq_for_report = exon_seq[:relative_pos] + alt_allele + exon_seq[relative_pos + len(ref_allele):]
    print_exon_report(desc_id, gene_name, gene_id, target_exon, exon_seq, mutated_exon_seq_for_report, variant_position, relative_pos, ref_allele, alt_allele)

def print_exon_derived_protein_segment(transcript_info: dict, target_exon: dict):
    """
    打印该转录本最终蛋白序列中，由目标外显子“编码部分”对应的蛋白片段。
    同时提供两种视图：
      - inclusive: 包含跨界氨基酸（该外显子至少贡献1或2个核苷酸）
      - strict: 仅包含完全由该外显子内部3个一组独立编码的氨基酸
    """
    cds_seq, protein_seq, segments = build_cds_and_segment_map(transcript_info)
    t_id = transcript_info["id"]
    prot_len = len(protein_seq)

    # 定位目标外显子在 segments 中（segments 只包含编码外显子）
    seg = next((s for s in segments if s["exon_id"] == target_exon["id"]), None)
    print("\n--- 外显子→蛋白片段映射 ---")
    print(f"转录本: {t_id}  (蛋白长度: {prot_len} aa)")
    if not seg:
        print("  - 该外显子不在CDS中（完全UTR），没有对应的蛋白片段。")
        return

    offset = seg["cds_offset"]      # 该外显子编码片段在CDS中的0-based nt偏移
    L = seg["length_nt"]            # 该外显子在CDS中的nt长度
    phase = offset % 3              # 外显子编码片段起始处的相位（0,1,2）

    # 视图1：inclusive（包含跨界氨基酸）
    aa_start_incl = (offset // 3) + 1
    aa_end_incl = ((offset + L - 1) // 3) + 1
    aa_start_incl = max(1, aa_start_incl)
    aa_end_incl = min(prot_len, aa_end_incl)

    incl_seg = protein_seq[aa_start_incl - 1: aa_end_incl] if aa_start_incl <= aa_end_incl else ""

    # 视图2：strict（只保留完全由该外显子内部3nt独立编码出的aa）
    n_skip = (3 - phase) % 3
    if L <= n_skip:
        strict_seg = ""
        aa_start_strict = None
        aa_end_strict = None
    else:
        start_nt = offset + n_skip
        n_full = (L - n_skip) // 3
        aa_start_strict = (start_nt // 3) + 1
        aa_end_strict = aa_start_strict + n_full - 1
        strict_seg = protein_seq[aa_start_strict - 1: aa_end_strict] if n_full > 0 else ""

    print(f"  - 外显子ID: {target_exon['id']}")
    print(f"  - 外显子CDS相对相位 phase = {phase}")
    print(f"  - inclusive 蛋白区间: aa {aa_start_incl}..{aa_end_incl}  (含跨界aa)")
    print(f"    序列: {incl_seg if incl_seg else '(空)'}")
    if strict_seg:
        print(f"  - strict    蛋白区间: aa {aa_start_strict}..{aa_end_strict}  (完全由此外显子内部3nt编码)")
        print(f"    序列: {strict_seg}")
    else:
        print("  - strict    蛋白区间: (该外显子的编码长度不足以形成完整密码子，strict 片段为空)")

def wrap_seq(seq: str, width: int = 60):
    for i in range(0, len(seq), width):
        yield seq[i:i+width]

def fetch_transcript_full(transcript_id: str) -> dict:
    """获取转录本的完整信息（包含 Exon 列表与 Translation 基因组坐标）"""
    r = requests.get(f"{ENSEMBL_SERVER}/lookup/id/{transcript_id}?expand=1", headers=JSON_HEADERS)
    r.raise_for_status()
    return r.json()

def build_cds_and_segment_map(transcript_info: dict):
    """
    拼出CDS并建立“编码外显子→CDS偏移”的映射。
    返回:
      - cds_seq: 拼接后的CDS序列（转录本方向）
      - protein_seq_official: 官方蛋白序列（Ensembl翻译）
      - segments: 每个编码外显子的映射列表（按阅读顺序）
    """
    strand = transcript_info["strand"]
    exons = transcript_info["Exon"]
    t = transcript_info.get("Translation")
    if not t:
        raise ValueError("该转录本没有 Translation（非编码）")
    cds_g_start, cds_g_end = t["start"], t["end"]
    cds_lo, cds_hi = min(cds_g_start, cds_g_end), max(cds_g_start, cds_g_end)

    # 官方蛋白
    r_prot = requests.get(f"{ENSEMBL_SERVER}/sequence/id/{transcript_info['id']}?type=protein", headers=JSON_HEADERS)
    r_prot.raise_for_status()
    protein_seq_official = r_prot.json()["seq"]

    # 编码外显子，按阅读方向排序
    coding_exons = [e for e in exons if max(e['start'], cds_lo) <= min(e['end'], cds_hi)]
    if strand == 1:
        coding_exons.sort(key=lambda x: x['start'])
    else:
        coding_exons.sort(key=lambda x: x['start'], reverse=True)

    cds_parts = []
    segments = []
    cum = 0  # CDS累计nt长度

    for e in coding_exons:
        e_start, e_end, e_strand = e['start'], e['end'], e['strand']
        ovl_start = max(e_start, cds_lo)
        ovl_end   = min(e_end, cds_hi)
        if ovl_start > ovl_end:
            continue

        # 外显子序列（已按外显子自身链方向）
        seq_json = requests.get(f"{ENSEMBL_SERVER}/sequence/id/{e['id']}", headers=JSON_HEADERS).json()
        exon_seq = seq_json["seq"]

        # 切片索引：正链相对 exon.start；负链相对 exon.end
        if e_strand == 1:
            i0 = ovl_start - e_start
            i1 = ovl_end   - e_start + 1
        else:
            i0 = e_end - ovl_end
            i1 = e_end - ovl_start + 1

        part = exon_seq[i0:i1]
        cds_parts.append(part)

        segments.append({
            "exon_id": e["id"],
            "exon_start": e_start,
            "exon_end": e_end,
            "exon_strand": e_strand,
            "ovl_start": ovl_start,
            "ovl_end": ovl_end,
            "i0": i0,  # 该编码片段在外显子序列中的起始索引
            "length_nt": len(part),
            "cds_offset": cum  # 在CDS中的0-based nt偏移
        })
        cum += len(part)

    cds_seq = "".join(cds_parts)
    return cds_seq, protein_seq_official, segments

def find_exon_on_transcript(transcript_info: dict, variant_position: int) -> dict or None:
    """在给定转录本的外显子中找到包含 variant_position 的那个外显子。"""
    exons = transcript_info.get("Exon", [])
    for i, e in enumerate(exons, start=1):
        if e["start"] <= variant_position <= e["end"]:
            ee = dict(e)
            ee["number"] = i
            ee["total_count"] = len(exons)
            ee["variant_pos"] = variant_position
            return ee
    return None

def report_full_protein_and_exon_effect(transcript_full: dict, target_exon: dict,
                                        ref: str, alt: str, variant_position: int):
    """
    打印：
      1) 整个基因（该转录本）所翻译的蛋白序列
      2) 当前外显子（严格strict）突变前的蛋白片段
      3) 当前外显子（严格strict）突变后的蛋白片段
      4) 只针对该外显子片段的影响分析
    """
    print("\n================= 蛋白与外显子片段报告 =================")

    # 1) 拼CDS并拿到官方全蛋白、segments
    cds_seq, protein_seq, segments = build_cds_and_segment_map(transcript_full)
    prot_len = len(protein_seq)
    print(f"[全蛋白] 长度: {prot_len} aa")
    for line in wrap_seq(protein_seq, 60):
        print(f"  {line}")

    print(f"Debug_Segments_Info:{segments}")
    # 2) 定位目标外显子的编码片段
    seg = next((s for s in segments if s["exon_id"] == target_exon["id"]), None)
    if not seg:
        print("\n[提示] 该外显子不在CDS中（UTR），因此没有对应的蛋白片段。")
        return

    e_start = seg["exon_start"]
    e_end   = seg["exon_end"]
    e_strand = seg["exon_strand"]
    i0 = seg["i0"]               # 编码片段在外显子序列中的起始索引
    L  = seg["length_nt"]        # 编码长度（nt）
    offset = seg["cds_offset"]   # 该片段在CDS中的0-based nt偏移
    phase = offset % 3
    n_skip = (3 - phase) % 3     # 丢弃到最近密码子边界的碱基数
    n_full = (L - n_skip) // 3

    # 外显子“严格”片段在全蛋白中的aa坐标
    if n_full <= 0:
        print("\n[提示] 该外显子在CDS中的编码长度不足以形成完整的密码子（strict 片段为空）。")
        return

    aa_start_strict = ((offset + n_skip) // 3) + 1  # 1-based
    aa_end_strict   = aa_start_strict + n_full - 1

    ref_exon_aa_strict = protein_seq[aa_start_strict - 1: aa_end_strict]

    # 3) 构建“突变后”的外显子编码片段核苷酸序列并翻译
    # 取外显子全序列（外显子自身方向）
    seq_json = requests.get(f"{ENSEMBL_SERVER}/sequence/id/{target_exon['id']}", headers=JSON_HEADERS).json()
    exon_seq = seq_json["seq"]
    # 编码部分
    part = exon_seq[i0:i0 + L]

    # 变异相对外显子序列的索引
    if not (e_start <= variant_position <= e_end):
        # 变异不在该外显子上（极少见于多转录本不一致时）
        idx_in_exon = None
    else:
        idx_in_exon = (variant_position - e_start) if e_strand == 1 else (e_end - variant_position)

    # 仅当变异落在该外显子的编码重叠区间内，才会影响 part
    mutated_part = part
    if idx_in_exon is not None:
        # 变异落在编码区间时的 part 内相对索引
        idx_in_part = idx_in_exon - i0
        if 0 <= idx_in_part <= len(part):
            # 将 REF/ALT 转成与外显子序列方向一致的等效碱基
            local_ref = ref
            local_alt = alt
            if e_strand == -1:
                local_ref = reverse_complement(ref)
                local_alt = reverse_complement(alt)

            # 检查是否真的在编码区间内
            if seg["ovl_start"] <= variant_position <= seg["ovl_end"]:
                mutated_part = part[:idx_in_part] + local_alt + part[idx_in_part + len(local_ref):]

    # 严格视图：从第一个完整密码子开始翻译
    ref_part_strict = part[n_skip:]
    mut_part_strict = mutated_part[n_skip:]

    ref_exon_aa_strict_manual = translate_dna(ref_part_strict)
    mut_exon_aa_strict = translate_dna(mut_part_strict)

    # 4) 打印“外显子片段（突变前/后）”并做影响分析（仅限该外显子片段）
    print("\n[外显子严格片段 | 突变前] (全蛋白位置: aa {}..{})".format(aa_start_strict, aa_end_strict))
    # 用官方全蛋白切出来的片段作为“参考真相”
    for line in wrap_seq(ref_exon_aa_strict, 60):
        print(f"  {line}")

    print("\n[外显子严格片段 | 突变后] (起点对齐到同一全蛋白位置 aa {})".format(aa_start_strict))
    for line in wrap_seq(mut_exon_aa_strict, 60):
        print(f"  {line}")

    # 影响类型判定（仅限该严格片段）
    print("\n[只看该外显子片段的影响]")
    if mut_exon_aa_strict == ref_exon_aa_strict:
        print("  - 变异类型: 同义 (Synonymous) / 片段未改变")
        return

    # 判断是否移码（看严格片段翻译所用的nt长度变化是否为3的倍数）
    # 注意：这里只比较该片段内部的变化（不跨外显子）
    delta_nt = len(mut_part_strict) - len(ref_part_strict)
    if delta_nt % 3 != 0:
        effect_type = "移码突变 (Frameshift)"
    else:
        # 看是否出现/消失终止（你的translate_dna以'_'为终止）
        if "_" in mut_exon_aa_strict and "_" not in ref_exon_aa_strict:
            effect_type = "无义突变 (Nonsense)"
        elif "_" not in mut_exon_aa_strict and "_" in ref_exon_aa_strict:
            effect_type = "终止丢失 (Stop-loss)"
        elif len(mut_exon_aa_strict) != len(ref_exon_aa_strict):
            effect_type = "插入/缺失（不移码）(In-frame indel)"
        else:
            effect_type = "错义突变 (Missense)"
    print(f"  - 变异类型: {effect_type}")

    # 找到片段内首个差异AA并报告对应的全蛋白位置
    i = 0
    min_len = min(len(mut_exon_aa_strict), len(ref_exon_aa_strict))
    while i < min_len and mut_exon_aa_strict[i] == ref_exon_aa_strict[i]:
        i += 1
    if i < min_len:
        aa_global = aa_start_strict + i
        print(f"  - 首个差异位点（片段内）: 第 {i+1} 个氨基酸（全蛋白位置 aa {aa_global}）")
        # 打印一个小窗口
        left = max(0, i-10)
        right = i+11
        ref_win = ref_exon_aa_strict[left:right]
        mut_win = mut_exon_aa_strict[left:right]
        # 用中括号标记差异位点
        if i-left < len(ref_win):
            ref_win = ref_win[:i-left] + "[" + ref_win[i-left:i-left+1] + "]" + ref_win[i-left+1:]
        if i-left < len(mut_win):
            mut_win = mut_win[:i-left] + "[" + mut_win[i-left:i-left+1] + "]" + mut_win[i-left+1:]
        print(f"  - 参考窗口: ...{ref_win}...")
        print(f"  - 突变窗口: ...{mut_win}...")
    else:
        # 长度不同但前缀一致的情况
        aa_global = aa_start_strict + min_len
        print(f"  - 片段长度改变，分歧点在末端（全蛋白位置起于 aa {aa_global} 之后）")

    print("=========================================================\n")

def process_variant_with_uniprot_workflow(gene_id: str, target_exon_info: dict, ref: str, alt: str,
                                          desc_id, gene_name, variant_position, relative_pos):
    # 步骤1~3保持不变（UniProt→全蛋白→匹配转录本）
    uniprot_id = get_uniprot_id_for_gene_debug(gene_id)
    if not uniprot_id:
        return

    protein_seq = get_protein_sequence_from_uniprot(uniprot_id)
    if not protein_seq:
        return

    print(f"\n--- UniProt ID ({uniprot_id}) 对应的蛋白质序列 ---")
    for i in range(0, len(protein_seq), 60):
        print(f"  {protein_seq[i:i+60]}")
    print("-" * 40)

    matched_transcript = find_transcript_by_protein_sequence(gene_id, protein_seq)
    if not matched_transcript:
        return

    t_id = matched_transcript["id"]
    print(f"\n--- 反向匹配到的转录本信息 ---")
    print(f"  - 转录本ID: {t_id}")
    print("-" * 40)

    # 用匹配的转录本做后续映射
    transcript_full = fetch_transcript_full(t_id)

    # 在该转录本内重新定位目标外显子（避免与canonical不一致）
    target_exon_on_matched = find_exon_on_transcript(transcript_full, variant_position)
    if not target_exon_on_matched:
        print("  > 提示：该变异不位于匹配转录本的任何外显子上（可能处于内含子/UTR）。")
        return

    # UTR/CDS 分析（保留）
    analyze_exon_location(target_exon_on_matched, transcript_full)

    # 新增：按你的最终需求打印“全蛋白 + 外显子片段（前/后）+ 仅该片段的影响”
    report_full_protein_and_exon_effect(transcript_full, target_exon_on_matched, ref, alt, variant_position)
    
def classify_variant(ref: str, alt: str) -> str:
    """根据REF和ALT碱基的长度判断基础的变异类型。"""
    if len(ref) == 1 and len(alt) == 1:
        return "SNP (单碱基多态性)"
    elif len(ref) > len(alt):
        return "Deletion (缺失)"
    elif len(ref) < len(alt):
        return "Insertion (插入)"
    else:
        return "Complex (复杂变异)"


def generate_variant_context_report(desc_id: str):
    """
    解析 "CHROM_POS_REF_ALT" 格式的ID，并根据变异位置生成上下文序列报告。
    """
    print("\n--- 开始处理该位点的ID上下文序列 ---")
    print(f"===== 开始处理ID: {desc_id} =====")
    
    try:
        parts = desc_id.split('_')
        chrom_code, pos_str, ref_allele, alt_allele = parts
        variant_position = int(pos_str)
        chromosome_map = {'24': 'Y', '23': 'X'}
        chromosome = chromosome_map.get(chrom_code, chrom_code)
    except (ValueError, IndexError):
        print(f"错误：ID '{desc_id}' 格式不正确，已跳过序列分析。")
        return
    
    gene_name, gene_id = "N/A", "N/A"
    try:
        ext_overlap = f"/overlap/region/human/{chromosome}:{variant_position}-{variant_position}?feature=gene"
        r_overlap = requests.get(ENSEMBL_SERVER + ext_overlap, headers=JSON_HEADERS)
        r_overlap.raise_for_status()
        overlap_data = r_overlap.json()
        if overlap_data:
            gene_info = overlap_data[0]
            gene_id, gene_name = gene_info.get('id'), gene_info.get('external_name', 'N/A')
    except requests.exceptions.RequestException:
        pass # 静默处理错误，后续流程会按基因间区处理
    
    '''寻找对应的外显子'''
    target_exon = None
    if gene_id != "N/A":
        try:
            ext_lookup = f"/lookup/id/{gene_id}?expand=1"
            r_lookup = requests.get(ENSEMBL_SERVER + ext_lookup, headers=JSON_HEADERS)
            r_lookup.raise_for_status()
            gene_data = r_lookup.json()
            transcript = next((t for t in gene_data.get('Transcript', []) if t.get('is_canonical')), 
                              gene_data.get('Transcript', [{}])[0])
            all_exons = transcript.get('Exon', [])
            for i, exon in enumerate(all_exons):
                if exon['start'] <= variant_position <= exon['end']:
                    target_exon = exon
                    target_exon['number'] = i + 1
                    target_exon['total_count'] = len(all_exons)
                    break
        except (requests.exceptions.RequestException, IndexError, KeyError):
            pass # 静默处理错误，后续流程会按内含子处理

    if target_exon:
        target_exon['variant_pos'] = variant_position
        exon_chr, exon_start, exon_end = target_exon['seq_region_name'], target_exon['start'], target_exon['end']
        relative_pos = variant_position - exon_start
        '''
        执行从基因到单外显子翻译的完整UniProt锚定工作流
        '''
        process_variant_with_uniprot_workflow(gene_id, target_exon, ref_allele, alt_allele,desc_id,gene_name,variant_position, relative_pos)
    else:
        flank_size = 30
        ref_len = len(ref_allele)
        start_fetch = variant_position - flank_size
        end_fetch = (variant_position + ref_len - 1) + flank_size
        try:
            ext_sequence = f"/sequence/region/human/{chromosome}:{start_fetch}-{end_fetch}"
            r_seq = requests.get(ENSEMBL_SERVER + ext_sequence, headers=JSON_HEADERS)
            r_seq.raise_for_status()
            sequence = r_seq.json().get('seq')
            print_flanking_report(desc_id, gene_name, gene_id, chromosome, variant_position, sequence, flank_size, ref_allele, alt_allele)
        except requests.exceptions.RequestException as err:
            print(f"错误：获取序列片段失败。原因: {err}")

def print_exon_report(desc_id, gene_name, gene_id, exon, sequence, mutated_exon_seq,variant_pos, relative_pos, ref_allele, alt_allele):
    """格式化打印外显子报告，并新增序列比对功能"""
    print(f"\n======================================================================")
    print(f"      外显子 (Exon) 上下文序列报告")
    print(f"      查询ID: {desc_id}")
    print(f"      所属基因: {gene_name} ({gene_id})")
    print(f"======================================================================")
    print(f"\n【目标外显子信息】")
    print(f"  - 编号:       外显子 {exon['number']} (共 {exon['total_count']} 个)")
    print(f"  - 坐标:       chr{exon['seq_region_name']}:{exon['start']}-{exon['end']}")
    
    # mutated_sequence = sequence[:relative_pos] + alt_allele + sequence[relative_pos + len(ref_allele):]

    print("\n【序列比对】")
    print(f"  - 原始外显子序列 (Reference):")
    for i in range(0, len(sequence), 60):
         print(f"    {sequence[i:i+60]}")
    print(f"  - 突变后外显子序列 (Mutated):")
    for i in range(0, len(mutated_exon_seq), 60):
         print(f"    {mutated_exon_seq[i:i+60]}")
    print(f"\n======================================================================")
    print("报告结束。")
    
def print_flanking_report(desc_id, gene_name, gene_id, chrom, var_pos, sequence, flank, ref_allele, alt_allele):
    """格式化打印内含子/基因间区报告，并新增序列比对功能"""
    context_type = "内含子" if gene_name != "N/A" else "基因间区"
    print(f"\n======================================================================")
    print(f"      {context_type} 上下文序列报告")
    print(f"      查询ID: {desc_id}")
    print(f"      所属基因: {gene_name} ({gene_id})")
    print(f"======================================================================")
    print(f"\n【目标区域信息】")
    print(f"  - 查询坐标:   chr{chrom}:{var_pos}")
    print(f"  - 侧翼长度:   每侧 {flank} bp")

    mutated_sequence = sequence[:flank] + alt_allele + sequence[flank + len(ref_allele):]
    
    print("\n【序列比对】")
    print(f"  - 原始片段序列 (Reference):")
    print(f"    {sequence}")
    print(f"  - 突变后片段序列 (Mutated):")
    print(f"    {mutated_sequence}")
    print(f"\n======================================================================")
    print("报告结束。")

def read_vcf_and_process(file_path: str, record_limit: int = 5):
    """
    读取VCF文件，按照指定格式打印每条记录，并调用序列上下文分析。
    """
    try:
        reader = vcfpy.Reader.from_path(file_path)
    except FileNotFoundError:
        print(f"错误：VCF文件 '{file_path}' 未找到。")
        return

    # --- 严格按照您的格式打印元信息 ---
    print("---VCF元信息（INFO描述）---")
    for line in reader.header.lines:
        if line.key == "INFO":
            print(f"{line.id}: {line.description}")

    print(f"\n---VCF数据记录（前{record_limit}条）---")
    for i, record in enumerate(reader):
        if record_limit is not None and i >= record_limit:
            break
        
        # --- 严格按照您的格式打印VCF记录 ---
        alt_str = record.ALT[0].value if record.ALT else ""
        variant_type = classify_variant(record.REF, alt_str)

        print(f"变异位点: {record.CHROM}:{record.POS}")
        print(f"ID: {record.ID}")
        print(f"REF: {record.REF}, ALT: {record.ALT}")
        print(f"突变类型: {variant_type}")  # <-- 新增信息
        print(f"质量 (QUAL): {record.QUAL}")
        print(f"过滤状态 (FILTER): {record.FILTER}")

        print("INFO字段:")
        for key, value in record.INFO.items():
            print(f"  {key}: {value}")

        # --- 调用序列分析模块 ---
        if record.ID:
            for single_id in record.ID:
                if '_' in single_id and len(single_id.split('_')) == 4:
                    generate_variant_context_report(single_id)
                else:
                    print(f"\nID '{single_id}' 不符合 'CHROM_POS_REF_ALT' 格式，已跳过序列分析。")
        else:
            print("\n该记录ID字段为空，跳过序列上下文分析。")
        
        print("\n" + "-"*20)

# if __name__ == "__main__":
#     vcf_file_path = "rgc_me_variant_frequencies_chrY_20231004.vcf"
    
#     # 默认只处理文件的前3条记录，您可以修改此数字
#     read_vcf_and_process(file_path=vcf_file_path, record_limit=3)

# 调试脚本
if __name__ == "__main__":
    
     generate_variant_context_report("24_13366334_T_C")
    #generate_variant_context_report("24_5737391_C_A")
