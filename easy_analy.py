import vcfpy
import requests
from datetime import datetime
import json

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

# =================================================================================
#  模块 5: UniProt锚定工作流
# =================================================================================

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
    ensembl_server = "https://rest.ensembl.org"
    
    # 尝试不同的查询方式
    endpoints = [
        f"/xrefs/id/{gene_id}?content-type=application/json",  # 获取所有外部引用
        f"/xrefs/id/{gene_id}?external_db=UniProt/SWISSPROT&content-type=application/json",
        f"/xrefs/id/{gene_id}?external_db=UniProt&content-type=application/json"
    ]
    
    for ext in endpoints:
        try:
            print(f"  尝试端点: {ext}")
            r = requests.get(ensembl_server + ext)
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
    ensembl_server = "https://rest.ensembl.org"
    
    # 获取基因下所有转录本
    ext_lookup = f"/lookup/id/{gene_id}?expand=1;content-type=application/json"
    try:
        r_gene = requests.get(ensembl_server + ext_lookup)
        r_gene.raise_for_status()
        transcripts = r_gene.json().get('Transcript', [])
    except requests.exceptions.RequestException as err:
        print(f"  > 获取基因信息失败: {err}")
        return None

    for transcript in transcripts:
        t_id = transcript['id']
        # 获取该转录本翻译的蛋白质序列
        ext_prot = f"/sequence/id/{t_id}?type=protein;content-type=application/json"
        try:
            r_prot = requests.get(ensembl_server + ext_prot)
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
    
    # print(f"Debug:TRAN_INFO--{transcript_info}")
    # 1. 计算外显子在CDS中的起始相位(phase)
    cds_len_before_exon = 0
    cds_exons = sorted([exon for exon in transcript_info.get('Exon', []) if exon.get('object_type') == 'Exon'], 
                       key=lambda x: x['start'])
    # print(f"Debug_CDS_EXONS:{cds_exons}")
    if strand == -1:
        cds_exons = sorted(cds_exons, key=lambda x: x['start'], reverse=True)
        
    for exon in cds_exons:
        if exon['id'] == target_exon['id']:
            break
        cds_len_before_exon += (exon['end'] - exon['start'] + 1)
    
    phase = cds_len_before_exon % 3
    print(f"  > 计算得到外显子起始相位(Phase): {phase}")

    # 2. 获取外显子序列
    ensembl_server = "https://rest.ensembl.org"
    ext_seq = f"/sequence/id/{target_exon['id']}?content-type=application/json"
    try:
        r_seq = requests.get(ensembl_server + ext_seq)
        r_seq.raise_for_status()
        exon_seq = r_seq.json().get('seq')
    except requests.exceptions.RequestException as err:
        print(f"  > 获取外显子序列失败: {err}")
        return

    # 3. 根据链向处理序列和变异
    if strand == -1:
        exon_seq = reverse_complement(exon_seq)
        ref_allele = reverse_complement(ref_allele)
        alt_allele = reverse_complement(alt_allele)

    # 4. 引入变异
    # 需要计算变异在外显子序列中的相对位置
    variant_pos_in_genome = int(target_exon['variant_pos']) # 假设我们已经把这个值传进来了
    
    if strand == 1:
        relative_pos = variant_pos_in_genome - target_exon['start']
    else: # 反向链
        relative_pos = target_exon['end'] - variant_pos_in_genome
        
    mutated_exon_seq = exon_seq[:relative_pos] + alt_allele + exon_seq[relative_pos + len(ref_allele):]

    # 5. 根据相位进行翻译
    # 去掉序列开头不构成完整密码子的部分碱基
    ref_to_translate = exon_seq[phase:]
    mut_to_translate = mutated_exon_seq[phase:]
    ref_protein_fragment = translate_dna(ref_to_translate)
    mut_protein_fragment = translate_dna(mut_to_translate)
    
    print("\n  --- 单外显子翻译结果 ---")
    print(f"  > 原始外显子蛋白片段: {ref_protein_fragment}")
    print(f"  > 突变后外显子蛋白片段: {mut_protein_fragment}")
    
    # --- 新增：详细的蛋白质突变分析和打印 ---
    print("\n  --- 蛋白质变异分析 ---")
    if ref_protein_fragment == mut_protein_fragment:
        print("  - 变异类型: 同义突变 (Synonymous)")
        print("  - 变异结论: 蛋白质序列未发生改变。")
    else:
        # 计算氨基酸位置
        aa_pos = (relative_pos - phase) // 3
        ref_aa = ref_protein_fragment[aa_pos] if aa_pos < len(ref_protein_fragment) else '?'
        mut_aa = mut_protein_fragment[aa_pos] if aa_pos < len(mut_protein_fragment) else '?'

        # 判断详细突变类型Debug_R_SEQ
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
        print(f"  - 影响位置: 外显子翻译片段的第 {aa_pos + 1} 个氨基酸")
        print(f"  - 具体改变: {ref_aa} -> {mut_aa}")
    
        print(f"  - 原始片段: ...{ref_protein_fragment[max(0, aa_pos-10) : aa_pos]}[{ref_protein_fragment[aa_pos]}]{ref_protein_fragment[aa_pos+1: aa_pos+11]}...")
        print(f"  - 突变后片段: ...{mut_protein_fragment[max(0, aa_pos-10) : aa_pos]}[{mut_protein_fragment[aa_pos]}]{mut_protein_fragment[aa_pos+1: aa_pos+11]}...")
        
    # 新增：对结果进行打印
    print_exon_report(desc_id, gene_name, gene_id, target_exon, exon_seq, mutated_exon_seq,variant_position, relative_pos, ref_allele, alt_allele)


def process_variant_with_uniprot_workflow(gene_id: str, target_exon_info: dict, ref: str, alt: str,desc_id,gene_name,variant_position, relative_pos):
    """
    执行从基因到单外显子翻译的完整UniProt锚定工作流
    target_exon_info 应该是一个包含了外显子信息的字典，例如从 `generate_variant_context_report` 中获取的 target_exon
    """
    # 步骤1
    uniprot_id = get_uniprot_id_for_gene_debug(gene_id)
    if not uniprot_id:
        return

    # 步骤2
    protein_seq = get_protein_sequence_from_uniprot(uniprot_id)
    if not protein_seq:
        return
        
    # --- 新增打印 1: UniProt蛋白质的氨基酸序列 ---
    print(f"\n--- UniProt ID ({uniprot_id}) 对应的蛋白质序列 ---")
    # 为了方便阅读，每60个氨基酸换一行打印
    for i in range(0, len(protein_seq), 60):
        print(f"  {protein_seq[i:i+60]}")
    print("-" * 40)
    # ---------------------------------------------
    
    # 步骤3
    matched_transcript = find_transcript_by_protein_sequence(gene_id, protein_seq)
    if not matched_transcript:
        return
    
    # --- 新增打印 2: 匹配到的转录本ID和内容 ---
    # 使用json.dumps让输出的字典格式更美观
    import json
    print(f"\n--- 反向匹配到的转录本信息 ---")
    print(f"  - 转录本ID: {matched_transcript.get('id')}")
    # print(f"  - 转录本内容 (JSON格式):")
    # print(json.dumps(matched_transcript, indent=4)) # indent=4 表示缩进4个空格
    # print("-" * 40)
    # ---------------------------------------------
    
    # 步骤4
    translate_single_exon(target_exon_info, matched_transcript, ref, alt,variant_position, relative_pos,desc_id, gene_name, gene_id)

    
# =================================================================================
#  功能模块 1: 变异类型判断
# =================================================================================
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

    '''根据ID获取变异位置所在的gene和外显子碱基序列'''
    ensembl_server = "https://rest.ensembl.org"
    
    gene_name, gene_id = "N/A", "N/A"
    try:
        ext_overlap = f"/overlap/region/human/{chromosome}:{variant_position}-{variant_position}?feature=gene;content-type=application/json"
        r_overlap = requests.get(ensembl_server + ext_overlap)
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
            ext_lookup = f"/lookup/id/{gene_id}?content-type=application/json;expand=1"
            r_lookup = requests.get(ensembl_server + ext_lookup)
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
        # print("--Debug--target_exon:",target_exon)
        # print(f"外显子ID:{target_exon['id']}")
        exon_chr, exon_start, exon_end = target_exon['seq_region_name'], target_exon['start'], target_exon['end']
        relative_pos = variant_position - exon_start
        '''
        执行从基因到单外显子翻译的完整UniProt锚定工作流
        '''
        process_variant_with_uniprot_workflow(gene_id, target_exon, ref_allele, alt_allele,desc_id,gene_name,variant_position, relative_pos)
        try:
            ext_sequence = f"/sequence/region/human/{exon_chr}:{exon_start}-{exon_end}?content-type=application/json"
            r_seq = requests.get(ensembl_server + ext_sequence)
            r_seq.raise_for_status()
            sequence = r_seq.json().get('seq')
            """
            格式化打印外显子报告，并新增序列比对功能
            """
            # print_exon_report(desc_id, gene_name, gene_id, target_exon, sequence, variant_position, relative_pos, ref_allele, alt_allele)
        except requests.exceptions.RequestException as err:
            print(f"错误：获取外显子序列失败。原因: {err}")
    else:
        flank_size = 30
        ref_len = len(ref_allele)
        start_fetch = variant_position - flank_size
        end_fetch = (variant_position + ref_len - 1) + flank_size
        try:
            ext_sequence = f"/sequence/region/human/{chromosome}:{start_fetch}-{end_fetch}?content-type=application/json"
            r_seq = requests.get(ensembl_server + ext_sequence)
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

# =================================================================================
#  功能模块 3: VCF文件读取与主流程控制 (已优化)
# =================================================================================

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

# --- 执行脚本 ---
# if __name__ == "__main__":
#     vcf_file_path = "rgc_me_variant_frequencies_chrY_20231004.vcf"
    
#     # 默认只处理文件的前3条记录，您可以修改此数字
#     read_vcf_and_process(file_path=vcf_file_path, record_limit=3)

# 调试脚本
if __name__ == "__main__":
    generate_variant_context_report("24_13366334_T_C")