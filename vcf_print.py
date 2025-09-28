import vcfpy
import requests
import json

# =================================================================================
#  Biopython 库导入 (无变动)
# =================================================================================
try:
    from Bio.Seq import Seq
    from Bio.Data import CodonTable
except ImportError:
    print("错误：缺少 'biopython' 库。请先通过 'pip install biopython' 命令进行安装。")
    exit()

# =================================================================================
#  功能模块 1 & 2 (无变动)
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

def translate_dna_to_protein(dna_sequence: str) -> (str, str):
    """将DNA序列翻译成蛋白质序列。"""
    warning = ""
    if len(dna_sequence) % 3 != 0:
        remainder = len(dna_sequence) % 3
        warning = f"注意：CDS序列长度({len(dna_sequence)})不是3的倍数，这可能由一个移码突变导致。"
        dna_sequence = dna_sequence[:-remainder]
    if not dna_sequence:
        return "", warning
    try:
        dna_seq_obj = Seq(dna_sequence)
        protein_sequence = dna_seq_obj.translate(to_stop=True) 
        return str(protein_sequence), warning
    except CodonTable.TranslationError as e:
        return f"翻译错误: {e}", warning

# =================================================================================
#  核心功能: (已修复错误处理逻辑)
# =================================================================================
ENSEMBL_SERVER = "https://rest.ensembl.org"

def generate_variant_cds_report(desc_id: str, variant_type: str):
    """
    根据正确的CDS信息，生成变异位点的翻译报告。
    """
    print("\n--- 开始处理该位点的CDS翻译分析 ---")
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

    # 1. & 2. 查找基因和标准转录本
    try:
        ext_overlap = f"/overlap/region/human/{chromosome}:{variant_position}-{variant_position}?feature=gene;content-type=application/json"
        r_overlap = requests.get(ENSEMBL_SERVER + ext_overlap, timeout=15)
        r_overlap.raise_for_status()
        overlap_data = r_overlap.json()
        if not overlap_data:
            print(f"信息：变异位点 chr{chromosome}:{variant_position} 不在任何已知基因内，判定为基因间区。")
            return
        gene_info = overlap_data[0]
        gene_id = gene_info.get('id')
        gene_name = gene_info.get('external_name', 'N/A')

        ext_lookup = f"/lookup/id/{gene_id}?content-type=application/json;expand=1"
        r_lookup = requests.get(ENSEMBL_SERVER + ext_lookup, timeout=15)
        r_lookup.raise_for_status()
        gene_data = r_lookup.json()
        transcript = next((t for t in gene_data.get('Transcript', []) if t.get('is_canonical')), None)
        if not transcript:
            print(f"错误：在基因 {gene_name} ({gene_id}) 中未找到标准转录本。")
            return
        transcript_id = transcript['id']
    except requests.exceptions.RequestException as e:
        print(f"错误：查找基因或转录本信息失败。原因: {e}")
        return

    # 3. 获取该转录本的完整CDS序列
    try:
        ext_cds = f"/sequence/id/{transcript_id}?type=cds;content-type=application/json"
        r_cds = requests.get(ENSEMBL_SERVER + ext_cds, timeout=15)
        r_cds.raise_for_status()
        reference_cds = r_cds.json().get('seq')
        if not reference_cds:
            print(f"信息：转录本 {transcript_id} 没有对应的CDS序列 (可能为非编码RNA)。")
            return
    except requests.exceptions.RequestException as e:
        print(f"错误：获取转录本 {transcript_id} 的CDS序列失败。原因: {e}")
        return

    # 4. 【已修复】使用官方API映射坐标，并正确处理API返回的400错误
    variant_start = variant_position
    variant_end = variant_position + len(ref_allele) - 1
    cds_pos = -1
    
    try:
        ext_map = f"/map/cdna/{transcript_id}/{chromosome}:{variant_start}..{variant_end}?content-type=application/json"
        r_map = requests.get(ENSEMBL_SERVER + ext_map, timeout=15)
        r_map.raise_for_status() # 对4xx或5xx错误码，主动抛出异常
        map_data = r_map.json()
        
        if not map_data or not map_data.get('mappings'):
            # 这一步理论上不会走到，因为无效映射会直接触发400错误
            print(f"信息：变异位点 chr{chromosome}:{variant_position} 无法映射到转录本 {transcript_id} 的CDS区域。")
            return
            
        mapping = map_data['mappings'][0]
        cds_pos = mapping['start'] - 1 # 1-based to 0-based

    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 400:
            # 这是我们预期的“失败”，意味着坐标不在CDS内
            error_details = e.response.json().get('error', e.response.text)
            print(f"信息：该位点无法映射到CDS。这通常意味着它位于内含子、UTR或基因区域之外。")
            print(f"      Ensembl API 返回信息: {error_details}")
        else:
            # 其他网络错误（如500服务器内部错误）
            print(f"错误：坐标映射时遇到网络错误。原因: {e}")
        return
    except requests.exceptions.RequestException as e:
        # 其他连接问题（如超时）
        print(f"错误：坐标映射时遇到网络连接问题。原因: {e}")
        return

    # 5. 构建突变CDS并进行翻译
    ref_in_cds = reference_cds[cds_pos : cds_pos + len(ref_allele)]
    if ref_in_cds.upper() != ref_allele.upper():
        print(f"严重错误：REF碱基不匹配！基因组REF '{ref_allele}' 与CDS中对应位置的碱基 '{ref_in_cds}' 不符。")
        return

    mutated_cds = reference_cds[:cds_pos] + alt_allele + reference_cds[cds_pos + len(ref_allele):]
    
    # 6. 打印最终报告
    print_cds_translation_report(desc_id, gene_name, transcript_id, reference_cds, mutated_cds, variant_type)

def print_cds_translation_report(desc_id, gene_name, transcript_id, ref_cds, mut_cds, variant_type):
    """格式化打印基于CDS的翻译报告。"""
    print(f"\n======================================================================")
    print(f"      CDS 翻译分析报告")
    print(f"      查询ID: {desc_id}")
    print(f"      所属基因: {gene_name} (转录本: {transcript_id})")
    print(f"======================================================================")
    print("\n【CDS序列比对】")
    print(f"  - 突变类型: {variant_type}")
    print(f"  - 原始CDS序列 (Reference CDS): (长度: {len(ref_cds)})")
    print(f"    {ref_cds[:60]}...")
    print(f"  - 突变后CDS序列 (Mutated CDS): (长度: {len(mut_cds)})")
    print(f"    {mut_cds[:60]}...")
    print("\n【蛋白质序列翻译】")
    ref_protein, ref_warn = translate_dna_to_protein(ref_cds)
    mut_protein, mut_warn = translate_dna_to_protein(mut_cds)
    if ref_warn: print(f"  - {ref_warn}")
    print(f"  - 原始蛋白序列 (Reference Protein): (长度: {len(ref_protein)})")
    print(f"    {ref_protein}")
    if mut_warn: print(f"  - {mut_warn}")
    print(f"  - 突变后蛋白序列 (Mutated Protein): (长度: {len(mut_protein)})")
    print(f"    {mut_protein}")
    print("\n【翻译结果分析】")
    if ref_protein == mut_protein:
        print("  - 结论: 蛋白质序列未发生改变 (同义突变)。")
    else:
        print("  - 结论: 蛋白质序列发生改变 (可能为错义、无义或移码突变)。")
    print(f"\n======================================================================")
    print("报告结束。")

# =================================================================================
#  主流程控制 (无变动)
# =================================================================================
def read_vcf_and_process(file_path: str, record_limit: int = 5):
    """读取VCF文件，调用CDS序列上下文分析。"""
    try:
        reader = vcfpy.Reader.from_path(file_path)
    except FileNotFoundError:
        print(f"错误：VCF文件 '{file_path}' 未找到。")
        return

    print(f"\n---VCF数据记录（前{record_limit}条）---")
    for i, record in enumerate(reader):
        if record_limit is not None and i >= record_limit:
            break
        
        alt_str = record.ALT[0].value if record.ALT else ""
        variant_type = classify_variant(record.REF, alt_str)

        print(f"变异位点: {record.CHROM}:{record.POS}")
        print(f"ID: {record.ID}")

        if record.ID:
            for single_id in record.ID:
                if '_' in single_id and len(single_id.split('_')) == 4:
                    generate_variant_cds_report(single_id, variant_type)
                else:
                    print(f"\nID '{single_id}' 不符合 'CHROM_POS_REF_ALT' 格式，已跳过序列分析。")
        else:
            print("\n该记录ID字段为空，跳过序列上下文分析。")
        
        print("\n" + "-"*20)

# --- 执行脚本 ---
if __name__ == "__main__":
    vcf_file_path = "rgc_me_variant_frequencies_chrY_20231004.vcf"
    read_vcf_and_process(file_path=vcf_file_path, record_limit=3)