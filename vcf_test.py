import vcfpy
import requests
from datetime import datetime
import json

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

# =================================================================================
#  功能模块 2: 序列获取与报告生成 (已优化)
# =================================================================================

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
        exon_chr, exon_start, exon_end = target_exon['seq_region_name'], target_exon['start'], target_exon['end']
        try:
            ext_sequence = f"/sequence/region/human/{exon_chr}:{exon_start}-{exon_end}?content-type=application/json"
            r_seq = requests.get(ensembl_server + ext_sequence)
            r_seq.raise_for_status()
            sequence = r_seq.json().get('seq')
            relative_pos = variant_position - exon_start
            print_exon_report(desc_id, gene_name, gene_id, target_exon, sequence, variant_position, relative_pos, ref_allele, alt_allele)
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

def print_exon_report(desc_id, gene_name, gene_id, exon, sequence, variant_pos, relative_pos, ref_allele, alt_allele):
    """格式化打印外显子报告，并新增序列比对功能"""
    print(f"\n======================================================================")
    print(f"      外显子 (Exon) 上下文序列报告")
    print(f"      查询ID: {desc_id}")
    print(f"      所属基因: {gene_name} ({gene_id})")
    print(f"======================================================================")
    print(f"\n【目标外显子信息】")
    print(f"  - 编号:       外显子 {exon['number']} (共 {exon['total_count']} 个)")
    print(f"  - 坐标:       chr{exon['seq_region_name']}:{exon['start']}-{exon['end']}")
    
    mutated_sequence = sequence[:relative_pos] + alt_allele + sequence[relative_pos + len(ref_allele):]

    print("\n【序列比对】")
    print(f"  - 原始外显子序列 (Reference):")
    for i in range(0, len(sequence), 60):
         print(f"    {sequence[i:i+60]}")
    print(f"  - 突变后外显子序列 (Mutated):")
    for i in range(0, len(mutated_sequence), 60):
         print(f"    {mutated_sequence[i:i+60]}")
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
if __name__ == "__main__":
    vcf_file_path = "rgc_me_variant_frequencies_chrY_20231004.vcf"
    
    # 默认只处理文件的前3条记录，您可以修改此数字
    read_vcf_and_process(file_path=vcf_file_path, record_limit=3)