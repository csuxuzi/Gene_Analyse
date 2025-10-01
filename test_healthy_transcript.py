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

if __name__ == "__main__":
    # 我们来验证之前匹配到的那个“健康”的转录本
    healthy_transcript_id = "ENST00000383070"
    verify_translation_logic_fixed(healthy_transcript_id)