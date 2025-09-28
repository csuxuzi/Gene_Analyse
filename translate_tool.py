# -*- coding: utf-8 -*-

# =================================================================================
#  功能：对给定的两个DNA序列（变异前/后），从指定位置开始翻译并比较蛋白质
# =================================================================================

# 标准密码子表（包含终止密码子 '*'）
CODON_MAP = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

def translate_dna(dna_seq: str) -> str:
    """将DNA序列翻译成蛋白质序列。"""
    protein = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3].upper()
        amino_acid = CODON_MAP.get(codon, 'X') # 'X' for unknown codons
        protein.append(amino_acid)
    return "".join(protein)

def run_translation_analysis(ref_sequence: str, mut_sequence: str, translation_start_base: int):
    """
    对变异前后的两个DNA序列进行翻译分析。
    """
    print("======================================================================")
    print("                      DNA 序列翻译比对工具")
    print("======================================================================")
    print(f"\n【输入参数】")
    print(f"  - 变异前DNA序列: {ref_sequence[:60]}...")
    print(f"  - 变异后DNA序列: {mut_sequence[:60]}...")
    print(f"  - 翻译起始位置:  第 {translation_start_base} 个碱基")

    # 1. 自动检测第一个发生变异的碱基位置
    mutation_index = -1
    # 比较共同长度的部分
    for i, (ref_base, mut_base) in enumerate(zip(ref_sequence, mut_sequence)):
        if ref_base.upper() != mut_base.upper():
            mutation_index = i
            break
    # 如果共同长度部分完全一样，但总长度不同（插入/缺失），则变异点在共同长度的末尾
    if mutation_index == -1 and len(ref_sequence) != len(mut_sequence):
        mutation_index = min(len(ref_sequence), len(mut_sequence))

    # 2. 参数预处理 (将1-based转为0-based索引)
    start_index = translation_start_base - 1

    # 3. 根据指定的起始位置翻译蛋白
    ref_protein = translate_dna(ref_sequence[start_index:])
    mut_protein = translate_dna(mut_sequence[start_index:])

    # 4. 计算并标记受影响的氨基酸
    aa_pos_to_mark = -1
    if mutation_index != -1 and mutation_index >= start_index:
        aa_pos = (mutation_index - start_index) // 3
        if aa_pos < len(ref_protein):
            aa_pos_to_mark = aa_pos
    
    marked_ref_protein, marked_mut_protein = ref_protein, mut_protein
    if aa_pos_to_mark != -1:
        p_ref, idx = ref_protein, aa_pos_to_mark
        marked_ref_protein = f"{p_ref[:idx]}[{p_ref[idx]}]{p_ref[idx+1:]}"
        if idx < len(mut_protein):
            p_mut = mut_protein
            marked_mut_protein = f"{p_mut[:idx]}[{p_mut[idx]}]{p_mut[idx+1:]}"

    # 5. 打印最终结果
    print("\n【翻译结果】")
    print(f"  - 变异前的蛋白序列: {marked_ref_protein}")
    print(f"  - 变异后的蛋白序列: {marked_mut_protein}")

    # 6. 分析并打印变异影响
    effect, change = "无变化", "N/A"
    if ref_protein != mut_protein:
        # 优先判断是否为移码突变
        if len(ref_protein) != len(mut_protein):
            effect = "移码突变 (Frameshift)"
            change = "蛋白质序列长度发生改变"
        elif aa_pos_to_mark != -1: # 如果是等长替换
            ref_aa = ref_protein[aa_pos_to_mark]
            mut_aa = mut_protein[aa_pos_to_mark]
            if ref_aa != mut_aa:
                 if mut_aa == '*':
                    effect = "无义突变 (Nonsense)"
                 else:
                    effect = "错义突变 (Missense)"
                 change = f"氨基酸从 {ref_aa} 变成了 {mut_aa}"
            else:
                effect = "同义突变 (Synonymous)"
                change = f"氨基酸 {ref_aa} 未发生改变"

    print(f"\n【变异分析】")
    print(f"  - 变异的影响: {effect}")
    print(f"  - 具体改变:   {change}")
    print("\n======================================================================")


# --- 执行脚本 ---
if __name__ == "__main__":
    # --- 请在这里输入您的参数 ---
    # 1. 变异前的DNA序列
    ref_dna = "CTTCCAAGAAAATACCACGATTGGCCAGAATTAGGATCTGCCTCCAAAGACTTTTGGAGATACTGAATAGCATAGCTTTCCTTTGTGGCTTTGTCTCCTACTAGATCCATATTATGATGCATCCAAC"
    
    # 2. 变异后的DNA序列
    mut_dna = "CTTCCAAGAAAATACCACGATTGGCCAGAATTAGGATCTGCCTCCAAAGACTTTTGGAGATACGGAATAGCATAGCTTTCCTTTGTGGCTTTGTCTCCTACTAGATCCATATTATGATGCATCCAAC"

    # 3. 从第几个碱基开始翻译
    start_pos = 2
    # --- 参数输入结束 ---

    run_translation_analysis(ref_dna, mut_dna, start_pos)
    