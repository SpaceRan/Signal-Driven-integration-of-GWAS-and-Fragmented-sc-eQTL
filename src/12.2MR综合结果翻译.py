"""
MR-Dt综合显著的对象，作为root，在另外一个文件中进行寻找，并寻找如下信息
1.细胞类型
2.基因名称
3.共病名称
4.保护共病or引起共病
5.暴露gwas信息来源
6.结局gwas信息来源
7.QTL数据来源
"""

import pandas as pd
import os
import re

def extract_prefix_and_base(filename):
    """
    切除前缀后返回
    """
    prefixes = [
        "finngen_R12_C3_LUNG_NONSMALL_EXALLC_",
        "b-4954_"
    ]
    for prefix in prefixes:
        if filename.startswith(prefix):
            base_part = filename[len(prefix):]
            return prefix.rstrip('_'), base_part
    return None, filename  # 如果没有匹配前缀，保留原名

def get_exposure_gwas_source(filename):
    """
    从filename中提取“暴露gwas信息来源”：
    '_LD_matrix_gwas_snplist_' 之后、'.csv' 之前的内容
    """
    pattern = r'_LD_matrix_gwas_snplist_(.+)\.csv$'
    match = re.search(pattern, filename)
    if match:
        return match.group(1)
    return None

def get_base(filename):
    """
    从filename中提取用于匹配gene信息的base，
    """
    delimiter = "_LD_matrix_gwas_snplist_"
    return filename.split(delimiter, 1)[0]

def extract_comorbidity_name(filename):
    """
    - 去掉 .csv
    - 按 '_' 分割，推测共病
    """
    if filename.lower().endswith('.csv'):
        base = filename[:-4]
    else:
        base = filename  # 如果没有 .csv，直接使用
    parts = [part for part in base.split('_') if part]
    if not parts:
        return None
    last_word = parts[-1]
    if last_word == "EARLYANDLATER":
        return "COPD"
    return last_word

def infer_comorbidity_effect(row):
    """
    输入一个MR综合分析那里的row
    根据这个row的 IVW_beta 推导“这个locus的共病效应”
    """
    D_t = row.get('beta_ivw')
    try:
        D_t = float(D_t)
    except (TypeError, ValueError):
        return None
    if D_t > 0:
        return "Loci with Risk-increasing effect"
    elif D_t < 0:
        return "Loci with Risk-decreasing effect"
    return None

def extract_cell_type(base_part):
    """
    从base_part 中提取 cell_type：
    - 找第一个 '-' 之前的内容
    - 删除最后一个 "_单词"
    """
    if '-' not in base_part:
        return None
    before_dash = base_part.split('-', 1)[0]
    if '_' in before_dash:
        parts = before_dash.rsplit('_', 1)
        if len(parts) == 2:
            return parts[0]
    return before_dash

def extract_scqtl_source(base_part):
    """
    提取 scQTL 来源：
    - 以第一个 '-' 为锚点
    - 向前检索直到遇到 '_'，提取这部分
    - 向后检索直到遇到 '_'，提取这部分
    - 拼接两部分
    """
    if '-' not in base_part:
        return None
    dash_pos = base_part.find('-')

    before_start = base_part.rfind('_', 0, dash_pos)
    before_part = base_part[before_start + 1:dash_pos] if before_start != -1 else base_part[:dash_pos]
    
    after_end = base_part.find('_', dash_pos + 1)
    after_part = base_part[dash_pos + 1:after_end] if after_end != -1 else base_part[dash_pos + 1:]
    
    return before_part + after_part


def process_csv(input_file,  output_file): 
    df = pd.read_csv(input_file)
    result_final= []
    for _, row in df.iterrows():
        filename = row['filename']  ## 最基础的root信息
        prefix, base_part = extract_prefix_and_base(filename) ##得到初步base
        outcome_gwas = prefix   ## 得到结局GWAS
        exposure_gwas = get_exposure_gwas_source(filename)  ## 得到暴露GWAS
        Loci_comorbidity_effect = infer_comorbidity_effect(row) ## 得到loci效果
        commondisease = extract_comorbidity_name(filename) ##得到共病名称
        cell_type = extract_cell_type(base_part)
        scqtl_source = extract_scqtl_source(base_part)
        ### 得到了结局GWAS，暴露GWAS，loci效果，细胞名称，QTL数据来源，loci_gwas方向
        result_final.append({
            'commondisease':commondisease,
            'outcome_gwas_source': outcome_gwas,
            'exposure_gwas_source': exposure_gwas,
            'comorbidity_effect': Loci_comorbidity_effect,
            'cell_type': cell_type,
            'scQTL_source': scqtl_source,
            'filename': filename,
        })

    result_df = pd.DataFrame(result_final)
    result_df.to_csv(output_file, index=False)
    print(f"Processing complete. Output saved to {output_file}")

if __name__ == "__main__":
    INPUT_CSV = r"D:\desk\SMR\9-1Aresult\ZNEW_unit_SOP&MR\AAAA_MR_fliter.csv"         
    OUTPUT_CSV = r"D:\desk\SMR\9-1Aresult\ZNEW_unit_SOP&MR\AAAA_MR_fliter_translate.csv"
    process_csv(INPUT_CSV, OUTPUT_CSV)
