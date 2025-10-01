import pandas as pd
import os
import re

'''
根据每个基因计算的Dt，返回其综合效果
细胞-基因-抑癌或促癌
'''

def lookup_genes_and_functions(matched_rows):
    """
    输入多行，返回这个行的对应的gene和推测的gene_function
    """
    results = []
    for _, row in matched_rows.iterrows():
        gene = row.get('gene_QTL')
        gwas_dir = str(row.get('gwas_direction', '')).strip().lower()
        D_t_val = row.get('D_t')                
        if D_t_val is None or pd.isna(D_t_val):
            continue
        try:
            D_t = float(D_t_val)
        except (ValueError, TypeError):
            continue        
        if not (D_t > 0.8 or D_t < 0.2):
            continue
        gene_func = None
        if D_t > 0.8:
            gene_func = "tumor suppressor gene"
        elif D_t < 0.2:
            gene_func = "oncogene"                
        if gene and gene_func:
            results.append((gene, gene_func))
    return results


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

def get_base(filename):
    """
    从filename中提取用于匹配gene信息的base，
    """
    delimiter = "_LD_matrix_gwas_snplist_"
    return filename.split(delimiter, 1)[0]


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


def lookup_genes_and_functions(matched_rows):
    gene = matched_rows.get('gene_QTL')
    D_t_val = matched_rows.get('D_t')
    if D_t_val is None or pd.isna(D_t_val):
        return None
    if D_t_val > 0.8:
        return gene, "tumor suppressor gene"
    elif D_t_val < 0.2:
        return gene, "oncogene"
    else:
        return None            

def process_csv(input_file, output_file): 
    df = pd.read_csv(input_file)
    result_final= []
    for _, row in df.iterrows():
        filename = row['filename']  
        prefix, base_part = extract_prefix_and_base(filename) 
        outcome_gwas = prefix 
        cell_type = extract_cell_type(base_part)
        resultA = lookup_genes_and_functions(row) 
        if resultA is not None:
            gene, gene_type = resultA
            result_final.append({
                'outcome_gwas_source': outcome_gwas,
                'gene': gene,
                'gene_function': gene_type,
                'cell_type': cell_type,
                'filename': filename
                })
            
    result_df = pd.DataFrame(result_final)
    result_df.to_csv(output_file, index=False)
    print(f"Processing complete. Output saved to {output_file}")

if __name__ == "__main__":
    INPUT_CSV = r"D:\desk\study5_COPDxLC_SMR\Aresult_9_1\ZNEW_unit_SOP&MR\AAAA_gene_Dt.csv"         
    OUTPUT_CSV = r"D:\desk\study5_COPDxLC_SMR\Aresult_9_1\ZNEW_unit_SOP&MR\AAAA_gene_Dt_cellfunction.csv"
    process_csv(INPUT_CSV, OUTPUT_CSV)
