'''
单纯用于计算Dt的函数
对于一个SOP合格的locus文件/*snplist文件，计算一个locus内部每个gene内部dt： QTL--结局 的方向一致性
符合方法1.3节
输出AAAA_Dt.csv文件，每个locus每个gene都有dt（QTL--结局）
'''
import numpy as np
import pandas as pd
import os
################方向掉转函数
def is_subset_np(a_arr, b_arr):
    result = []
    for a, b in zip(a_arr, b_arr):
        if pd.isna(a) or pd.isna(b):
            result.append(False)
            continue
        str_a = str(a).upper()
        str_b = str(b).upper()
        set_a = set(str_a.split(','))
        set_b = set(str_b.split(','))
        if any(len(allele) > 1 and allele != '-' for allele in set_a | set_b):
            result.append(False)
            continue
        result.append(set_a <= set_b or set_b <= set_a)
    return np.array(result)

def classify_and_adjust_beta_vectorized(df):
    ref_exp = df['REF_QTL'].values      # 注意这里要改
    alt_exp = df['ALT_QTL'].values      
    ref_gwas = df['REF_GWAS'].values    
    alt_gwas = df['ALF_GWAS'].values    
    beta_exp = df['beta_QTL'].values   ## 被调整的对象
    cond1 = is_subset_np(alt_gwas, alt_exp) & is_subset_np(ref_gwas, ref_exp)  # 方向一致
    cond2 = is_subset_np(alt_gwas, ref_exp) & is_subset_np(ref_gwas, alt_exp)  # 方向相反
    valid_mask = cond1 | cond2
    invalid_count = (~valid_mask).sum()
    print(f"共 {invalid_count} 行被丢弃（无法归类），模式和原始计算一模一样")
    adjusted_beta = np.where(cond2, -beta_exp , beta_exp)
    df.loc[valid_mask, 'beta_QTL'] = adjusted_beta[valid_mask]  ## 如果满足调整条件，就调整exp的内容
    return df[valid_mask].copy()                                ## 返回调整后的完整df

################
def compute_direction_agreement(beta_nsclc, beta_trait):
    beta_nsclc = np.array(beta_nsclc)
    beta_trait = np.array(beta_trait)
    if len(beta_nsclc) != len(beta_trait):
        raise ValueError("must have the same length.")
    weights = np.abs(beta_nsclc)
    sign_same = (np.sign(beta_nsclc) == np.sign(beta_trait))
    numerator = np.sum(weights * sign_same)
    denominator = np.sum(weights)
    if denominator == 0:
        return 0.0
    D_t = numerator / denominator
    return D_t
################

def process_all_csv_in_folder(folder_path):
    csv_files = [
        f for f in os.listdir(folder_path)
        if f.endswith('snplist.csv') 
    ]
    results_list = []
    for filename in csv_files:
        filepath = os.path.join(folder_path, filename)
        print(f"\n正在处理: {filename}")
        df = pd.read_csv(filepath)
        try:
            df_adjusted = classify_and_adjust_beta_vectorized(df.copy())
        except:
            continue  # 直接跳过啊

        unique_gene_qtls = df_adjusted['gene_QTL'].unique()
        for gene_qtl in unique_gene_qtls:
            subset = df_adjusted[df_adjusted['gene_QTL'] == gene_qtl]
            beta_qtl_clean = subset['beta_QTL'].values
            beta_gwas_clean = subset['BETA_GWAS'].values
            D_t = compute_direction_agreement(beta_qtl_clean, beta_gwas_clean)
            results_list.append({
            'filename': filename,
            'gene_QTL': gene_qtl,
            'D_t': D_t })

    results_df = pd.DataFrame(results_list)
    results_df.to_csv(r"D:\desk\study5_COPDxLC_SMR\Aresult_9_1\ZNEW_unit_SOP&MR\AAAA_gene_Dt.csv", index=False)

#=======================
if __name__ == "__main__":
    folder_path = r"D:\desk\study5_COPDxLC_SMR\Aresult_9_1\ZNEW_unit_SOP&MR"
    process_all_csv_in_folder(folder_path)
