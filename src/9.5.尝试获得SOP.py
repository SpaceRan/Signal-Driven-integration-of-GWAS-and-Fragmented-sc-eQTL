'''
我们现在需要原始parquet然后再做打算
反正r-clean是从pkl中直接获得的，我认为到时候再传入r-clean也可以
'''
import os
import pickle
import pandas as pd
import numpy as np
from pathlib import Path

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
    required_cols = {'REF_GWAS', 'ALF_GWAS', 'REF_QTL', 'ALT_QTL', 'beta_QTL'}
    if not required_cols.issubset(df.columns):
        missing = required_cols - set(df.columns)
        raise ValueError(f"缺少必需列: {missing}")
    ref_qtl = df['REF_QTL'].values      # 大写
    alt_qtl = df['ALT_QTL'].values      # 大写
    ref_gwas = df['REF_GWAS'].values    # 大写
    alt_gwas = df['ALF_GWAS'].values    # 根据你的列名是 ALF_GWAS
    beta_qtl = df['beta_QTL'].values    # 小写保持不变
    cond1 = is_subset_np(alt_gwas, alt_qtl) & is_subset_np(ref_gwas, ref_qtl)  # 方向一致
    cond2 = is_subset_np(alt_gwas, ref_qtl) & is_subset_np(ref_gwas, alt_qtl)  # 方向相反
    valid_mask = cond1 | cond2
    invalid_count = (~valid_mask).sum()
    
    print(f"共 {invalid_count} 行被丢弃（无法归类），模式和原始计算一模一样")
    adjusted_beta = np.where(cond2, -beta_qtl, beta_qtl)
    df.loc[valid_mask, 'beta_QTL'] = adjusted_beta[valid_mask]
    return df[valid_mask].copy()
###########################

def estimate_sigma_ire(z_cond, tol=1e-2, max_iter=10):
    z = np.asarray(z_cond).flatten()
    sigma2_initial =np.median(z**2) / 0.454936448
    sigma2 = sigma2_initial
    if len(z) == 0:
        return 1.0
    for iter_idx in range(max_iter):
        w = np.exp(-z**2 / (2 * sigma2 + 1e-8))        
        z2 = z**2
        sorted_idx = np.argsort(z2)
        z2_sorted = z2[sorted_idx]
        w_sorted = w[sorted_idx]
        cumw = np.cumsum(w_sorted)
        total_weight = cumw[-1]
        target = 0.5 * total_weight
        weighted_median_z2 = z2_sorted[np.searchsorted(cumw, target, side='right') - 1]
        sigma2_new = weighted_median_z2 / 0.45493644
        sigma2_new = np.clip(sigma2_new, 0.8, 5.0)        
        if abs(sigma2_new - sigma2) < tol:
            sigma2 = sigma2_new
            break
        sigma2 = sigma2_new

    return sigma2

####谱截断空间
def apply_spectral_truncation(R: np.ndarray, threshold: float = None):
    eigenvals, eigenvecs = np.linalg.eigh(R)    
    idx = np.argsort(eigenvals)[::-1]
    eigenvals = eigenvals[idx]
    eigenvecs = eigenvecs[:, idx]    
    if threshold is None:
        threshold = max(0.2, 1e-6 * eigenvals[0])  # 使用合理的默认值    
    keep = eigenvals > threshold
    if not np.any(keep):
        keep[0] = True     
    U_trunc = eigenvecs[:, keep]
    Lambda_trunc = eigenvals[keep]
    return U_trunc, Lambda_trunc

### 计算z-con
def compute_conditional_z(
    z: np.ndarray,
    R: np.ndarray,
    selected_indices: list,
    U_trunc: np.ndarray,      # 初始全局谱截断特征向量
    Lambda_trunc: np.ndarray  # 初始全局谱截断特征值
):
    if len(selected_indices) == 0:
        return z.copy()
    S = selected_indices
    R_full_sub = R[:, S]      # (p, |S|) 从所有 SNP 到已选 SNP 的 LD
    z_selected = z[S]         # (|S|,) 已选 SNP 的 Z 分数
    R_sub = R[np.ix_(S, S)]   # (|S|, |S|) 已选 SNP 之间的 LD
    try:
        U_global_S = U_trunc[S, :]  # (|S|, k) 已选 SNP 在主成分空间中的表示
        R_sub_inv = U_global_S @ np.diag(1.0 / Lambda_trunc) @ U_global_S.T
        beta = R_sub_inv @ z_selected
    except Exception as e:
        print(f"compute_conditional_z 谱截断求逆失败: {e}")
        try:
            beta = np.linalg.solve(R_sub, z_selected)
        except:
            beta = np.linalg.pinv(R_sub) @ z_selected
    proj_mean = R_full_sub @ beta
    z_cond = z - proj_mean
    return z_cond

def pivot_ld_to_matrix(ld_df):
    """将三列格式的 LD 数据转换为实对称矩阵"""
    ld_df = ld_df.drop_duplicates(subset=['ID_A', 'ID_B'])
    all_snps = sorted(set(ld_df['ID_A']) | set(ld_df['ID_B']))
    snp_map = {snp: i for i, snp in enumerate(all_snps)}
    n = len(all_snps)
    matrix = np.eye(n)  # 默认对角线为1
    for _, row in ld_df.iterrows():
        i = snp_map[row['ID_A']]
        j = snp_map[row['ID_B']]
        matrix[i, j] = matrix[j, i] = row['R']
    return pd.DataFrame(matrix, index=all_snps, columns=all_snps)


def gain_value(directory_path):
    directory = Path(directory_path)   
    pkl_files = list(directory.glob("*.pkl"))
    success_pkl = []
    for pkl_file in pkl_files:
        print(f"处理文件: {pkl_file.name}")
        base_name = pkl_file.stem        ### 根文件名称
        LD_matrix_mom_file = pkl_file.parent / ( base_name+ '.parquet')
        ld_long = pd.read_parquet(LD_matrix_mom_file, engine='fastparquet')
        ld_df = pivot_ld_to_matrix(ld_long)

        with open(pkl_file, 'rb') as f:
            data = pickle.load(f)
        keys = list(data.keys())
        fifth_key = keys[4]            
        block_data = data[fifth_key]

        to_ana_ids =  block_data['remaining_snp_idx']  ## 第二层,拉取剩余snp_id，很正确的拉取了
        to_ana_ids_list = to_ana_ids.tolist()          ## 转成 Python list

        csv_file = pkl_file.parent / (base_name.replace('_LD_matrix', '') + '.csv')
        df = pd.read_csv(csv_file)

        ## 然后去重，然后校对方向
        df_full = df.loc[df.groupby('SNP')['P_GWAS'].idxmin()]      ## 这里去重
        df_adjusted = classify_and_adjust_beta_vectorized(df_full)  ## 然后校对方向
        df_adjusted = df_adjusted.set_index('SNP') ##设置索引
        
        ## 然后要根据commsnp再进行一波调整和筛选
        ## 总之是需要把原始的parquet拉过来才行
        snps_common = ld_df.index.intersection(ld_df.columns).intersection(df_adjusted.index)
        df_sub = df_adjusted.loc[snps_common]   # 然后这个就是完全处理好了的真实df索引列表
        
        ## 处理后的内容原样进入，检索相关内容
        ## 这个是掩码后的，理论上没有问题，就是拉取了相关小局部内容
        selected_df = df_sub.iloc[to_ana_ids_list] 
        
        # 掩码+调整后的相关内容
        # Z 分数计算
        selected_df = selected_df.copy()  # 明确创建独立副本
        selected_df['z_gwas'] = selected_df['BETA_GWAS'] / selected_df['SE_GWAS']
        selected_df['z_qtl']  = selected_df['beta_QTL']  / selected_df['SE_QTL']
        # z 分数计算完成

        # 获得block的计算信息，用于传入空间进行计算
        blocks = block_data['blocks']  ## 第二层，但是需要进去再构造
        block_info_list = []

        for _, block_dict in enumerate(blocks):
            snp_ids = block_dict.get('snp_ids', [])
            first_snp = snp_ids[0]  # 第一个 SNP ID
            block_id = f"block|{first_snp}"
            z_gwas_block = block_dict.get('z_gwas_block', None)
            z_qtl_block = block_dict.get('z_qtl_block', None)
        
        ### 这里确实是找到了block内容的
            block_info_list.append({
                'block_id': block_id,
                'z_gwas_block': z_gwas_block,
                'z_qtl_block': z_qtl_block,
            })    

        block_df = pd.DataFrame(block_info_list)

        ### 然后就开始传入了，这里的逻辑应该是掩码后的
        z_all_gwas = pd.concat([
                selected_df['z_gwas'],
                block_df['z_gwas_block']
            ], ignore_index=True)     ##这个是掩码后所有准备计算的原始能量list

        z_all_qtl = pd.concat([
            selected_df['z_qtl'],
            block_df['z_qtl_block']
        ], ignore_index=True)         ##这个是掩码后所有准备计算的原始能量list

        list_snporblock = pd.concat([
            selected_df.index.to_series().reset_index(drop=True),  # ← 把索引转成 Series
            block_df['block_id']
        ], ignore_index=True)
        ### 这个是最终的内容

        ## 下面是获取相关的parquet文件信息，获得纯粹的矩阵内容
        parquet  = pkl_file.parent /(base_name + '_r_clean.parquet') #这个是很纯粹的B+P格式的LD矩阵
        df_parquet = pd.read_parquet(parquet)
        LD_matrix_extend = df_parquet.values.astype(np.float32)  # 或 np.float64

        ##由于这个LD内容是P+B模式的，因此还需要基于这个逻辑做清洗
        ## 拉取目标索引的内容，以及最后block个内容
        N = len(block_df)
        total_size = LD_matrix_extend.shape[0]  # 假设是方阵        
        selected_ids_part1 = to_ana_ids_list 
        last_N_indices = list(range(total_size - N, total_size))  # e.g. 如果 N=2, total_size=100 → [98, 99]
        combined_ids = list(dict.fromkeys(selected_ids_part1 + last_N_indices))
        R_clean = LD_matrix_extend[np.ix_(combined_ids, combined_ids)]
        ## 终于获得了拼接的R-clean

        ### 开始获得stable_id，是列名
        qtl_target_csv = pkl_file.parent / (base_name + '_qtl.csv')
        gwas_target_csv = pkl_file.parent / (base_name + '_gwas.csv')
        qtl_first_col = pd.read_csv(qtl_target_csv, usecols=[0]).iloc[:, 0]
        gwas_first_col = pd.read_csv(gwas_target_csv, usecols=[0]).iloc[:, 0]
        ## 这里是核对正确的，就是可以获得

        ## 知道stable_id分别在全部整体中的index
        mapper = pd.Series(list_snporblock.index, index=list_snporblock.values).drop_duplicates()
        index_in_list_gwas = gwas_first_col.map(mapper)
        index_in_list_qtl = qtl_first_col.map(mapper)
        ### 我们需要的是“小局部中的索引”，似乎这里是对的

        # 检查有多少没找到，如果不出内容就说明是真的
        if index_in_list_gwas.isna().any() or index_in_list_qtl.isna().any():
            missing_in_gwas = gwas_first_col[index_in_list_gwas.isna()]
            missing_in_qtl = qtl_first_col[index_in_list_qtl.isna()]
            print("GWAS 中未找到的数量:", index_in_list_gwas.isna().sum())
            if len(missing_in_gwas) > 0:
                print("GWAS 中未找到的 SNP/block:", list(missing_in_gwas))
            print("QTL 中未找到的数量:", index_in_list_qtl.isna().sum())
            if len(missing_in_qtl) > 0:
                print("QTL 中未找到的 SNP/block:", list(missing_in_qtl))
            raise KeyError("某些 SNP 或 block 在 GWAS 或 QTL 数据中未找到，说明索引不对")
        # 开始正式进行计算

        ## 求谱截断，再还原新的LD空间
        U_trunc, Lambda_trunc = apply_spectral_truncation(R_clean)
        P = U_trunc @ np.diag(1.0 / Lambda_trunc) @ U_trunc.T   
        
        ## gwas系列计算，qtl的能量在gwas中是否合理
        ## 首先估计sigma，然后传入sigma后的内容，然后能量计算
        sigma_gwas =  estimate_sigma_ire(z_all_gwas, tol=1e-2, max_iter=10)
        z_gwas_cond = z_all_gwas / np.sqrt(sigma_gwas)
        z_cond_gwas = compute_conditional_z(
                            z_gwas_cond,   ##传入的是校正的z
                            R_clean,
                            index_in_list_qtl,  ## 有趣的地方来了
                            U_trunc,            ## 全局谱截断特征向量
                            Lambda_trunc)
        
        ### 得到去除了相关内容的剩下的z值，去除已选 SNP 的 LD 贡献后剩下的独立信号
        E_total_gwas = z_gwas_cond @ P @ z_gwas_cond
        E_residual_gwas = z_cond_gwas @ P @ z_cond_gwas
        E_explained_gwas= E_total_gwas - E_residual_gwas
        Eg= E_explained_gwas/E_total_gwas

        ### qtl系列计算，gwas的能量在qtl中是否合理
        sigma_qtl =   estimate_sigma_ire(z_all_qtl , tol=1e-2, max_iter=10)
        z_qtl_cond = z_all_qtl / np.sqrt(sigma_qtl)
        z_cond_qtl = compute_conditional_z(
                    z_qtl_cond,
                    R_clean,
                    index_in_list_gwas, ## 有趣的地方来了，我传入另外一个的稳定内核
                    U_trunc,            # 全局谱截断特征向量
                    Lambda_trunc )
        E_total_qtl = z_qtl_cond @ P @ z_qtl_cond
        E_residual_qtl = z_cond_qtl @ P @ z_cond_qtl
        E_explained_qtl= E_total_qtl - E_residual_qtl
        Eq= E_explained_qtl/E_total_qtl

        SOP = (Eg+Eq)/2
        label = 1 if SOP >= 0.7 else 0
        success_pkl.append([str(pkl_file), label, SOP])  
        
    df = pd.DataFrame(success_pkl, columns=['file_path', 'state','SOP'])
    directory_path = Path(directory_path)
    df.to_csv(directory_path / 'AAA_SOP_results.csv', index=False, encoding='utf-8')

if __name__ == "__main__":
    directory_path = r"D:\desk\study5_COPDxLC_SMR\Aresult_9_1\ZNEW_unit_pass1" 
    gain_value(directory_path)