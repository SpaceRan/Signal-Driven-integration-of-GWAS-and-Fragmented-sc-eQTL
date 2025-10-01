'''
思路比较简单：
1.SNP筛选F
2.进行clump
3.如果符合，则执行标准MR分析
'''
import numpy as np
import pandas as pd
import os
import scipy.stats as stats

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
    ref_exp = df['ref_exp'].values      # 注意这里要改
    alt_exp = df['alt_exp'].values      
    ref_gwas = df['REF_GWAS'].values    
    alt_gwas = df['ALF_GWAS'].values    
    beta_exp = df['beta_exp'].values
    cond1 = is_subset_np(alt_gwas, alt_exp) & is_subset_np(ref_gwas, ref_exp)  # 方向一致
    cond2 = is_subset_np(alt_gwas, ref_exp) & is_subset_np(ref_gwas, alt_exp)  # 方向相反
    valid_mask = cond1 | cond2
    invalid_count = (~valid_mask).sum()
    print(f"共 {invalid_count} 行被丢弃（无法归类），模式和原始计算一模一样")
    adjusted_beta = np.where(cond2, -beta_exp , beta_exp)
    df.loc[valid_mask, 'beta_exp'] = adjusted_beta[valid_mask]  ## 如果满足调整条件，就调整exp的内容
    return df[valid_mask].copy() ## 返回调整后的完整df

################ MR分析函数
def mr_ivw(beta_exp, se_exp, beta_out, se_out):
        beta_exp = np.array(beta_exp)
        se_exp = np.array(se_exp)
        beta_out = np.array(beta_out)
        se_out = np.array(se_out)
        n = len(beta_exp)
        if n < 2:
            return np.nan, np.nan, np.nan
        weights = 1 / (se_out ** 2)
        numerator = np.sum(weights * beta_out * beta_exp)
        denominator = np.sum(weights * beta_exp ** 2)
        if denominator <= 0:
            return np.nan, np.nan, np.nan
        beta_ivw = numerator / denominator        
        residuals = beta_out - beta_ivw * beta_exp
        Q = np.sum(weights * residuals ** 2)
        Q_pval = 1 - stats.chi2.cdf(Q, df=n-1) if Q >= 0 else np.nan
        se_ivw = np.sqrt(1 / denominator)
        z = beta_ivw / se_ivw
        IVW_pval = 2 * (1 - stats.norm.cdf(abs(z)))
        return IVW_pval, Q_pval, beta_ivw

def mr_egger(beta_exp, se_exp, beta_out, se_out):
    beta_exp = np.array(beta_exp)
    se_exp = np.array(se_exp)
    beta_out = np.array(beta_out)
    se_out = np.array(se_out)
    if not len(beta_exp) == len(se_exp) == len(beta_out) == len(se_out):
        raise ValueError("All inputs must have same length.")
    n = len(beta_exp)
    weights = 1 / (se_out ** 2)
    X = np.vstack([np.ones(n), beta_exp]).T  # [1, X] for intercept + slope
    W = np.diag(weights)
    try:
        XtW = X.T @ W
        XtWX = XtW @ X
        XtWy = XtW @ beta_out
        coef = np.linalg.solve(XtWX, XtWy)
        intercept, beta = coef
        residuals = beta_out - X @ coef
        df = n - 2
        sigma2 = np.sum(weights * residuals**2) / df
        var_cov = sigma2 * np.linalg.inv(XtWX)
        se_intercept = np.sqrt(var_cov[0, 0])
        se_beta = np.sqrt(var_cov[1, 1])
        t_intercept = intercept / se_intercept
        t_beta = beta / se_beta
        pval_intercept = 2 * (1 - stats.t.cdf(abs(t_intercept), df=df))
        pval_beta = 2 * (1 - stats.t.cdf(abs(t_beta), df=df))        
        return pval_beta, pval_intercept
    
    except np.linalg.LinAlgError:
        return np.nan, np.nan

################ 留一法敏感性分析
'''这个敏感性分析函数先放在这里，但是我不准备用'''
def mr_leave_one_out(method_func, beta_exp, se_exp, beta_out, se_out):
    n = len(beta_exp)
    results = []
    for i in range(n):
        idx = [j for j in range(n) if j != i]
        res = method_func(
            [beta_exp[j] for j in idx],
            [se_exp[j] for j in idx],
            [beta_out[j] for j in idx],
            [se_out[j] for j in idx]
        )
        res['left_out_snp_index'] = i
        results.append(res)
    return results

################ 第一波过滤函数，数量过滤、调整方向、并且计算F
def process_all_csv_in_folder(folder_path, output_summary_name="AAAA_MR.csv"):
    results = []
    csv_files = [
        f for f in os.listdir(folder_path)
        if f.endswith('.csv') and not f.endswith('snplist.csv')
    ]
    for filename in csv_files:
        filepath = os.path.join(folder_path, filename)
        print(f"\n正在处理: {filename}")
        try:
            df = pd.read_csv(filepath)
            first_column = df.columns[0]
            df = df.drop_duplicates(subset=[first_column]) ##去重
            total_rows = len(df)
            if total_rows == 0:
                results.append({'filename': filename, 'state': 'skip (empty)'})
                continue
            ref_exp_missing_rate = df['ref_exp'].isna().mean()
            if ref_exp_missing_rate > 0.3:
                print(f"  ❌ ref_exp 缺失率 {ref_exp_missing_rate:.1%} > 30%，跳过")
                results.append({'filename': filename, 'state': 'skip (missing >30%)'})
                continue
            df_adjusted = classify_and_adjust_beta_vectorized(df.copy())
            retained_ratio = len(df_adjusted) / total_rows
            if retained_ratio < 0.7:
                print(f"  ❌ 调整后保留比例 {retained_ratio:.1%} < 70%，跳过")
                results.append({'filename': filename, 'state': 'skip (retained <70%)'})
                continue

            # 提取用于 分析的 clean 数据
            beta_exp_clean = df_adjusted['beta_exp'].values
            se_exp_clean = df_adjusted['se_exp'].values
            beta_gwas_clean = df_adjusted['BETA_GWAS'].values
            se_gwas_clean = df_adjusted['SE_GWAS'].values 

            # 过滤掉任何 NaN
            mask = ~(
                np.isnan(beta_exp_clean) | 
                np.isnan(se_exp_clean) | 
                np.isnan(beta_gwas_clean) | 
                np.isnan(se_gwas_clean)
            )
            beta_exp_clean = beta_exp_clean[mask]
            se_exp_clean = se_exp_clean[mask]
            beta_gwas_clean = beta_gwas_clean[mask]
            se_gwas_clean = se_gwas_clean[mask]

            if len(beta_exp_clean) == 0:
                print("  ❌ 调整后无有效数据用于计算")
                results.append({'filename': filename, 'state': 'skip (no valid data after adjust)'})
                continue
        
            # ========== 计算平均 F 和最小 F ==========
            F_stats = (beta_exp_clean / se_exp_clean) ** 2
            mean_F = np.mean(F_stats)
            min_F = np.min(F_stats)
            # ========== 调用 MR 函数 ==========
            ivw_pval, ivw_q_pval, beta_ivw = mr_ivw(beta_exp_clean, se_exp_clean, beta_gwas_clean, se_gwas_clean)
            egger_pval, egger_intercept_pval = mr_egger(beta_exp_clean, se_exp_clean, beta_gwas_clean, se_gwas_clean)
            # ========== 记录结果 ==========
            result_entry = {
                'filename': filename,
                'state': 'success',
                'total_snps': total_rows,
                'retained_snps': len(beta_exp_clean),
                'mean_F': mean_F,
                'min_F': min_F,
                'IVW_pval': ivw_pval,
                'IVW_Q_pval': ivw_q_pval,
                "beta_ivw" :beta_ivw,
                'Egger_pval': egger_pval,
                'Egger_intercept_pval': egger_intercept_pval
            }
            results.append(result_entry)
        except Exception as e:
            print(f"  ❌ 处理失败: {e}")
            results.append({'filename': filename, 'state': f'error: {str(e)}'})
        results_df = pd.DataFrame(results)
        output_path = os.path.join(folder_path, output_summary_name)
        results_df.to_csv(output_path, index=False, encoding='utf-8-sig')

if __name__ == "__main__":
    folder = r"D:\desk\SMR\9-1Aresult\ZNEW_unit_SOP&MR" 
    process_all_csv_in_folder(folder)