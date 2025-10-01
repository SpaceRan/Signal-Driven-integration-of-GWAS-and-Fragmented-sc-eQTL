import os
import pandas as pd
from pathlib import Path

folder_path = r'D:\desk\SMR\9-1Aresult\ZNEW_unit_SOP&MR' 
files = list(Path(folder_path).glob("*.csv"))
file_names = {f.stem: f for f in files}

for stem in list(file_names.keys()):
    if stem.endswith('_LD_matrix_gwas_snplist'):
        base_prefix = stem.replace('_LD_matrix_gwas_snplist', '')
        base_file = file_names.get(base_prefix)
        snplist_file = file_names.get(stem)
        if not base_file or not snplist_file:
            print(f"跳过不完整配对: {base_prefix}.csv 或 {stem}.csv 缺失")
            continue
        print(f"处理配对: {base_file.name} -> 输出到 {snplist_file.name}")
        try:
            snplist_df = pd.read_csv(snplist_file)
            snp_list = snplist_df.iloc[:, 0].dropna().astype(str).tolist()
            snp_set = set(snp_list)
            base_df = pd.read_csv(base_file)
            base_snp_col = base_df.iloc[:, 0].astype(str)
            matched_df = base_df[base_snp_col.isin(snp_set)]
            found_snps = set(matched_df.iloc[:, 0].astype(str))
            missing_snps = snp_set - found_snps

            if missing_snps:
                print(f"  ❌ 错误: 有 {len(missing_snps)} 个 SNP 在 {base_file.name} 中未找到，例如: {list(missing_snps)[:3]}")
                print(f"  跳过此文件对，未覆盖 {snplist_file.name}")
                continue
            matched_df.to_csv(snplist_file, index=False)
            print(f"  ✅ 成功覆盖 {snplist_file.name}，共写入 {len(matched_df)} 行")
        except Exception as e:
            print(f"  ❌ 处理失败 {snplist_file.name}: {e}")
