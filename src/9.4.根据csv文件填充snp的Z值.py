from pathlib import Path
import pandas as pd
import numpy as np

## 寻找配对文件
def find_matching_files():
    '''
    matched_pairs 是这样的结构
    [
        (Path('trait1_gwas.csv'), Path('trait1.csv')),
        (Path('trait2_qtl.csv'), Path('trait2.csv')),
    ]
    '''
    directory = r"D:\desk\study5_COPDxLC_SMR\Aresult_9_1\ZNEW_unit_pass1"
    dir_path = Path(directory)
    if not dir_path.exists():
        print(f"目录 {directory} 不存在")
        return []
    csv_son_files = []
    csv_mom_files = []
    for file_path in dir_path.iterdir():
        if file_path.is_file():
            if file_path.name.endswith('_LD_matrix_gwas.csv') or file_path.name.endswith('_LD_matrix_qtl.csv'):
                csv_son_files.append(file_path)
            elif '_LD_matrix_' not in file_path.name and file_path.name.endswith('.csv'): 
                csv_mom_files.append(file_path)
    
    matched_pairs = []
    for csv_file in csv_son_files:
        csv_name = csv_file.name
        if csv_name.endswith('_gwas.csv'):
            csv_base = (csv_name.replace('_LD_matrix_gwas.csv', ''))
        elif csv_name.endswith('_qtl.csv'):
            csv_base = (csv_name.replace('_LD_matrix_qtl.csv', ''))
        else:
            continue
        for mom_file in csv_mom_files:
            mom_name = mom_file.name
            mom_base = mom_name[:-4] 
            if csv_base == mom_base:
                matched_pairs.append((csv_file, mom_file))
                break
    return matched_pairs


def lala_to_csv(csv_file, mom_file):
    if csv_file.name.endswith('_LD_matrix_gwas.csv'):
        file_type = 'gwas'
    elif csv_file.name.endswith('_LD_matrix_qtl.csv'):
        file_type = 'qtl'
    
    csv_df = pd.read_csv(csv_file)
    mom_df = pd.read_csv(mom_file)     
    required_columns = ['SNP']

    if file_type == 'gwas':
        required_columns.extend(['SE_GWAS', 'BETA_GWAS', 'P_GWAS'])
    else:  # qtl
        required_columns.extend(['SE_QTL', 'beta_QTL', 'p-value_QTL'])

    missing_columns = [col for col in required_columns if col not in mom_df.columns]
    if missing_columns:
        raise ValueError(f"mom_file缺少必要列: {missing_columns}")
    
    if file_type == 'gwas':
        mom_df = mom_df.loc[mom_df.groupby('SNP')['P_GWAS'].idxmin()]
        mom_df['Z_gwas'] = mom_df['BETA_GWAS'] / mom_df['SE_GWAS']
        z_column = 'Z_gwas'
    else:  
        mom_df = mom_df.loc[mom_df.groupby('SNP')['p-value_QTL'].idxmin()]
        mom_df['Z_qtl'] = mom_df['beta_QTL'] / mom_df['SE_QTL']
        z_column = 'Z_qtl'

    snp_z_dict = dict(zip(mom_df['SNP'], mom_df[z_column]))  ## 组合成小文件
    first_col = csv_df.columns[0]
    csv_df['Z_SNP'] = np.nan
    for idx, snp in enumerate(csv_df[first_col]):
        if snp in snp_z_dict and not pd.isna(snp_z_dict[snp]):  ## 只要是，且非控制就添加进去
            csv_df.loc[idx, 'Z_SNP'] = snp_z_dict[snp]
    csv_df.to_csv(csv_file, index=False)

if __name__ == "__main__":
    pairs = find_matching_files()    
    for i, (csv_file, mom_file) in enumerate(pairs, 1):
        lala_to_csv(csv_file, mom_file)
