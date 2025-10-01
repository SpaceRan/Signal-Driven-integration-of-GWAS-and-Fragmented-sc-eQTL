import pickle
import os
from pathlib import Path
import pandas as pd

## 寻找配对文件
def find_matching_files():
    '''
    matched_pairs 是这样的结构
    [
        (Path('trait1_gwas.csv'), Path('trait1_LD_matrix.pkl')),
        (Path('trait2_qtl.csv'), Path('trait2_LD_matrix.pkl')),
        (Path('geneA_gwas.csv'), Path('geneA_LD_matrix.pkl'))
    ]
    '''
    directory = r"D:\desk\study5_COPDxLC_SMR\Aresult_9_1\ZNEW_unit_pass1"
    dir_path = Path(directory)
    if not dir_path.exists():
        print(f"目录 {directory} 不存在")
        return []    
    csv_files = []
    pkl_files = []
    for file_path in dir_path.iterdir():
        if file_path.is_file():
            if file_path.name.endswith('_gwas.csv') or file_path.name.endswith('_qtl.csv'):
                csv_files.append(file_path)
            elif file_path.name.endswith('.pkl'):
                pkl_files.append(file_path)
    matched_pairs = []
    for csv_file in csv_files:
        csv_name = csv_file.name
        if csv_name.endswith('_gwas.csv'):
            csv_base = csv_name[:-9]  # 去掉'_gwas.csv'
        elif csv_name.endswith('_qtl.csv'):
            csv_base = csv_name[:-8]  # 去掉'_qtl.csv'
        else:
            continue            
        for pkl_file in pkl_files:
            pkl_name = pkl_file.name
            if pkl_name.endswith('.pkl'):
                pkl_base = pkl_name[:-4] 
                if csv_base == pkl_base:
                    matched_pairs.append((csv_file, pkl_file))
                    break
    return matched_pairs

## 输入block|rsxxxx, 拉取相关Z值
def get_snp_block_info(pkl_file_path, block_query):
    _, query_rsid = block_query.split('|', 1)
    with open(pkl_file_path, 'rb') as f:
        data = pickle.load(f)
    blocks_list = data['block']['blocks']
    for block_item in blocks_list:
        if isinstance(block_item, dict):
            if ('snp_ids' in block_item and 
                len(block_item['snp_ids']) > 0 and 
                block_item['snp_ids'][0] == query_rsid):
                return {
                    'snp_ids': block_item.get('snp_ids', []),
                    'z_gwas_block': block_item.get('z_gwas_block', None),
                    'z_qtl_block': block_item.get('z_qtl_block', None)
                }    
    return None

## 输入 csv_file_path, pkl_file_path，直接配对式获得相关信息，并直接写回源文件
def append_z_values_to_csv(csv_file_path, pkl_file_path):
    file_name = csv_file_path.name
    if file_name.endswith('_gwas.csv'):
        z_key = 'z_gwas_block'
    elif file_name.endswith('_qtl.csv'):
        z_key = 'z_qtl_block'
    df = pd.read_csv(csv_file_path)
    df['z_block_value'] = None
    if df.shape[1] == 0:
        print("⚠️ 文件无列，跳过")
        return
    first_col = df.iloc[:, 0]  # 第一列
    for idx, cell in enumerate(first_col):
        if isinstance(cell, str) and cell.startswith('block|rs'):
            block_query = cell
            block_info = get_snp_block_info(pkl_file_path, block_query)
            if block_info and z_key in block_info and block_info[z_key] is not None:
                df.at[idx, 'z_block_value'] = block_info[z_key]
    df.to_csv( csv_file_path, index=False)
    print(f"✅ 已保存: { csv_file_path}")

if __name__ == "__main__":
    pairs = find_matching_files()    
    for i, (csv_file, pkl_file) in enumerate(pairs, 1):
        append_z_values_to_csv(csv_file, pkl_file)
    print("\n🎉 所有文件处理完成！")
