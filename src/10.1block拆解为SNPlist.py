import os
import pandas as pd
import pickle
import re


def process_pkl_csv_pairs(folder_path):
    """
    处理文件夹中的pkl和csv文件对
    """
    pkl_files = [f for f in os.listdir(folder_path) if f.endswith('.pkl')]    
    for pkl_file in pkl_files:
        csv_file = pkl_file.replace('_LD_matrix.pkl', '_LD_matrix_gwas.csv')
        pkl_path = os.path.join(folder_path, pkl_file)
        csv_path = os.path.join(folder_path, csv_file)
     
        df_csv = pd.read_csv(csv_path)
        snp_column = df_csv.iloc[:, 0]
        
        all_snps = []
        with open(pkl_path, 'rb') as f:
            pkl_data = pickle.load(f)

        block_data = pkl_data['block']
        block_mapping = {}
        
        if 'blocks' in block_data:
            blocks_list = block_data['blocks']
            for block_info in blocks_list:
                if isinstance(block_info, dict) and 'snp_ids' in block_info:
                    snp_ids = block_info['snp_ids']
                    if snp_ids and len(snp_ids) > 0:
                        first_rs = snp_ids[0]  # 直接取第一个元素
                        block_key = f"block|{first_rs}"
                        block_mapping[block_key] = snp_ids
        
        missing_blocks = []  
        for snp in snp_column:
            snp_str = str(snp).strip()
            if snp_str.startswith('rs'):
                all_snps.append(snp_str)
            elif snp_str.startswith('block|'):
                if snp_str in block_mapping:
                    all_snps.extend(block_mapping[snp_str])
                else:
                    missing_blocks.append(snp_str)
                    print(f"警告: 在PKL中找不到block: {snp_str}")
            else:
                all_snps.append(snp_str)
        if missing_blocks:
            print(f"警告: {len(missing_blocks)} 个blocks在PKL中未找到")   
            raise KeyError         
        result_df = pd.DataFrame({'SNP': all_snps})
        output_file = csv_file.replace('_gwas.csv', '_gwas_snplist.csv')
        output_path = os.path.join(folder_path, output_file)
        result_df.to_csv(output_path, index=False)
        print("-" * 50)

if __name__ == "__main__":
    folder_path = "."  # 当前目录，可以根据需要修改
    print("开始处理PKL和CSV文件对...")
    process_pkl_csv_pairs(r"D:\desk\study5_COPDxLC_SMR\Aresult_9_1\ZNEW_unit_passMR")    
    print("\n处理完成！")