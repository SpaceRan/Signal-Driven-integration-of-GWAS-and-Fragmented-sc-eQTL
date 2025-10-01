import os
import pandas as pd
import pickle

def process_pkl_files(input_directory):
    """
    处理目录中的所有pkl文件，提取gwas_bootstrap和qtl_bootstrap中的stable_snp_id
    并保存为新的csv文件
    """
    # 遍历指定目录中的所有pkl文件
    for filename in os.listdir(input_directory):
        if filename.endswith('.pkl'):
            file_path = os.path.join(input_directory, filename)
            print(f"正在处理文件: {filename}")
            
            try:
                with open(file_path, 'rb') as f:
                    data = pickle.load(f)                
                if 'gwas_bootstrap' not in data or 'qtl_bootstrap' not in data:
                    print(f"警告: 文件 {filename} 缺少必要的键")
                    continue
                
                if 'stable_snp_id' in data['gwas_bootstrap']:
                    gwas_stable_snp_ids = data['gwas_bootstrap']['stable_snp_id']
                    gwas_df = pd.DataFrame({
                        'gwas_bootstrap_stable': gwas_stable_snp_ids
                    })
                    gwas_output_filename = filename.replace('.pkl', '_gwas.pkl')
                    gwas_output_path = os.path.join(input_directory, gwas_output_filename)
                    gwas_df.to_pickle(gwas_output_path)
                    print(f"已保存GWAS文件: {gwas_output_filename}")
                else:
                    print(f"警告: {filename} 中gwas_bootstrap缺少stable_snp_id")
                
                # 提取qtl_bootstrap中的stable_snp_id
                if 'stable_snp_id' in data['qtl_bootstrap']:
                    qtl_stable_snp_ids = data['qtl_bootstrap']['stable_snp_id']
                    # 创建DataFrame
                    qtl_df = pd.DataFrame({
                        'qtl_bootstrap_stable': qtl_stable_snp_ids
                    })
                    # 保存为CSV文件
                    qtl_output_filename = filename.replace('.pkl', '_qtl.pkl')
                    qtl_output_path = os.path.join(input_directory, qtl_output_filename)
                    qtl_df.to_pickle(qtl_output_path)
                    print(f"已保存QTL文件: {qtl_output_filename}")
                else:
                    print(f"警告: {filename} 中qtl_bootstrap缺少stable_snp_id")
                    
            except Exception as e:
                print(f"处理文件 {filename} 时出错: {str(e)}")
                continue

def process_pkl_files_as_csv(input_directory):
    """
    处理目录中的所有pkl文件，提取stable_snp_id并保存为CSV格式
    """
    # 遍历指定目录中的所有pkl文件
    for filename in os.listdir(input_directory):
        if filename.endswith('.pkl'):
            file_path = os.path.join(input_directory, filename)
            print(f"正在处理文件: {filename}")

            try:
                # 读取pkl文件
                with open(file_path, 'rb') as f:
                    data = pickle.load(f)
                
                # 检查必要的键是否存在
                if 'gwas_bootstrap' not in data or 'qtl_bootstrap' not in data:
                    print(f"警告: 文件 {filename} 缺少必要的键")
                    continue
                
                # 提取gwas_bootstrap中的stable_snp_id
                if 'stable_snp_id' in data['gwas_bootstrap']:
                    gwas_stable_snp_ids = data['gwas_bootstrap']['stable_snp_id']
                    # 创建DataFrame
                    gwas_df = pd.DataFrame({
                        'gwas_bootstrap_stable': gwas_stable_snp_ids
                    })
                    # 保存为CSV文件
                    gwas_output_filename = filename.replace('.pkl', '_gwas.csv')
                    gwas_output_path = os.path.join(input_directory, gwas_output_filename)
                    gwas_df.to_csv(gwas_output_path, index=False)
                    print(f"已保存GWAS文件: {gwas_output_filename}")
                else:
                    print(f"警告: {filename} 中gwas_bootstrap缺少stable_snp_id")
                
                # 提取qtl_bootstrap中的stable_snp_id
                if 'stable_snp_id' in data['qtl_bootstrap']:
                    qtl_stable_snp_ids = data['qtl_bootstrap']['stable_snp_id']
                    # 创建DataFrame
                    qtl_df = pd.DataFrame({
                        'qtl_bootstrap_stable': qtl_stable_snp_ids
                    })
                    # 保存为CSV文件
                    qtl_output_filename = filename.replace('.pkl', '_qtl.csv')
                    qtl_output_path = os.path.join(input_directory, qtl_output_filename)
                    qtl_df.to_csv(qtl_output_path, index=False)
                    print(f"已保存QTL文件: {qtl_output_filename}")
                else:
                    print(f"警告: {filename} 中qtl_bootstrap缺少stable_snp_id")
                    
            except Exception as e:
                print(f"处理文件 {filename} 时出错: {str(e)}")
                continue

# 使用示例
if __name__ == "__main__":
    # 设置您的pkl文件所在目录
    input_directory = r"D:\desk\study5_COPDxLC_SMR\Aresult_9_1\ZNEW_unit_pass1"  # 请修改为您的实际目录路径
    if not os.path.exists(input_directory):
        print(f"目录 {input_directory} 不存在，请检查路径")
    else:
        process_pkl_files_as_csv(input_directory)
        
        print("所有文件处理完成！")