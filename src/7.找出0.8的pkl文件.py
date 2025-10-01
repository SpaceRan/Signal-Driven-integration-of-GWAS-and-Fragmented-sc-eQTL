import os
import shutil
import pickle

def filter_and_move_pkl_files(source_dir, target_dir, threshold=0.8):
    """
    筛选满足条件的pkl文件并移动到指定文件夹
    Args:
        source_dir (str): 源文件夹路径
        target_dir (str): 目标文件夹路径
        threshold (float): avg_completion阈值，默认0.8
    """
    # 创建目标文件夹（如果不存在）
    os.makedirs(target_dir, exist_ok=True)
    
    # 统计信息
    processed_count = 0
    moved_count = 0
    error_count = 0
    
    # 遍历源文件夹中的所有pkl文件
    for filename in os.listdir(source_dir):
        if filename.endswith('.pkl'):
            file_path = os.path.join(source_dir, filename)
            
            try:
                # 读取pkl文件
                with open(file_path, 'rb') as f:
                    data = pickle.load(f)
                
                # 检查数据是否为字典类型
                if not isinstance(data, dict):
                    print(f"警告: {filename} 不是字典格式，跳过")
                    error_count += 1
                    continue
                
                # 检查是否存在 gwas_bootstrap 键
                if 'gwas_bootstrap' not in data:
                    print(f"警告: {filename} 中未找到 'gwas_bootstrap' 键，跳过")
                    processed_count += 1
                    continue
                
                gwas_data = data['gwas_bootstrap']
                
                # 检查 gwas_bootstrap 是否为字典
                if not isinstance(gwas_data, dict):
                    print(f"警告: {filename} 中 'gwas_bootstrap' 不是字典格式，跳过")
                    processed_count += 1
                    continue
                
                # 检查是否存在 'avg_completion' 键
                if 'avg_completion' not in gwas_data:
                    print(f"警告: {filename} 中 'gwas_bootstrap' 下未找到 'avg_completion' 键，跳过")
                    processed_count += 1
                    continue
                
                # 获取 avg_completion 值
                avg_completion = gwas_data['avg_completion']
                
                # 检查是否为数字类型
                if not isinstance(avg_completion, (int, float)):
                    print(f"警告: {filename} 中 'avg_completion' 不是数字类型，跳过")
                    processed_count += 1
                    continue
                
                # 检查是否大于阈值
                if avg_completion > threshold:
                    # 移动文件
                    target_path = os.path.join(target_dir, filename)
                    shutil.move(file_path, target_path)
                    print(f"已移动: {filename} (avg_completion = {avg_completion})")
                    moved_count += 1
                else:
                    print(f"跳过: {filename} (avg_completion = {avg_completion} <= {threshold})")
                
                processed_count += 1
                
            except Exception as e:
                print(f"错误处理文件 {filename}: {str(e)}")
                error_count += 1
    
    # 打印统计信息
    print(f"\n处理完成!")
    print(f"总共处理文件数: {processed_count}")
    print(f"移动文件数: {moved_count}")
    print(f"错误文件数: {error_count}")

if __name__ == "__main__":
    source_directory = r"D:\desk\study5_COPDxLC_SMR\Aresult_9_1\ZNEW_SMR_SNP_unit"  # 源文件夹路径
    target_directory = r"D:\desk\study5_COPDxLC_SMR\Aresult_9_1\ZNEW_unit_topass"  # 目标文件夹路径
    threshold_value = 0.8
    
    # 执行筛选和移动操作
    filter_and_move_pkl_files(source_directory, target_directory, threshold_value)