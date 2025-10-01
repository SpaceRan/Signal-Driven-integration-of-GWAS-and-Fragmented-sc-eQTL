import os
import re
from datetime import datetime
import shutil

def get_base_name(filename):
    """获取文件的基础名称（去除结尾的差异部分）"""
    if filename.endswith('_LD_matrix.pkl'):
        return filename[:-14]  # 去除 '_LD_matrix.pkl'
    elif filename.endswith('.csv'):
        return filename[:-4]   # 去除 '.csv'
    else:
        return filename

def find_perfect_matches(pkl_folder, csv_folder):
    """找出完全匹配的文件对"""
    # 获取所有pkl和csv文件
    pkl_files = [f for f in os.listdir(pkl_folder) if f.endswith('_LD_matrix.pkl')]
    csv_files = [f for f in os.listdir(csv_folder) if f.endswith('.csv')]
    
    print(f"PKL文件夹中有 {len(pkl_files)} 个pkl文件")
    print(f"CSV文件夹中有 {len(csv_files)} 个csv文件")
    
    # 创建基础名称到完整文件名的映射
    pkl_base_map = {}
    csv_base_map = {}
    
    for pkl_file in pkl_files:
        base_name = get_base_name(pkl_file)
        pkl_base_map[base_name] = pkl_file
    
    for csv_file in csv_files:
        base_name = get_base_name(csv_file)
        csv_base_map[base_name] = csv_file
    
    # 查找完美匹配
    matches = []
    pkl_matched = set()
    csv_matched = set()
    
    print("\n=== 匹配检查 ===")
    
    for base_name, pkl_file in pkl_base_map.items():
        if base_name in csv_base_map:
            csv_file = csv_base_map[base_name]
            
            matches.append({
                'pkl_file': pkl_file,
                'csv_file': csv_file,
                'base_name': base_name
            })
            pkl_matched.add(pkl_file)
            csv_matched.add(csv_file)
            
            print(f"✓ 匹配成功:")
            print(f"  PKL: {pkl_file}")
            print(f"  CSV: {csv_file}")
            print(f"  基础名称: {base_name}")
            print("-" * 50)
        else:
            print(f"✗ 未找到匹配: {pkl_file}")
    
    # 检查完整性
    print(f"\n=== 完整性检查 ===")
    print(f"PKL文件总数: {len(pkl_files)}")
    print(f"已匹配PKL文件: {len(pkl_matched)}")
    
    if len(pkl_matched) == len(pkl_files):
        print("✓ 所有PKL文件都有唯一匹配")
        all_matched = True
    else:
        print("✗ 存在未匹配的PKL文件")
        all_matched = False
        # 显示未匹配的文件
        for pkl_file in pkl_files:
            if pkl_file not in pkl_matched:
                print(f"  未匹配: {pkl_file}")
    
    # 显示未被匹配的CSV文件
    unmatched_csv = [f for f in csv_files if f not in csv_matched]
    if unmatched_csv:
        print(f"\n未被匹配的CSV文件 ({len(unmatched_csv)} 个):")
        for csv_file in unmatched_csv:
            print(f"  {csv_file}")
    
    return matches, all_matched

def move_csv_files(matches, pkl_folder, csv_folder):
    """将匹配的csv文件移动到pkl文件夹"""
    moved_count = 0
    
    print("\n=== 移动文件 ===")
    for match in matches:
        csv_file = match['csv_file']
        source_path = os.path.join(csv_folder, csv_file)
        dest_path = os.path.join(pkl_folder, csv_file)
        
        try:
            if os.path.exists(source_path):
                shutil.move(source_path, dest_path)
                print(f"✓ 已移动: {csv_file}")
                moved_count += 1
            else:
                print(f"✗ 源文件不存在: {csv_file}")
        except Exception as e:
            print(f"✗ 移动失败 {csv_file}: {e}")
    
    return moved_count

def main():
    # 固定文件夹路径
    pkl_folder = r"D:\desk\study5_COPDxLC_SMR\Aresult_9_1\ZNEW_unit_pass1"
    csv_folder = r"D:\desk\study5_COPDxLC_SMR\Aresult_9_1\ZNEW_SMR_SNP_unit"
    
    # 检查文件夹是否存在
    if not os.path.exists(pkl_folder):
        print(f"错误: PKL文件夹 '{pkl_folder}' 不存在")
        return
    
    if not os.path.exists(csv_folder):
        print(f"错误: CSV文件夹 '{csv_folder}' 不存在")
        return
    
    print(f"PKL文件夹: {pkl_folder}")
    print(f"CSV文件夹: {csv_folder}")
    print("=" * 80)
    
    # 查找完美匹配
    matches, all_matched = find_perfect_matches(pkl_folder, csv_folder)
    
    print(f"\n总共找到 {len(matches)} 个完美匹配的文件对")
    
    if matches and all_matched and len(matches) > 0:
        print("✓ 所有PKL文件都有唯一匹配，通过检查！")
        
        # 询问是否移动文件
        response = input(f"\n是否将 {len(matches)} 个匹配的csv文件移动到pkl文件夹? (y/n): ").strip().lower()
        
        if response in ['y', 'yes', '是']:
            moved_count = move_csv_files(matches, pkl_folder, csv_folder)
            print(f"\n成功移动 {moved_count} 个文件")
            
            # 保存匹配结果
            result_file = os.path.join(pkl_folder, 'perfect_matches.txt')
            with open(result_file, 'w', encoding='utf-8') as f:
                f.write("完美匹配结果报告\n")
                f.write("=" * 50 + "\n")
                f.write(f"匹配时间: {datetime.now()}\n")
                f.write(f"匹配数: {len(matches)}\n\n")
                for i, match in enumerate(matches, 1):
                    f.write(f"{i}. 完美匹配:\n")
                    f.write(f"   PKL文件: {match['pkl_file']}\n")
                    f.write(f"   CSV文件: {match['csv_file']}\n")
                    f.write(f"   基础名称: {match['base_name']}\n\n")
            print(f"匹配报告已保存到: {result_file}")
        else:
            print("取消移动操作")
    elif len(matches) == 0:
        print("⚠ 没有找到匹配的文件")
    else:
        print("✗ 虽然找到了匹配，但不是所有PKL文件都有匹配")

if __name__ == "__main__":
    main()