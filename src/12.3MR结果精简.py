import pandas as pd
import re
import os

def extract_rs_from_filename(filename):
    # 使用正则表达式匹配rs+数字的模式
    pattern = r'rs\d{4,}'  # rs后跟至少4位数字
    match = re.search(pattern, filename)
    if match:
        return match.group(0)
    else:
        return None  # 如果没有找到匹配项，返回None

def process_and_group_csv(input_file, output_file):
    # 检查输入文件是否存在
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"输入文件不存在: {input_file}")

    df = pd.read_csv(input_file)
    
    if 'filename' not in df.columns:
        raise ValueError("CSV文件中没有找到filename列")
    
    df['Lead SNP'] = df['filename'].apply(extract_rs_from_filename)
    
    required_columns = ['commondisease', 'comorbidity_effect', 'cell_type', 'scQTL_source', 'Lead SNP']
    
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"CSV文件中缺少以下列: {missing_columns}")
    
    df_subset = df[required_columns].copy()
    
    df_subset = df_subset.dropna(subset=['Lead SNP'])
    
    group_columns = ['commondisease', 'comorbidity_effect', 'cell_type', 'scQTL_source']
    grouped = df_subset.groupby(group_columns, dropna=False).agg({
        'Lead SNP': lambda x: ','.join(sorted(set(x.astype(str))))
    }).reset_index()
    print(f"原始数据行数: {len(df_subset)}")
    print(f"分组后行数: {len(grouped)}")
    print(f"减少了 {len(df_subset) - len(grouped)} 行数据")
    grouped.to_csv(output_file, index=False)
    print(f"处理完成！汇总文件已保存为: {output_file}")
    print(f"\n汇总统计:")
    print(f"- 总组合数: {len(grouped)}")
    snp_counts = grouped['Lead SNP'].str.count(',').apply(lambda x: x + 1 if pd.notna(x) else 1)
    print(f"- 平均每个组合包含 {snp_counts.mean():.2f} 个SNP")
    print(f"- 最少SNP数: {snp_counts.min()}")
    print(f"- 最多SNP数: {snp_counts.max()}")
    
    # 显示前几行示例
    print(f"\n前5行汇总结果:")
    for i, row in grouped.head(5).iterrows():
        print(f"{i+1}. {row['commondisease']} | {row['comorbidity_effect']} | "
              f"{row['cell_type']} | {row['scQTL_source']} | {row['Lead SNP']}")

if __name__ == "__main__":
    # 请替换为你的输入和输出文件路径
    input_file_path = "D:\desk\SMR\9-1Aresult\ZNEW_unit_SOP&MR\AAAA_MR_fliter_translate.csv"  # 请修改为你的实际输入文件路径
    output_file_path = "D:\desk\SMR\9-1Aresult\ZNEW_unit_SOP&MR\AAAA_MR_fliter_translate-clean.csv"  # 输出文件路径
    try:
        process_and_group_csv(input_file_path, output_file_path)
    except FileNotFoundError:
        print(f"错误: 找不到输入文件 {input_file_path}")
        print("请确保文件路径正确")
    except ValueError as e:
        print(f"错误: {e}")
        print("请检查CSV文件是否包含所需的列")
    except Exception as e:
        print(f"处理过程中出现错误: {e}")