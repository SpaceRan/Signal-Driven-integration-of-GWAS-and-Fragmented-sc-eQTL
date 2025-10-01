import gzip
import pandas as pd
import os
import io
from tqdm import tqdm
import glob

"根据文件夹中的list文件，原地产生对应的MR交叉文件"

# ========================================# 
'''格式检查函数'''
def explore_finngen_gz(finngen_gz_path, n_lines=2, encoding='ascii'):
    """
    探索 FinnGen 格式的 .gz 文件结构。
    """
    print(f"🔍 正在探索文件: {finngen_gz_path}")
    print("=" * 80)

    BUFFER_SIZE = 64 * 1024
    line_count = 0
    header = None
    sample_lines = []

    try:
        with gzip.open(finngen_gz_path, 'rb') as f_raw:
            f = io.BufferedReader(f_raw, buffer_size=BUFFER_SIZE)
            for line_bytes in f:
                if not line_bytes.strip():
                    continue
                try:
                    line_str = line_bytes.decode(encoding).rstrip('\n\r')
                except UnicodeDecodeError:
                    print(f"⚠️  第 {line_count + 1} 行解码失败（尝试改 encoding='utf-8'）")
                    line_str = "<DECODE ERROR>"

                line_count += 1
                if line_count <= n_lines:
                    sample_lines.append(line_str)
                if line_count == 1 and line_str.startswith('#'):
                    print(f"📌 第1行以 '#' 开头 → 可能是注释或带#的表头")
                if '\t' in line_str:
                    delimiter = '\\t (制表符)'
                    parts = line_str.split('\t')
                elif ' ' in line_str:
                    delimiter = '空格'
                    parts = line_str.split(' ')
                else:
                    delimiter = '无明显分隔符'
                    parts = [line_str]
                if header is None:
                    header = parts
                    print(f"📋 推测表头 (第 {line_count} 行):")
                    for i, col in enumerate(header):
                        print(f"    [{i}] {col}")
                    print(f"    → 共 {len(header)} 列，分隔符推测为: {delimiter}")
                    print("-" * 80)

                # 如果已收集足够样本，跳出
                if line_count >= n_lines:
                    break

    except Exception as e:
        print(f"❌ 读取失败: {e}")
        return

    # 打印样本行
    print(f"📄 前 {min(n_lines, line_count)} 行内容:")
    for i, line in enumerate(sample_lines, 1):
        print(f"  {i:2d} | {line}")

def check_vcf_format(vcf_gz_path, n_lines=1):
    """检查 VCF 文件前 n_lines 行的结构
    打印关键字段用于人工核对"""
    print(f"=== 检查文件: {vcf_gz_path} ===\n")
    
    BUFFER_SIZE = 256 * 1024
    bytes_newline = b'\n'
    bytes_tab = b'\t'
    line_count = 0
    with gzip.open(vcf_gz_path, 'rb') as f:
        f = io.BufferedReader(f, buffer_size=BUFFER_SIZE)
        for line_bytes in f:
            if line_bytes.startswith(b'#'):
                if line_bytes.startswith(b'#CHROM'):
                    print(">>> 列头（参考）:")
                    print(line_bytes.decode('ascii', errors='replace').strip())
                    print("-" * 80)
                continue

            if line_count >= n_lines:
                break

            line_count += 1
            if line_bytes.endswith(bytes_newline):
                line_bytes = line_bytes[:-1]

            parts = line_bytes.split(bytes_tab)
            if len(parts) < 10:
                print(f"[行{line_count}] 字段不足10列，跳过")
                continue

            # 解码基本字段
            try:
                chrom = parts[0].decode('ascii')
                pos = parts[1].decode('ascii')
                rsid = parts[2].decode('ascii')
                ref = parts[3].decode('ascii')
                alt = parts[4].decode('ascii')
                fmt_str = parts[8].decode('ascii')
                sample_str = parts[9].decode('ascii')
            except UnicodeDecodeError as e:
                print(f"[行{line_count}] 解码错误: {e}")
                continue

            print(f"--- 行 {line_count} ---")
            print(f"RSID: {rsid}")
            print(f"REF:  {ref}")
            print(f"ALT:  {alt}")
            print(f"FORMAT: {fmt_str}")
            print(f"SAMPLE: {sample_str}")

            # 尝试解析 FORMAT 和 SAMPLE
            fmt_keys = fmt_str.split(':')
            sample_vals = sample_str.split(':')
            sample_dict = dict(zip(fmt_keys, sample_vals))

            print("解析的键值:")
            for k, v in sample_dict.items():
                print(f"  {k} = {v}")

            # 特别检查 ES 和 SE
            es_val = sample_dict.get('ES', 'N/A')
            se_val = sample_dict.get('SE', 'N/A')
            try:
                es_float = float(es_val)
                se_float = float(se_val)
                print(f"  [数值] ES = {es_float}, SE = {se_float}")
            except ValueError:
                print(f"  [警告] ES 或 SE 无法转为数值: ES={es_val}, SE={se_val}")

            print()  # 空行分隔

    print(f"=== 共检查 {line_count} 行有效数据 ===")

def check_txt_file_structure(file_path, delimiter='\t', encoding='utf-8'):
    """
    检查指定txt文件的列名以及前5行
    dict: 包含列名和前5行数据的字典
    """

    df = pd.read_csv(file_path, delimiter=delimiter, encoding=encoding, header=0,nrows=6)
    columns = df.columns.tolist()        
    head_data = df.head(5)    
    print(f"文件路径: {file_path}")
    print(f"文件形状: {df.shape}")
    print(f"列名: {columns}")
    print("\n前5行数据:")
    print(head_data)

# ========================================# 
"""新思路辅助函数，读取文件，再做决定"""

def load_vcf_data(vcf_gz_path):
    """加载单个 VCF 文件到内存"""
    vcf_dict = {}
    BUFFER_SIZE = 256 * 1024
    bytes_newline = b'\n'
    bytes_tab = b'\t'
    bytes_rs = b'rs'
    with gzip.open(vcf_gz_path, 'rb') as f:
        f = io.BufferedReader(f, buffer_size=BUFFER_SIZE)
        for line_bytes in f:
            if line_bytes.startswith(b'#'):
                continue
            if line_bytes.endswith(bytes_newline):
                line_bytes = line_bytes[:-1]
            parts = line_bytes.split(bytes_tab)
            if len(parts) < 10:
                continue
            rsid_b = parts[2]
            if not rsid_b.startswith(bytes_rs):
                continue
            try:
                rsid = rsid_b.decode('ascii')
            except UnicodeDecodeError:
                continue
            ref_b, alt_b = parts[3], parts[4]
            fmt_str_b, sample_str_b = parts[8], parts[9]
            try:
                fmt_str = fmt_str_b.decode('ascii')
                sample_str = sample_str_b.decode('ascii')
            except UnicodeDecodeError:
                continue
            fmt_keys = fmt_str.split(':')
            sample_vals = sample_str.split(':')
            if len(fmt_keys) != len(sample_vals):
                continue
            sample_dict = dict(zip(fmt_keys, sample_vals))
            try:
                es = float(sample_dict['ES'])
                se = float(sample_dict['SE'])
            except (ValueError, KeyError):
                continue
            try:
                ref = ref_b.decode('ascii')
                alt = alt_b.decode('ascii')
            except UnicodeDecodeError:
                continue
            vcf_dict[rsid] = {
                'beta_exp': es,
                'se_exp': se,
                'alt_exp': alt,
                'ref_exp': ref
            }
    return vcf_dict

def load_finngen_data(finngen_gz_path):
    finn_dict = {}
    BUFFER_SIZE = 256 * 1024
    with gzip.open(finngen_gz_path, 'rb') as f_raw:
        f = io.BufferedReader(f_raw, buffer_size=BUFFER_SIZE)
        header_line = f.readline()
        header_str = header_line.decode('ascii').lstrip('#').strip()
        header = header_str.split('\t')
        col_index = {name: idx for idx, name in enumerate(header)}
        rs_col = col_index['rsids']
        ref_col = col_index['ref']
        alt_col = col_index['alt']
        beta_col = col_index['beta']
        se_col = col_index['sebeta']
        pval_col = col_index['pval']
        bytes_tab = b'\t'
        for line_bytes in f:
            if not line_bytes.strip():
                continue
            parts = line_bytes.rstrip(b'\n').split(bytes_tab)
            if len(parts) <= max(rs_col, ref_col, alt_col, beta_col, se_col, pval_col):
                continue
            rsids_b = parts[rs_col]
            if not rsids_b:
                continue
            try:
                rsids_str = rsids_b.decode('ascii')
            except UnicodeDecodeError:
                continue
            rs_ids = rsids_str.split(',')
            try:
                ref = parts[ref_col].decode('ascii')
                alt = parts[alt_col].decode('ascii')
                beta = float(parts[beta_col].decode('ascii'))
                se = float(parts[se_col].decode('ascii'))
                pval = float(parts[pval_col].decode('ascii'))
            except (UnicodeDecodeError, ValueError, IndexError):
                continue
            for rs in rs_ids:
                finn_dict[rs] = {
                    'beta_exp': beta,
                    'se_exp': se,
                    'pval_exp': pval,
                    'alt_exp': alt,
                    'ref_exp': ref
                }
    return finn_dict

def enrich_csv_from_cache(csv_file, data_dict, output_path, data_type):
    df = pd.read_csv(csv_file)
    snp_series = df['SNP'].astype(str)
    mask = snp_series.isin(data_dict)
    matched_snps = snp_series[mask]
    if data_type == 'vcf':
        df['beta_exp'] = None
        df['se_exp'] = None
        df['alt_exp'] = None
        df['ref_exp'] = None
        if not matched_snps.empty:
            df.loc[mask, 'beta_exp'] = matched_snps.map(lambda x: data_dict[x]['beta_exp'])
            df.loc[mask, 'se_exp'] = matched_snps.map(lambda x: data_dict[x]['se_exp'])
            df.loc[mask, 'alt_exp'] = matched_snps.map(lambda x: data_dict[x]['alt_exp'])
            df.loc[mask, 'ref_exp'] = matched_snps.map(lambda x: data_dict[x]['ref_exp'])
    elif data_type == 'finngen':
        df['beta_exp'] = None
        df['se_exp'] = None
        df['pval_exp'] = None
        df['alt_exp'] = None
        df['ref_exp'] = None
        if not matched_snps.empty:
            df.loc[mask, 'beta_exp'] = matched_snps.map(lambda x: data_dict[x]['beta_exp'])
            df.loc[mask, 'se_exp'] = matched_snps.map(lambda x: data_dict[x]['se_exp'])
            df.loc[mask, 'pval_exp'] = matched_snps.map(lambda x: data_dict[x]['pval_exp'])
            df.loc[mask, 'alt_exp'] = matched_snps.map(lambda x: data_dict[x]['alt_exp'])
            df.loc[mask, 'ref_exp'] = matched_snps.map(lambda x: data_dict[x]['ref_exp'])
    df.to_csv(output_path, index=False)
    print(f"💾 已保存增强文件至: {output_path}")

#========================================#
'''主控流程函数：对于所有的vcfgz或gz文件，分别实施数据填充'''
def process_csv_folder(csv_folder, data_folder):
    csv_files = glob.glob(os.path.join(csv_folder, "*snplist.csv"))
    all_data_files = glob.glob(os.path.join(data_folder, "*"))
    vcf_files = [f for f in all_data_files if f.endswith('.vcf.gz')]
    finngen_files = [f for f in all_data_files if f.endswith('.gz') and not f.endswith('.vcf.gz')]
    from collections import defaultdict
    vcf_tasks = defaultdict(list)  # {vcf_file: [csv_file1, csv_file2, ...]}
    finngen_tasks = defaultdict(list)  # {finngen_file: [csv_file1, csv_file2, ...]}
    for vcf_file in vcf_files:
        for csv_file in csv_files:
            csv_base = os.path.splitext(os.path.basename(csv_file))[0]
            vcf_base = os.path.splitext(os.path.splitext(os.path.basename(vcf_file))[0])[0]
            expected_output = os.path.join(csv_folder, f"{csv_base}_{vcf_base}.csv")
            if not os.path.exists(expected_output):
                vcf_tasks[vcf_file].append(csv_file)
    for finngen_file in finngen_files:
        for csv_file in csv_files:
            csv_base = os.path.splitext(os.path.basename(csv_file))[0]
            finngen_base = os.path.splitext(os.path.basename(finngen_file))[0]
            expected_output = os.path.join(csv_folder, f"{csv_base}_{finngen_base}.csv")
            if not os.path.exists(expected_output):
                finngen_tasks[finngen_file].append(csv_file)

    total_tasks = sum(len(files) for files in vcf_tasks.values()) + sum(len(files) for files in finngen_tasks.values())
    print(f"🎯 总共需要处理 {total_tasks} 个任务")

    for vcf_file, csv_list in tqdm(vcf_tasks.items(), desc="处理 VCF 文件"):
        if not csv_list:
            continue
        print(f"\n📄 正在处理 VCF 文件: {os.path.basename(vcf_file)}")
        vcf_data = load_vcf_data(vcf_file)
        for csv_file in csv_list:
            csv_base = os.path.splitext(os.path.basename(csv_file))[0]
            vcf_base = os.path.splitext(os.path.splitext(os.path.basename(vcf_file))[0])[0]
            output_path = os.path.join(csv_folder, f"{csv_base}_{vcf_base}.csv")
            enrich_csv_from_cache(csv_file, vcf_data, output_path, data_type='vcf')
        del vcf_data
    for finngen_file, csv_list in tqdm(finngen_tasks.items(), desc="处理 FinnGen 文件"):
        if not csv_list:
            continue
        print(f"\n📄 正在处理 FinnGen 文件: {os.path.basename(finngen_file)}")
        finngen_data = load_finngen_data(finngen_file)
        for csv_file in csv_list:
            csv_base = os.path.splitext(os.path.basename(csv_file))[0]
            finngen_base = os.path.splitext(os.path.basename(finngen_file))[0]
            output_path = os.path.join(csv_folder, f"{csv_base}_{finngen_base}.csv")
            enrich_csv_from_cache(csv_file, finngen_data, output_path, data_type='finngen')
        del finngen_data

    print("\n✅ 所有任务处理完成！")


#========================================# 
if __name__ == "__main__":
    csv_folder = r"D:\desk\SMR\9-1Aresult\ZNEW_unit_SOP&MR"
    data_folder = r"D:\desk\SMR\9-1Aresult\GWAS_exp"
    process_csv_folder(csv_folder, data_folder)