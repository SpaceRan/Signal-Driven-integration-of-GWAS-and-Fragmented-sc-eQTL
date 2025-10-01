import pandas as pd
import os
'''
含义：拉取符合条件的MR结果
生成AAAA_MR_fliter文件：这个里面的都是符合严苛条件的MR分析结果
'''

def lala_printer(csv2_path):
    output_path = r"D:\desk\SMR\9-1Aresult\ZNEW_unit_SOP&MR\AAAA_MR_fliter.csv"
    output_dir = os.path.dirname(output_path)
    os.makedirs(output_dir, exist_ok=True)
    df2 = pd.read_csv(csv2_path)
    df2_filtered = df2[
        (df2["IVW_pval"] < 0.05) &                # 主效应显著
        (df2["IVW_Q_pval"] > 0.01) &              # 允许轻微显著异质性
        (df2["Egger_intercept_pval"] > 0.01)      # 允许轻微水平多效性
    ]
    
    # 保存过滤后的df2_filtered（不是df2）
    df2_filtered.to_csv(output_path, index=False)
    print(f"结果已保存到: {output_path}")
    return df2_filtered  # 返回过滤后的数据框
    

if __name__ == "__main__":
    csv2_path = r"D:\desk\SMR\9-1Aresult\ZNEW_unit_SOP&MR\AAAA_MR.csv"
    lala_printer(csv2_path)
    
