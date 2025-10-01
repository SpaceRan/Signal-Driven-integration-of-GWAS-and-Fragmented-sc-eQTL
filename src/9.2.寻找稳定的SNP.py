import pickle

def analyze_pkl_file(file_path):
    with open(file_path, 'rb') as f:
        data = pickle.load(f)
    
    if not isinstance(data, dict):
        print("数据不是字典类型")
        return
    
    # 获取所有键
    keys = list(data.keys())
    print(f"字典总共有 {len(keys)} 个键")
    
    # 显示所有键的名称
    print("所有键的名称：")
    for i, key in enumerate(keys):
        print(f"  {i+1}. {key}")
    
    # 检查是否有第五个键
    if len(keys) < 5:
        print("字典中没有第五个键")
        return
    
    fifth_key = keys[4]  # 索引4对应第五个元素
    print(f"\n第五个键的名称: {fifth_key}")
    print(f"第五个键的类型: {type(data[fifth_key])}")    
    fifth_value = data[fifth_key]
    
    # 如果第五个键的值也是字典，查看它的第一个键
    if isinstance(fifth_value, dict):
        sub_keys = list(fifth_value.keys())
        print(f"第五个键对应的字典有 {len(sub_keys)} 个键")
        if len(sub_keys) > 0:
            first_sub_key = sub_keys[0]
            first_sub_value = fifth_value[first_sub_key]            
            print(f"\n第五个键的第一个键的名称: {first_sub_key}")
            print(f"第五个键的第一个键的类型: {type(first_sub_value)}")
            print(f"第五个键的第一个键的内容:")

            # 根据不同类型显示内容
            if isinstance(first_sub_value, (str, int, float, bool)):
                print(first_sub_value)
            elif isinstance(first_sub_value, (list, tuple)):
                print(f"列表/元组长度: {len(first_sub_value)}")
                if len(first_sub_value) > 0:
                    print("前5个元素:")
                    for i, item in enumerate(first_sub_value[:5]):
                        print(f"  [{i}]: {item}")
                    if len(first_sub_value) > 5:
                        print(f"  ... 还有 {len(first_sub_value) - 5} 个元素")
            elif isinstance(first_sub_value, dict):
                print("字典内容:")
                for k, v in list(first_sub_value.items())[:5]:
                    print(f"  {k}: {v}")
                if len(first_sub_value) > 5:
                    print(f"  ... 还有 {len(first_sub_value) - 5} 个键值对")
            else:
                print(first_sub_value)
        else:
            print("第五个键对应的字典为空")
    else:
        print("第五个键对应的值不是字典类型")

# 使用示例
if __name__ == "__main__":
    file_path = r'D:\desk\study5_COPDxLC_SMR\Aresult_9_1\ZNEW_unit_pass1\b-4954_CD4_Memory_after_acquiring_effector_functions_(5_d)_Soskic-2022-Natp_Genet_dedup_rs2036534_LD_matrix.pkl'
    analyze_pkl_file(file_path)