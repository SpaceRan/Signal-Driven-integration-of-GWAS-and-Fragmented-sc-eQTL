import os
import pickle
from collections import OrderedDict

def analyze_block_content(file_path):
    """
    详细分析pkl文件中[block]键的详细内容
    
    Args:
        file_path (str): pkl文件路径
    """
    try:
        # 读取pkl文件
        with open(file_path, 'rb') as f:
            data = pickle.load(f)
        
        print(f"📁 文件: {os.path.basename(file_path)}")
        print(f"📊 文件大小: {os.path.getsize(file_path)} bytes")
        print(f"📦 顶层数据类型: {type(data).__name__}")
        
        # 检查是否为字典
        if not isinstance(data, dict):
            print("❌ 顶层数据不是字典格式")
            return
        
        print(f"딕 顶层键数: {len(data)}")
        print("📋 顶层键列表:")
        for key in data.keys():
            print(f"  • {repr(key)} ({type(data[key]).__name__})")
        
        # 检查是否存在[block]键
        if 'block' not in data:
            print("\n❌ 未找到 'block' 键")
            # 显示所有顶层键的详细信息
            analyze_all_top_level_keys(data)
            return
        
        print(f"\n🎯 找到 'block' 键!")
        block_data = data['block']
        print(f"📦 block数据类型: {type(block_data).__name__}")
        
        # 分析block内容
        if isinstance(block_data, dict):
            analyze_block_dict(block_data)
        elif isinstance(block_data, list):
            analyze_block_list(block_data)
        elif isinstance(block_data, tuple):
            analyze_block_tuple(block_data)
        else:
            print(f"📝 block值: {block_data}")
            print(f"📏 block值长度: {len(str(block_data)) if hasattr(block_data, '__len__') else 'N/A'}")
            
    except Exception as e:
        print(f"❌ 错误读取文件: {str(e)}")
        import traceback
        traceback.print_exc()

def analyze_block_dict(block_dict):
    """
    详细分析block字典内容
    """
    print(f"\n딕 block字典键数: {len(block_dict)}")
    print("📋 block键详细信息:")
    
    for i, (key, value) in enumerate(block_dict.items()):
        key_type = type(key).__name__
        value_type = type(value).__name__
        
        print(f"\n  🔑 键 {i+1}: {repr(key)} [{key_type}]")
        print(f"     📦 值类型: {value_type}")
        
        if isinstance(value, dict):
            print(f"     📊 字典大小: {len(value)} 个键")
            if len(value) > 0:
                print("     📋 子键列表:")
                for sub_key in list(value.keys())[:10]:  # 只显示前10个子键
                    sub_value = value[sub_key]
                    print(f"       • {repr(sub_key)} ({type(sub_value).__name__})")
                if len(value) > 10:
                    print(f"       ... (还有{len(value)-10}个键)")
                    
        elif isinstance(value, list):
            print(f"     📊 列表长度: {len(value)}")
            if len(value) > 0:
                print("     📌 前几个元素:")
                for j, item in enumerate(value[:5]):  # 只显示前5个元素
                    print(f"       [{j}] {type(item).__name__}: {repr(item)[:100]}{'...' if len(repr(item)) > 100 else ''}")
                if len(value) > 5:
                    print(f"       ... (还有{len(value)-5}个元素)")
                    
        elif isinstance(value, tuple):
            print(f"     📊 元组长度: {len(value)}")
            if len(value) > 0:
                print("     🎯 前几个元素:")
                for j, item in enumerate(value[:5]):  # 只显示前5个元素
                    print(f"       [{j}] {type(item).__name__}: {repr(item)[:100]}{'...' if len(repr(item)) > 100 else ''}")
                if len(value) > 5:
                    print(f"       ... (还有{len(value)-5}个元素)")
        else:
            value_str = repr(value)
            if len(value_str) > 200:
                print(f"     📝 值: {value_str[:200]}...")
            else:
                print(f"     📝 值: {value_str}")
            
            # 如果是数字类型，显示统计信息
            if isinstance(value, (int, float)):
                print(f"     📈 数值: {value}")
            elif isinstance(value, str):
                print(f"     📏 字符串长度: {len(value)}")

def analyze_block_list(block_list):
    """
    详细分析block列表内容
    """
    print(f"\n📋 block列表长度: {len(block_list)}")
    
    for i, item in enumerate(block_list):
        if i >= 10:  # 只显示前10个元素
            print(f"\n  ... (还有{len(block_list)-10}个元素)")
            break
            
        item_type = type(item).__name__
        print(f"\n  📌 元素 {i}: [{item_type}]")
        
        if isinstance(item, dict):
            print(f"     📊 字典大小: {len(item)} 个键")
            if len(item) > 0:
                print("     📋 键列表:")
                for sub_key in list(item.keys())[:5]:
                    sub_value = item[sub_key]
                    print(f"       • {repr(sub_key)} ({type(sub_value).__name__})")
                if len(item) > 5:
                    print(f"       ... (还有{len(item)-5}个键)")
                    
        elif isinstance(item, list):
            print(f"     📊 列表长度: {len(item)}")
            if len(item) > 0:
                print("     📌 前几个元素:")
                for j, sub_item in enumerate(item[:3]):
                    print(f"       [{j}] {type(sub_item).__name__}: {repr(sub_item)[:80]}{'...' if len(repr(sub_item)) > 80 else ''}")
        else:
            item_str = repr(item)
            if len(item_str) > 150:
                print(f"     📝 值: {item_str[:150]}...")
            else:
                print(f"     📝 值: {item_str}")

def analyze_block_tuple(block_tuple):
    """
    详细分析block元组内容
    """
    print(f"\n🔘 block元组长度: {len(block_tuple)}")
    
    for i, item in enumerate(block_tuple):
        if i >= 10:  # 只显示前10个元素
            print(f"\n  ... (还有{len(block_tuple)-10}个元素)")
            break
            
        item_type = type(item).__name__
        print(f"\n  🎯 元素 {i}: [{item_type}]")
        
        if isinstance(item, dict):
            print(f"     📊 字典大小: {len(item)} 个键")
        elif isinstance(item, (list, tuple)):
            print(f"     📊 长度: {len(item)}")
        else:
            item_str = repr(item)
            if len(item_str) > 100:
                print(f"     📝 值: {item_str[:100]}...")
            else:
                print(f"     📝 值: {item_str}")

def analyze_all_top_level_keys(data):
    """
    分析所有顶层键的内容（当没有找到block键时）
    """
    print(f"\n📋 分析所有顶层键内容:")
    
    for i, (key, value) in enumerate(data.items()):
        if i >= 5:  # 只分析前5个键
            print(f"\n  ... (还有{len(data)-5}个键未显示)")
            break
            
        print(f"\n  🔑 键: {repr(key)}")
        print(f"     📦 类型: {type(value).__name__}")
        
        if isinstance(value, dict):
            print(f"     📊 大小: {len(value)} 个键")
            if len(value) > 0:
                print("     📋 子键:")
                for sub_key in list(value.keys())[:5]:
                    print(f"       • {repr(sub_key)}")
                if len(value) > 5:
                    print(f"       ... (还有{len(value)-5}个键)")
        elif isinstance(value, (list, tuple)):
            print(f"     📊 长度: {len(value)}")
        else:
            value_str = repr(value)
            if len(value_str) > 100:
                print(f"     📝 值: {value_str[:100]}...")
            else:
                print(f"     📝 值: {value_str}")

def find_first_pkl_file(directory):
    """
    在指定目录中查找第一个pkl文件
    
    Args:
        directory (str): 目录路径
    
    Returns:
        str: 第一个pkl文件的路径，如果没找到返回None
    """
    if not os.path.exists(directory):
        print(f"❌ 目录不存在: {directory}")
        return None
    
    # 按文件名排序，确保一致性
    files = sorted([f for f in os.listdir(directory) if f.endswith('.pkl')])
    
    if files:
        return os.path.join(directory, files[0])
    
    print(f"❌ 在目录 {directory} 中未找到pkl文件")
    return None

def main():
    # 请修改为你的实际目录路径
    target_directory = r"D:\desk\study5_COPDxLC_SMR\Aresult_9_1\ZNEW_unit_passSOP\ceshi"  # 修改为你的pkl文件所在目录
    
    print(f"🔍 正在目录 '{target_directory}' 中查找第一个pkl文件...")
    
    # 查找第一个pkl文件
    pkl_file = find_first_pkl_file(target_directory)
    
    if pkl_file:
        print(f"✅ 找到文件: {os.path.basename(pkl_file)}")
        print("=" * 80)
        print("🎯 详细分析 [block] 键内容:")
        print("=" * 80)
        
        # 分析block内容
        analyze_block_content(pkl_file)
        
        print("\n" + "=" * 80)
        print("✅ 分析完成!")
    else:
        print("❌ 未找到可分析的pkl文件")

if __name__ == "__main__":
    main()

'''
所以，我们知道block的顺序内容中包含的SNP的顺序，我们也知道它们之间的R
block名字是用其第一个snp的名称指代的
block的Z是要现场再计算的

我们需要对于每个这样的block，计算出来相关值，再回去原本的csv文件去求相互解释力
'''