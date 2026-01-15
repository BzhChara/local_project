import codecs
import re
def remove_invisible_characters(text):
    # 移除所有 Unicode 定义的不可见字符
    return re.sub(r'[\u200B-\u200D\uFEFF\u2060]+', '', text)


def remove_numbers_and_sp(text):
    # 使用正则表达式删除所有数字及跟随的 S 或 P 字符，可能出现多次
    return re.sub(r'\d+[SP]+', '', text)
def process_line(line):
    # 移除不可见字符
    line = remove_invisible_characters(line)
    # Split the line by semicolons
    parts = line.strip().split(';')

    # The first part contains the glycan structure
    number = parts[0]
    # The secend part contains the glycan structure
    structure = parts[1]
    # The tride part contains the components and their counts
    components = parts[2].strip()
    # Split the components by commas
    components_list = []

    components_list = [c.strip() for c in components.split(',')]
    components_list = list(filter(lambda item: 'Neu5' not in item, components_list))
    components_list = [
        item.replace('9Ac', 'Neu5,9Ac') if '9Ac' in item else item
        for item in components_list
    ]
    # output_file_path = './output/components.txt'
    # # 一次性写入所有数据
    # with open(output_file_path, 'a') as output_file:
    #     for item in components_list:
    #         output_file.write(f"{item}\n")

    return number,structure, components_list

def combine_array(file_path):
    single_set=set()
    double_set = set()
    trip_set=set()
    with codecs.open(file_path, 'r', encoding='utf-8') as file:
        for line in file:
            number,structure, components = process_line(line)
            if number=='single':
                for component in components:
                    # 提取冒号前的内容
                    if ':' in component:
                        glycan = component.split(':')[0].strip()
                        glycan=remove_numbers_and_sp(glycan)
                        # 添加到集合中
                        single_set.add(glycan)
            if number=='double':
                for component in components:
                    if ':' in component:
                        glycan = component.split(':')[0].strip()
                        double_set.add(glycan)
            if number == 'trip':
                for component in components:
                    if ':' in component:
                        glycan = component.split(':')[0].strip()
                        trip_set.add(glycan)
        # 转换为列表并排序
        single_list = sorted(list(single_set))
        double_list = sorted(list(double_set))
        trip_list = sorted(list(trip_set))
    encode_data=[]
    count_components(file_path,single_list,double_list,trip_list,encode_data)
    output_file_path = './output/glycan_encode2.txt'
    # 一次性写入所有数据
    with open(output_file_path, 'w') as output_file:
        for structure, number, encode_list in encode_data:
            if number == 'single':
                output_file.write(f"{structure}:{number}:{encode_list}\n")
            elif number == 'double':
                output_file.write(f"{structure}:{number}:{encode_list}\n")
            elif number == 'trip':
                output_file.write(f"{structure}:{number}:{encode_list}\n")
    # output_file_path = './output/combined_glycansout3.txt'
    # with open(output_file_path, 'w') as output_file:
    #     output_file.write("Single Set:\n")
    #     for glycan in single_list:
    #         output_file.write(f"{glycan}\n")
    #     output_file.write("\nDouble Set:\n")
    #     for glycan in double_list:
    #         output_file.write(f"{glycan}\n")
    #     output_file.write("\nTrip Set:\n")
    #     for glycan in trip_list:
    #         output_file.write(f"{glycan}\n")
def count_components(file_path,single_list,double_list,trip_list,encode_data):
    with (codecs.open(file_path, 'r', encoding='utf-8') as file):

        for line in file:
            encode_list_single = [0] * len(single_list)
            encode_list_double = [0] * len(double_list)
            encode_list_trip = [0] * len(trip_list)
            number, structure, components = process_line(line)
            for component in components:
                # 提取糖名和数量
                if ':' not in component:  # 如果没有冒号，则跳过这次迭代
                    continue
                glycan,count = component.split(':')
                count = int(count)  # 将数量转换为整数
                if number == 'single':
                    glycan=remove_numbers_and_sp(glycan)
                    index=single_list.index(glycan)
                    encode_list_single[index]+=count
                if number == 'double':
                    index = double_list.index(glycan)
                    encode_list_double[index] += count
                if number == 'trip':
                    index=trip_list.index(glycan)
                    encode_list_trip[index]+=count

            if number == 'single':
                encode_data.append((structure, 'single', encode_list_single))
            elif number == 'double':
                encode_data.append((structure, 'double', encode_list_double))
            elif number == 'trip':
                encode_data.append((structure, 'trip', encode_list_trip))

# Assuming the file is named 'glycans.txt' and is located in the current working directory.
file_path = './output/sugarlist2.txt'
combine_array(file_path)