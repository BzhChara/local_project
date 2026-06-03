import codecs
import re
import torch
import torch.nn as nn
import pandas as pd
import torch.nn.functional as F

def remove_invisible_characters(text):
    return re.sub(r'[\u200B-\u200D\uFEFF\u2060]+', '', text)
def remove_numbers_and_sp(text):
    return re.sub(r'\d+[SP]+', '', text)
def process_line(line):
    line = remove_invisible_characters(line)
    parts = line.strip().split(';')
    number = parts[0]
    structure = parts[1]
    components = parts[2].strip()
    components_list = []

    components_list = [c.strip() for c in components.split(',')]
    components_list = list(filter(lambda item: 'Neu5' not in item, components_list))
    components_list = [
        item.replace('9Ac', 'Neu5,9Ac') if '9Ac' in item else item
        for item in components_list
    ]
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
                    if ':' in component:
                        glycan = component.split(':')[0].strip()
                        glycan=remove_numbers_and_sp(glycan)
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
        single_list = sorted(list(single_set))
        double_list = sorted(list(double_set))
        trip_list = sorted(list(trip_set))
    encode_data = []
    count_components(file_path, single_list, double_list, trip_list, encode_data)
    output_file_path = 'glycan_encode.csv'
    with open(output_file_path, 'w+') as output_file:
        encodelist = []
        for structure, number, encode_list in encode_data:
            if number == 'single':
                encodelist.extend(encode_list)
            elif number == 'double':
                encodelist.extend(encode_list)
            elif number == 'trip':
                encodelist.extend(encode_list)
                encode_str = ','.join(map(str, encodelist))
                output_file.write(f'{encode_str}\n')
                encodelist = []
def count_components(file_path,single_list,double_list,trip_list,encode_data):
    with (codecs.open(file_path, 'r', encoding='utf-8') as file):

        for line in file:
            encode_list_single = [0] * len(single_list)
            encode_list_double = [0] * len(double_list)
            encode_list_trip = [0] * len(trip_list)
            number, structure, components = process_line(line)
            for component in components:
                if ':' not in component:
                    continue
                glycan,count = component.split(':')
                count = int(count)
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
file_path = 'sugarlist.txt'
combine_array(file_path)