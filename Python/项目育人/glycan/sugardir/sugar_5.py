# -*- coding: utf-8 -*-
# @Time    : 2023/6/14 20:54
# @Author  : lijianxin
# @File    : sugar_5.py
# @Software: PyCharm
import sys
class glycan5:
    def __init__(self, label):
        self.label = label
        self.anomer = {'a': None, 'b': None,'c':None}
        self.parent = []

        self.children_4 = []
        self.children_7 = []
        self.children_8 = []
        self.children_9 = []

        self.childflagnum = 0
        self.childrenlistdic = {"4": self.children_4,
                                "7": self.children_7,
                                "8": self.children_8,
                                "9": self.children_9
                                }
def builtGlycan5Node(label1,sugarlabvalue):
    Glycan5Node = glycan5(label1[1:-2])
    if label1[-2] == "a":
        Glycan5Node.anomer["a"] = "a"
    elif label1[-2] == "b":
        Glycan5Node.anomer["b"] = "b"
    elif label1[-2] == "c":
        Glycan5Node.anomer["c"] = "c"
    else:
        print("annomer，不合法，请检查语法")
        sys.exit()

    # 需要对根节点处理
    try:
        if label1[-1] == "2":
            Glycan5Node.parent.append("2")
        elif label1[-1] == "0":
            Glycan5Node.parent = "root"
    except IndexError:
        print("父节点出现错误")
        sys.exit()
    # 这地方可以统计出标记的数量

    Modification_groups = label1[1:len(label1[1:-2]) - len(sugarlabvalue) + 1]
    search_index = []
    if 'S' or 's' in Modification_groups:
        for index, letter in enumerate(Modification_groups):
            if letter == 'S':
                search_index.append(Modification_groups[index - 1])  # [1]
            for SO3_index in search_index:
                for chidrenkey, chidrenvalue in Glycan5Node.childrenlistdic.items():
                    if SO3_index == chidrenkey:
                        chidrenvalue.append("SO3")
                search_index = []
        # 判断是否含P
    if 'P' or 'p' in Modification_groups:
        for index, letter in enumerate(Modification_groups):
            if letter == 'P':
                search_index.append(Modification_groups[index - 1])  # [1]
                for p_indx in search_index:
                    for chidrenkey, chidrenvalue in Glycan5Node.childrenlistdic.items():
                        if p_indx == chidrenkey:
                            chidrenvalue.append("P")
                    search_index = []
    if label1[0]!='0' and label1[0] not in Glycan5Node.childrenlistdic.keys():
        print(label1[0] + ':有问题')
        sys.exit()
    for childkey, childvalue in Glycan5Node.childrenlistdic.items():
        if label1[0] == childkey and len(childvalue) == 0:
            childvalue.append(childkey)
            Glycan5Node.childflagnum += 1
    # if Glycan5Node.childflagnum == 0:
    #     print("此节点是叶子节点，没有链接孩子")

    return Glycan5Node
if __name__ == '__main__':
    builtGlycan5Node()