import sys
from sugardir import createNodeFun as nodebuilt
from sugardir import iterate_tree
from sugardir import outputfuc as output


# 开始处理带有括弧形状的节点,一下情况全部可以处理。
#testsugar = "Gala1-3GalNAcb-Sp8"
#testsugar = "Galb1-3-(Neu5Aca2-3Galb1-4GlcNacb1-6)GalNAc-T"
# testsugar = "GlcNAca1-3(4S6SGlcNacb1-4)GlcNAc"
# testsugar = "GlcNAcb1-3(4S6SGlcNAcb1-4)2SGlcA"
# testsugar = "Fuca1-2Galb1-3GlcNAcb1-3(Fuca1-3(Galb1-4)GlcNAcb1-6)Galb1-4Glc"
# testsugar = "Mana1-3(GalNAca1-3(Fuca1-6)4SGlcNSb1-46SGlcNAcb1-2Mana1-6)Manb1-43SGlcNAc"
# testsugar = "Fuca1-2Galb1-4(Fuca1-3)GlcNAcb1-2Mana1-2(Fuca1-26SGalb1-4(Fuca1-3)GlcNAcb1-2Mana1-4)3SGlcAb1-4GlcNAcb1-4(Fuca1-6)3SGlcNAc"
# # testsugar = "GlcNAcb1-64SGlcNSa1-6(4S6SGlcNAcb1-4)3SGlcNS"
# testsugar = "Fuca1-2Galb1-4(Fuca1-3)GlcNAcb1-3(Neu5Aca2-6Galb1-4GlcNAcb1-4)IdoAb1-3Glc"
# testsugar = "Fuca1-2Galb1-46SGlcNAc-Sp8"



label = []
num = {}
depth = 0
num[f'num{depth}'] = 0
# 创建一个存放糖节点的一个列表
glycan = ['Gal', 'Glc', 'Glu', 'Man', 'HexA', 'Xyl', 'Ara', 'Fuc', 'Rha', 'GlcA', 'IdoA', 'GalNAc', 'GlcNAc',
          'GlcNS', 'Neu5Ac', 'Neu5Gc', 'KDN', 'KDO', 'Neu5,9Ac2']
extra = ["Sp",'Gly','T','LVaNKT','MDPlys','6AA','N','Asn','(OCH2CH2)6NH2','NLTAVL','MDPLys','LVANKT']
with open(r"D:\glycan1.csv", "r") as sugarfile:
    for sugar in sugarfile:
        testsugar = sugar.replace('"', "").strip('\n')
        print(testsugar)
        # 首先对括弧的数量进行统计，左右括弧的数量相等，或者直接程序退出
        sugarnodelist = []
        count = 0
        for char in testsugar:
            if char == ")":
                count += 1
            elif char == "(":
                count -= 1
        if count != 0:
            print("糖链结构中的括弧数量不相等")
        glycan_list = testsugar.split('-')
        if 'Sp' in glycan_list[-1]:
            glycan_list[-1] = 'Sp'
        for i in extra:
            if i == glycan_list[-1]:
                testsugar = '-'.join(glycan_list[0:-1])

        index = testsugar[-3:len(testsugar)]
        if index == 'Rha':
            testsugar = "-0" + testsugar + "c0"  # 人为添加
        elif index == 'Ara':
            testsugar = "-0" + testsugar + "c0"  # 人为添加
        elif testsugar[-1] == "a" or testsugar[-1] == "b":
            testsugar = "-0" + testsugar + "0"
        else:
            testsugar = "-0" + testsugar + "c0"
        #print(testsugar)# 没有a和b，自行添加

        for s in testsugar[::-1]:
            label.append(s)
            if s == "-":
                sugarlab = label[:-1][::-1]
                sugarlab = "".join(sugarlab)
                # (Neu5Aca2-3Galb1-3)((Neu5Aca2-3Galb1-4(Fuca1-3))GlcNAcb1-6)GalNAca1-4Glc
                if "(" not in sugarlab:
                    sugarlabreplace = sugarlab.replace(")", "").replace("(", "")
                    monosaccharide = ''
                    for i in sugarlabreplace[0:-2][::-1]:  # 对43SGalb0中的43SGal从右向左遍历
                        monosaccharide = i + monosaccharide  # 将字符串不断添加到新的空字符串中
                        if monosaccharide in glycan:  # 判断单糖是否出现在聚糖字典中
                            sugarlabvalue = monosaccharide.lower()
                            monosaccharide = ''
                    nodebuilt.createNode(sugarlabreplace, sugarlabvalue, sugarnodelist, num, depth)
                    # 处理括号
                    for char in sugarlab[::-1]:  # 看是否有括号
                        if char == ")":
                            depth += 1  #
                            num[f'num{depth}'] = 0  # 这层所含的糖的个数
                    label = []  # 标签变为空
                else:
                    tempchild = str(sugarlab[0])  # 把连接位点变为0，说明此节点为叶子节点
                    sugarlab1 = "0" + sugarlab[1:]
                    sugarlabreplace = sugarlab1.replace(")", "").replace("(", "")
                    monosaccharide = ''
                    for i in sugarlabreplace[0:-2][::-1]:  # 对43SGalb0中的43SGal从右向左遍历
                        monosaccharide = i + monosaccharide  # 将字符串不断添加到新的空字符串中
                        if monosaccharide in glycan:  # 判断单糖是否出现在聚糖字典中
                            sugarlabvalue = monosaccharide.lower()
                            monosaccharide = ''
                    nodebuilt.createNode(sugarlabreplace, sugarlabvalue, sugarnodelist, num, depth)
                    for char in sugarlab[::-1]:
                        if char == "(":
                            while num[f'num{depth}'] > 0:  # 表示该层有糖
                                sugarnodelist.pop()
                                num[f'num{depth}'] -= 1  # 糖的数量减少1
                            depth -= 1  # 等糖全部删除，层次减1
                        elif char == ")":  # 测试')('这种情况
                            depth += 1
                            num[f'num{depth}'] = 0
                        else:
                            continue
                    if sugarnodelist[-1].childflagnum == 0:  # 处理括号外节点的糖之间的连接
                        for chidrenkey, chidrenvalue in sugarnodelist[-1].childrenlistdic.items():  # 遍历节点的孩子，childrenvalue指添加子节点
                            if tempchild == chidrenkey and len(chidrenvalue) != 0:  # 3，出现两个占位置，报错
                                print("请检查糖链结构", "".join(label)[::-1])
                                sys.exit()
                            elif tempchild == chidrenkey and len(chidrenvalue) == 0:  # 没有出现两个占位，3
                                chidrenvalue.append(tempchild)  # 等效 chidrenvalue = [tempchild]
                                sugarnodelist[-1].childflagnum += 1  # 为什么这边又开始标记，表示该层结束。
                    label = []
        # print("转化成功......")
        # 对数进行遍历，从根节点开始
        # print("开始对转成的树进行遍历......")
        singlelist = []
        doublelist = []
        triplelist = []
        NodeGrp = []

        root = sugarnodelist[0]  # 从根节点开始遍历
        iterate_tree.preorder_print(root, singlelist, doublelist, triplelist, testsugar)  # 扫描根节点


        def getNodeGrp():
            for sugar in doublelist:
                sugartem = sugar.split("-")
                anomernode = sugartem[0][-1] + "-" + sugartem[1][0]
                node1 = "(" + sugartem[0][:-1] + "," + anomernode + ")"
                node2 = "(" + anomernode + "," + sugartem[1][1:] + ")"
                NodeGrp.append(node1)
                NodeGrp.append(node2)
            #print(NodeGrp)
            print("\n")


        getNodeGrp()
        output.writesugar(sugar, singlelist, doublelist, triplelist, NodeGrp)
        #output.monosaccharide(sugar, singlelist)
        # output.disaccharide(sugar, doublelist)
        # output.trisaccharide(sugar, triplelist)
        # output.monosaccharide_glycosidilinkage(sugar, NodeGrp)
