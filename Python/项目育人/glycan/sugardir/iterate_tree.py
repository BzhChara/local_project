# -*- coding: utf-8 -*-
# @Time    : 2023/6/14 22:45
# @Author  : lijianxin
# @File    : iterate_tree.py
# @Software: PyCharm
from sugardir import outputfuc as otp
def preorder_print(sugarnode, singlelist, doublelist, triplelist,testsugar):#打印糖
    if sugarnode is None:
        return
    else:
        def getsinglesugar(sugarnodelable):#打印单糖
            single_anomer_keys = [key for key,value in sugarnodelable.anomer.items() if value]
            if len(single_anomer_keys) == 0:
                singlesugar = sugarnodelable.label
            else:
                singlesugar = f"{sugarnodelable.label}" + f"{single_anomer_keys[0]}"
            singlelist.append(singlesugar)
            # print("这是所有的单糖：", singlesugar)

        # 打印节点，如果是二糖就是节点以及孩子节点，如果是三糖就是节点，节点孩子，以及节点孩子的孩子。
        # 当前节点下的所有孩子，此时等于一个二糖，当前节点，以及当前节点的孩子节点。
        def getdoublesugar(sugarnodelable):
            childnum = 0
            for currentchildkey, currentchildvalue in sugarnodelable.childrenlistdic.items():
                if len(currentchildvalue) != 0 and currentchildvalue[0] != 'SO3' and currentchildvalue[0] != 'P':#表示有单糖连接
                    childnum += 1
                    # print("孩子节点：", currentchildvalue[0].label)
                    children_anomer_keys = [key for key, value in currentchildvalue[0].anomer.items() if value][0]# IF value不是空的，判断是a还是b
                    parent_anomer_keys = [key for key, value in sugarnodelable.anomer.items() if value]
                    if len(parent_anomer_keys) == 0:
                        doublesugar = f"{currentchildvalue[0].label}" + f"{children_anomer_keys}" + f"{currentchildvalue[0].parent[0]}" + "-" + f"{currentchildkey}" + f"{sugarnode.label}"
                    else:
                        doublesugar = f"{currentchildvalue[0].label}" + f"{children_anomer_keys}" + f"{currentchildvalue[0].parent[0]}" + "-" + f"{currentchildkey}" + f"{sugarnode.label}" + f"{parent_anomer_keys[0]}"#看单糖之间的连接:当前节点的标签+糖苷键类型+父亲节点连接+”-“+在父亲节点的连接方式+父亲节点
                    doublelist.append(doublesugar)
            #         print("组成的二糖:", doublesugar)
            # if childnum == 0:
            #     print("当前孩子没有节点")
            # else:
            #     print("当前节点的孩子节点输出结束")

        def gettriplesugar(sugarnodelable):#遍历三糖
            childnum = 0
            # 这是根据父亲节点得到的孩子节点
            for currentchildkey, currentchildvalue in sugarnodelable.childrenlistdic.items():#得到它连接的孩子节点
                if len(currentchildvalue) != 0 and currentchildvalue[0] != 'SO3' and currentchildvalue[0] != 'P':#当前孩子不是0
                    # 设置参数对孩子节点的孩子进行数量统计
                    childnum += 1#统计有几个孩子
                    # 这是得到的孩子节点,即父亲节点的孩子节点
                    # print("孩子节点：", currentchildvalue[0].label)#输出该糖得到的孩子
                    anomer_keys = [key for key, value in currentchildvalue[0].anomer.items() if value][0]#判断糖苷键的类型
                    childsecnum = 0
                    # 这是对孩子的节点进行遍历，即孩子节点的孩子，为了遍历三糖做准备
                    for currentchildseckey, currentchildsecvalue in currentchildvalue[0].childrenlistdic.items():#对当前孩子节点所连接的单糖进行遍历
                        if len(currentchildsecvalue) != 0 and currentchildsecvalue[0] != 'SO3' and currentchildsecvalue[0] != 'P':#判断连接的孩子的情况
                            childsecnum += 1#统计有几个孙子
                            # 这是得到的孩子节点的孩子，可以说是孙子节点，或者是节点孩子的孩子
                            # print("孩子节点的孩子(孙子节点)：", currentchildsecvalue[0].label)#输出它的孙子节点
                            anomersec_keys = [key for key, value in currentchildsecvalue[0].anomer.items() if value][0]#判断糖苷键类型
                            parent_anomer_keys = [key for key, value in sugarnodelable.anomer.items() if value]
                            if len(parent_anomer_keys) == 0:
                                triplesugar = f"{currentchildsecvalue[0].label}" + f"{anomersec_keys}" + f"{currentchildsecvalue[0].parent[0]}" + "-" + f"{currentchildseckey}" + f"{currentchildvalue[0].label}" + f"{anomer_keys}" + f"{currentchildvalue[0].parent[0]}" +"-"+ f"{currentchildkey}" + f"{sugarnode.label}"
                            else:
                                triplesugar = f"{currentchildsecvalue[0].label}" + f"{anomersec_keys}" + f"{currentchildsecvalue[0].parent[0]}" + "-" + f"{currentchildseckey}" + f"{currentchildvalue[0].label}" + f"{anomer_keys}" + f"{currentchildvalue[0].parent[0]}" + "-" + f"{currentchildkey}" + f"{sugarnode.label}" + f"{parent_anomer_keys[0]}"
                            triplelist.append(triplesugar)
            #                 print("组成的三糖：", triplesugar)
            #         if childsecnum == 0:
            #             print("当前孩子没有节点")
            #         else:
            #             print("当前节点的孩子节点输出结束")
            # if childnum == 0:
            #     print("当前孩子没有节点")
            # else:
            #     print("当前节点的孩子节点输出结束")
        # 选择打印需要糖的种类
        getsinglesugar(sugarnode)
        getdoublesugar(sugarnode)
        gettriplesugar(sugarnode)

    for childnode in sugarnode.childrenlistdic.values():#找该节点连接的孩子节点
        if len(childnode) != 0 and childnode[0] != 'SO3' and childnode[0] != 'P':
            preorder_print(childnode[0], singlelist, doublelist, triplelist,testsugar)
            
            
        #后添的代码
        root_index = testsugar.split('-')
        root_index.remove('')
        # print(len(root_index))
        # 找出根节点
        root_name = singlelist[0]
        same_rootname = []
        # 首先判断root是否是'Rha'和'Ara'，若是则则无糖苷键，统计出现次数，及和它相同糖苷键的次数
        if root_name == 'Ara' or root_name == 'Rha':
            singlelist.remove(root_name)
            rootcount = {root_name: 0}
            a_bcount = {'a': 0, 'b': 0}
            # 判断根节点与其他结点是否相同，统计相同的次数
            for i in singlelist:
                if root_name == i[0:-1]:
                    same_rootname.append(i)
                    rootcount[root_name] += 1
            # 统计具有相同名称的糖苷键出现次数
            for s in same_rootname:
                if s[-1] == 'a':
                    a_bcount['a'] += 1
                elif s[-1] == 'b':
                    a_bcount['b'] += 1
            # 限制
            rootcount_values = rootcount[root_name]

            if rootcount_values / int(len(root_index)) >= 0.5:
                if a_bcount['a'] >= a_bcount['b']:
                    root_name = root_name + 'a'
                elif a_bcount['a'] < a_bcount['b']:
                    root_name = root_name + 'b'
            singlelist.append(root_name)
        else:
            if root_name[-1] == 'a' or root_name[-1] == 'b':
                singlelist = singlelist
            else:  # 添加糖苷键
                singlelist.remove(root_name)
                rootcount = {root_name: 0}
                a_bcount = {'a': 0, 'b': 0}
                # 判断根节点和其他结点的关系
                for i in singlelist:
                    if root_name == i[0:-1]:
                        same_rootname.append(i)
                        rootcount[root_name] += 1
                for s in same_rootname:
                    if s[-1] == 'a':
                        a_bcount['a'] += 1
                    elif s[-1] == 'b':
                        a_bcount['b'] += 1
                # 限制
                rootcount_values = rootcount[root_name]
                if rootcount_values / int(len(root_index)) >= 0.5:
                    if a_bcount['a'] >= a_bcount['b']:
                        root_name = root_name + 'a'
                    elif a_bcount['a'] < a_bcount['b']:
                        root_name = root_name + 'b'
                singlelist.append(root_name)
        # print(singlelist)
        otp.countsugar(singlelist)
        # print(doublelist)
        otp.countsugar(doublelist)
        # print(triplelist)
        otp.countsugar(triplelist)
if __name__ == '__main__':
    preorder_print()