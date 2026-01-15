# -*- coding: utf-8 -*-
# @Time    : 2023/6/14 20:58
# @Author  : lijianxin
# @File    : createNodeFun.py
# @Software: PyCharm
import sys
import sugardir.sugar_1 as sugar1
import sugardir.sugar_2 as sugar2
import sugardir.sugar_3 as sugar3
import sugardir.sugar_4 as sugar4
import sugardir.sugar_5 as sugar5
import sugardir.sugar_6 as sugar6
import sugardir.sugar_7 as sugar7
import sugardir.sugar_8 as sugar8

sugatlistdic = {"glycan_1": ['Gal', 'Glc', 'Gul', 'Man','HexA'],
                "glycan_2": ['Xyl', 'Ara', 'Fuc', 'Rha', 'GlcA', 'IdoA'],
                "glycan_3": ['GalNAc', 'GlcNAc','GlcNS'],
                "glycan_4": ['MurNAc'],
                "glycan_5": ['Neu5Ac', 'Neu5Gc'],
                "glycan_6": ['KDN'],
                "glycan_7": ['KDO'],
                "glycan_8": ['Neu5,9Ac2']}
# 其中所有单糖变成小写，做到不区分大小写都可以匹配
sugatlistdic = {key: [value.lower() for value in values] for key, values in sugatlistdic.items()}
def createNode(sugarlab,sugarlabvalue,sugarnodelist,num,depth):
    # 对得到的标签值更出全部转成小写
    # 根据单糖分子 创建合适的节点，然后存放节点列表，
    for key, value in sugatlistdic.items():
        for value_key in value:
            if value_key == sugarlabvalue:
                sugarkey = key
                if sugarkey == "glycan_1":
                    # 首先对糖节点列表进行判断，要是列表是空的是需要创建节点做为第一个根节点root,同时把第一个节点的父亲节点标记为root
                    # 判断是不是为空
                    if not sugarnodelist:
                        root = sugar1.builtGlycan1Node(sugarlab,sugarlabvalue)
                        # 创建根节点，根节点的parent 是root,
                        root.anomer.items()
                        root.parent = "root"
                        root.anomer["c"] = None
                        sugarnodelist.append(root)
                        num[f'num{depth}'] += 1

                    # 如果节点列表中不空以下操作

                    # 2、需要对新创建的节点指定父亲节点，即列表中的最后一个节点就是新创建节点的父亲节点。
                    # 3、需要对列表中的最后一个节点创建几号碳原子链接的孩子节点，即新创建的节点需要加入列表最后节点的孩子当中
                    # 4、要注意当一个父亲节点下面的几号碳原子后面链接的孩子节点，
                    else:
                        # 创建新的节点
                        Node1 = sugar1.builtGlycan1Node(sugarlab,sugarlabvalue)
                        # 添加到节点列表
                        sugarnodelist.append(Node1)
                        num[f'num{depth}'] += 1
                        # 节点列表中的元素增加

                        # 首先是对新创建的节点的指定父亲节点，当前节点sugarnodelist[sugarnodenum]，父亲节点sugarnodelist[sugarnodenum-1]
                        # 需要对父亲节点判断，如果父亲节点是正常的，才可以加入，或者需要报错或者异常处理。
                        if len(sugarnodelist[-1].parent) != 0 and sugarnodelist[-1].parent[0] == "1":#连接父节点
                            sugarnodelist[-1].parent.append(sugarnodelist[-2])
                        else:
                            print("此节点是非根节点，但此节点没有父节点，因此不能创建父亲")


                        # 目前这里出现问题，创建当前节点后，需要给上一个节点添加孩子，因此上个节点中有哪些孩子，目前无法直接获得，下面的写法，如果就不存在某些孩子，就会报错。
                        # ********************************************************************************************

                        if sugarnodelist[-2].childflagnum <= 1:#连接孩子节点
                            for chidrenkey, chidrenvalue in sugarnodelist[-2].childrenlistdic.items():
                                if len(chidrenvalue) != 0 and chidrenvalue[0] == chidrenkey:
                                    del chidrenvalue[0]
                                    chidrenvalue.append(sugarnodelist[-1])#给父节点添加子节点
                                    sugarnodelist[-2].childflagnum -= 1
                        else:
                            print("其中有孩子节点添加失败,退出程序")
                            sys.exit()
                    # 节点创建以及链接关系创建完成，继续下一个节点的创建，把label变成空
                elif sugarkey == "glycan_2":
                    # 首先对糖节点列表进行判断，要是列表是空的是需要创建节点做为第一个根节点root,同时把第一个节点的父亲节点标记为root
                    # 判断是不是为空
                    if not sugarnodelist:
                        root = sugar2.builtGlycan2Node(sugarlab,sugarlabvalue)
                        # 创建根节点，根节点的parent 是root,
                        root.parent = "root"
                        root.anomer["c"] = None
                        sugarnodelist.append(root)
                        num[f'num{depth}'] += 1

                    else:
                        # 创建新的节点
                        Node1 = sugar2.builtGlycan2Node(sugarlab,sugarlabvalue)
                        # 添加到节点列表
                        sugarnodelist.append(Node1)
                        # 当前层次添加的节点数
                        num[f'num{depth}'] += 1
                        # 节点列表中的元素增加
                        # 首先是对新创建的节点的指定父亲节点，当前节点sugarnodelist[sugarnodenum]，父亲节点sugarnodelist[sugarnodenum-1]
                        # 需要对父亲节点判断，如果父亲节点是正常的，才可以加入，或者需要报错或者异常处理。
                        if len(sugarnodelist[-1].parent) != 0 and sugarnodelist[-1].parent[0] == "1":
                            sugarnodelist[-1].parent.append(sugarnodelist[-2])
                        else:
                            print("此节点是非根节点，但此节点没有父节点，因此不能创建父亲")

                        if sugarnodelist[-2].childflagnum <= 1:
                            for chidrenkey, chidrenvalue in sugarnodelist[-2].childrenlistdic.items():
                                if len(chidrenvalue) != 0 and chidrenvalue[0] == chidrenkey:
                                    del chidrenvalue[0]
                                    chidrenvalue.append(sugarnodelist[-1])
                                    sugarnodelist[-2].childflagnum -= 1
                        else:
                            print("其中有孩子节点添加失败")
                            continue
                    # 节点创建以及链接关系创建完成，继续下一个节点的创建，把label变成空

                elif sugarkey == "glycan_3":
                    # 首先对糖节点列表进行判断，要是列表是空的是需要创建节点做为第一个根节点root,同时把第一个节点的父亲节点标记为root
                    # 判断是不是为空
                    if not sugarnodelist:#是空的，创建根节点
                        root = sugar3.builtGlycan3Node(sugarlab,sugarlabvalue)
                        # 创建根节点，根节点的parent 是root,
                        root.parent = "root"
                        root.anomer["c"] = None
                        sugarnodelist.append(root)
                        num[f'num{depth}'] += 1

                    # 如果节点列表中不空以下操作

                    # 2、需要对新创建的节点指定父亲节点，即列表中的最后一个节点就是新创建节点的父亲节点。
                    # 3、需要对列表中的最后一个节点创建几号碳原子链接的孩子节点，即新创建的节点需要加入列表最后节点的孩子当中
                    # 4、要注意当一个父亲节点下面的几号碳原子后面链接的孩子节点，
                    else:
                        # 创建新的节点
                        Node1 = sugar3.builtGlycan3Node(sugarlab,sugarlabvalue)
                        # 添加到节点列表
                        sugarnodelist.append(Node1)#添加的是一个类
                        # 当前层次添加的节点数
                        num[f'num{depth}'] += 1#第0层有几个元素
                        # 节点列表中的元素增加

                        # 首先是对新创建的节点的指定父亲节点，当前节点sugarnodelist[sugarnodenum]，父亲节点sugarnodelist[sugarnodenum-1]
                        # 需要对父亲节点判断，如果父亲节点是正常的，才可以加入，或者需要报错或者异常处理。
                        if len(sugarnodelist[-1].parent) != 0 and sugarnodelist[-1].parent[0] == "1":
                            sugarnodelist[-1].parent.append(sugarnodelist[-2])
                        else:
                            print("此节点是非根节点，但此节点没有父节点，因此不能创建父亲")


                        # 目前这里出现问题，创建当前节点后，需要给上一个节点添加孩子，因此上个节点中有哪些孩子，目前无法直接获得，下面的写法，如果就不存在某些孩子，就会报错。

                        if sugarnodelist[-2].childflagnum == 1:#等于0，没有被标记
                            for chidrenkey, chidrenvalue in sugarnodelist[-2].childrenlistdic.items():
                                if len(chidrenvalue) != 0 and chidrenvalue[0] == chidrenkey:
                                    del chidrenvalue[0]#给真正的子节点挪位置，删除之前标记的数字
                                    chidrenvalue.append(sugarnodelist[-1])#该节点位置添加一个对象孩子节点
                                    sugarnodelist[-2].childflagnum -= 1
                        else:
                            print("其中有孩子节点添加失败")
                            continue
                    # 节点创建以及链接关系创建完成，继续下一个节点的创建，把label变成空

                elif sugarkey == "glycan_4":
                    # 首先对糖节点列表进行判断，要是列表是空的是需要创建节点做为第一个根节点root,同时把第一个节点的父亲节点标记为root
                    # 判断是不是为空
                    if not sugarnodelist:
                        root = sugar4.builtGlycan4Node(sugarlab,sugarlabvalue)
                        # 创建根节点，根节点的parent 是root,
                        root.parent = "root"
                        root.anomer["c"] = None
                        sugarnodelist.append(root)
                        num[f'num{depth}'] += 1

                    # 如果节点列表中不空以下操作

                    # 2、需要对新创建的节点指定父亲节点，即列表中的最后一个节点就是新创建节点的父亲节点。
                    # 3、需要对列表中的最后一个节点创建几号碳原子链接的孩子节点，即新创建的节点需要加入列表最后节点的孩子当中
                    # 4、要注意当一个父亲节点下面的几号碳原子后面链接的孩子节点，
                    else:
                        # 创建新的节点
                        Node1 = sugar4.builtGlycan4Node(sugarlab,sugarlabvalue)
                        # 添加到节点列表
                        sugarnodelist.append(Node1)
                        # 当前层次添加的节点数
                        num[f'num{depth}'] += 1
                        # 节点列表中的元素增加
                        # 首先是对新创建的节点的指定父亲节点，当前节点sugarnodelist[sugarnodenum]，父亲节点sugarnodelist[sugarnodenum-1]
                        # 需要对父亲节点判断，如果父亲节点是正常的，才可以加入，或者需要报错或者异常处理。
                        if len(sugarnodelist[-1].parent) != 0 and sugarnodelist[-1].parent[0] == "1":
                            sugarnodelist[-1].parent.append(sugarnodelist[-2])
                        else:
                            print("此节点是非根节点，但此节点没有父节点，因此不能创建父亲")


                        if sugarnodelist[-2].childflagnum <= 1:
                            for chidrenkey, chidrenvalue in sugarnodelist[-2].childrenlistdic.items():
                                if len(chidrenvalue) != 0 and chidrenvalue[0] == chidrenkey:
                                    del chidrenvalue[0]
                                    chidrenvalue.append(sugarnodelist[-1])
                                    sugarnodelist[-2].childflagnum -= 1
                        else:
                            print("其中有孩子节点添加失败")
                            continue
                    # 节点创建以及链接关系创建完成，继续下一个节点的创建，把label变成空

                elif sugarkey == "glycan_5":
                    # 首先对糖节点列表进行判断，要是列表是空的是需要创建节点做为第一个根节点root,同时把第一个节点的父亲节点标记为root
                    # 判断是不是为空
                    if not sugarnodelist:
                        root = sugar5.builtGlycan5Node(sugarlab,sugarlabvalue)
                        # 创建根节点，根节点的parent 是root,
                        root.parent = "root"
                        root.anomer["c"] = None
                        sugarnodelist.append(root)
                        num[f'num{depth}'] += 1

                    # 如果节点列表中不空以下操作

                    # 2、需要对新创建的节点指定父亲节点，即列表中的最后一个节点就是新创建节点的父亲节点。
                    # 3、需要对列表中的最后一个节点创建几号碳原子链接的孩子节点，即新创建的节点需要加入列表最后节点的孩子当中
                    # 4、要注意当一个父亲节点下面的几号碳原子后面链接的孩子节点，
                    else:
                        # 创建新的节点
                        Node1 = sugar5.builtGlycan5Node(sugarlab,sugarlabvalue)
                        # 添加到节点列表
                        sugarnodelist.append(Node1)
                        # 当前层次添加的节点数
                        num[f'num{depth}'] += 1
                        # 节点列表中的元素增加
                        # 首先是对新创建的节点的指定父亲节点，当前节点sugarnodelist[sugarnodenum]，父亲节点sugarnodelist[sugarnodenum-1]
                        # 需要对父亲节点判断，如果父亲节点是正常的，才可以加入，或者需要报错或者异常处理。
                        if len(sugarnodelist[-1].parent) != 0 and sugarnodelist[-1].parent[0] == "2":
                            sugarnodelist[-1].parent.append(sugarnodelist[-2])
                        else:
                            print("此节点是非根节点，但此节点没有父节点，因此不能创建父亲")

                        if sugarnodelist[-2].childflagnum <= 1:
                            for chidrenkey, chidrenvalue in sugarnodelist[-2].childrenlistdic.items():
                                if len(chidrenvalue) != 0 and chidrenvalue[0] == chidrenkey:
                                    del chidrenvalue[0]
                                    chidrenvalue.append(sugarnodelist[-1])
                                    sugarnodelist[-2].childflagnum -= 1
                        else:
                            print("其中有孩子节点添加失败")
                            continue
                    # 节点创建以及链接关系创建完成，继续下一个节点的创建，把label变成空

                elif sugarkey == "glycan_6":
                    # 首先对糖节点列表进行判断，要是列表是空的是需要创建节点做为第一个根节点root,同时把第一个节点的父亲节点标记为root
                    # 判断是不是为空
                    if not sugarnodelist:
                        root = sugar6.builtGlycan6Node(sugarlab,sugarlabvalue)
                        # 创建根节点，根节点的parent 是root,
                        root.parent = "root"
                        root.anomer["c"] = None
                        sugarnodelist.append(root)
                        num[f'num{depth}'] += 1

                    # 如果节点列表中不空以下操作

                    # 2、需要对新创建的节点指定父亲节点，即列表中的最后一个节点就是新创建节点的父亲节点。
                    # 3、需要对列表中的最后一个节点创建几号碳原子链接的孩子节点，即新创建的节点需要加入列表最后节点的孩子当中
                    # 4、要注意当一个父亲节点下面的几号碳原子后面链接的孩子节点，
                    else:
                        # 创建新的节点
                        Node1 = sugar6.builtGlycan6Node(sugarlab,sugarlabvalue)
                        # 添加到节点列表
                        sugarnodelist.append(Node1)
                        # 当前层次添加的节点数
                        num[f'num{depth}'] += 1
                        # 节点列表中的元素增加
                        # 首先是对新创建的节点的指定父亲节点，当前节点sugarnodelist[sugarnodenum]，父亲节点sugarnodelist[sugarnodenum-1]
                        # 需要对父亲节点判断，如果父亲节点是正常的，才可以加入，或者需要报错或者异常处理。
                        if len(sugarnodelist[-1].parent) != 0 and sugarnodelist[-1].parent[0] == "2":
                            sugarnodelist[-1].parent.append(sugarnodelist[-2])
                        else:
                            print("此节点是非根节点，但此节点没有父节点，因此不能创建父亲")

                        if sugarnodelist[-2].childflagnum <= 1:
                            for chidrenkey, chidrenvalue in sugarnodelist[-2].childrenlistdic.items():
                                if len(chidrenvalue) != 0 and chidrenvalue[0] == chidrenkey:
                                    del chidrenvalue[0]
                                    chidrenvalue.append(sugarnodelist[-1])
                                    sugarnodelist[-2].childflagnum -= 1
                        else:
                            print("其中有孩子节点添加失败")
                            continue
                    # 节点创建以及链接关系创建完成，继续下一个节点的创建，把label变成空

                elif sugarkey == "glycan_7":
                    # 首先对糖节点列表进行判断，要是列表是空的是需要创建节点做为第一个根节点root,同时把第一个节点的父亲节点标记为root
                    # 判断是不是为空
                    if not sugarnodelist:
                        root = sugar7.builtGlycan7Node(sugarlab,sugarlabvalue)
                        # 创建根节点，根节点的parent 是root,
                        root.parent = "root"
                        root.anomer["c"] = None
                        sugarnodelist.append(root)
                        num[f'num{depth}'] += 1

                    # 如果节点列表中不空以下操作

                    # 2、需要对新创建的节点指定父亲节点，即列表中的最后一个节点就是新创建节点的父亲节点。
                    # 3、需要对列表中的最后一个节点创建几号碳原子链接的孩子节点，即新创建的节点需要加入列表最后节点的孩子当中
                    # 4、要注意当一个父亲节点下面的几号碳原子后面链接的孩子节点，
                    else:
                        # 创建新的节点
                        Node1 = sugar7.builtGlycan7Node(sugarlab,sugarlabvalue)
                        # 添加到节点列表
                        sugarnodelist.append(Node1)
                        # 当前层次添加的节点数
                        num[f'num{depth}'] += 1
                        # 节点列表中的元素增加
                        # 首先是对新创建的节点的指定父亲节点，当前节点sugarnodelist[sugarnodenum]，父亲节点sugarnodelist[sugarnodenum-1]
                        # 需要对父亲节点判断，如果父亲节点是正常的，才可以加入，或者需要报错或者异常处理。
                        if len(sugarnodelist[-1].parent) != 0 and sugarnodelist[-1].parent[0] == "2":
                            sugarnodelist[-1].parent.append(sugarnodelist[-2])
                        else:
                            print("此节点是非根节点，但此节点没有父节点，因此不能创建父亲")

                        if sugarnodelist[-2].childflagnum <= 1:
                            for chidrenkey, chidrenvalue in sugarnodelist[-2].childrenlistdic.items():
                                if len(chidrenvalue) != 0 and chidrenvalue[0] == chidrenkey:
                                    del chidrenvalue[0]
                                    chidrenvalue.append(sugarnodelist[-1])
                                    sugarnodelist[-2].childflagnum -= 1
                        else:
                            print("其中有孩子节点添加失败")
                            continue
                    # 节点创建以及链接关系创建完成，继续下一个节点的创建，把label变成空

                elif sugarkey == "glycan_8":
                    # 首先对糖节点列表进行判断，要是列表是空的是需要创建节点做为第一个根节点root,同时把第一个节点的父亲节点标记为root
                    # 判断是不是为空
                    if not sugarnodelist:
                        root = sugar8.builtGlycan8Node(sugarlab,sugarlabvalue)
                        # 创建根节点，根节点的parent 是root,
                        root.parent = "root"
                        root.anomer["c"] = None
                        sugarnodelist.append(root)
                        num[f'num{depth}'] += 1

                    # 如果节点列表中不空以下操作

                    # 2、需要对新创建的节点指定父亲节点，即列表中的最后一个节点就是新创建节点的父亲节点。
                    # 3、需要对列表中的最后一个节点创建几号碳原子链接的孩子节点，即新创建的节点需要加入列表最后节点的孩子当中
                    # 4、要注意当一个父亲节点下面的几号碳原子后面链接的孩子节点，
                    else:
                        # 创建新的节点
                        Node1 = sugar8.builtGlycan8Node(sugarlab,sugarlabvalue)
                        # 添加到节点列表
                        sugarnodelist.append(Node1)
                        # 当前层次添加的节点数
                        num[f'num{depth}'] += 1
                        # 节点列表中的元素增加
                        # 首先是对新创建的节点的指定父亲节点，当前节点sugarnodelist[sugarnodenum]，父亲节点sugarnodelist[sugarnodenum-1]
                        # 需要对父亲节点判断，如果父亲节点是正常的，才可以加入，或者需要报错或者异常处理。
                        if len(sugarnodelist[-1].parent) != 0 and sugarnodelist[-1].parent[0] == "2":
                            sugarnodelist[-1].parent.append(sugarnodelist[-2])
                        else:
                            print("此节点是非根节点，但此节点没有父节点，因此不能创建父亲")

                        if sugarnodelist[-2].childflagnum <= 1:
                            for chidrenkey, chidrenvalue in sugarnodelist[-2].childrenlistdic.items():
                                if len(chidrenvalue) != 0 and chidrenvalue[0] == chidrenkey:
                                    del chidrenvalue[0]
                                    chidrenvalue.append(sugarnodelist[-1])
                                    sugarnodelist[-2].childflagnum -= 1
                        else:
                            print("其中有孩子节点添加失败")
                            continue
                    # 节点创建以及链接关系创建完成，继续下一个节点的创建，把label变成空

if __name__ == '__main__':
    createNode()