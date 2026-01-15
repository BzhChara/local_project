# -*- coding: utf-8 -*-
# @Time    : 2023/6/15 14:57
# @Author  : lijianxin
# @File    : outputfuc.py
# @Software: PyCharm

outputfile = open("./output/sugarlist.txt","w")

def countsugar(sugarlist):
    sugardiccount = {}
    for sugar in sugarlist:
        if sugar in sugardiccount:
            sugardiccount[sugar] += 1
        else:
            sugardiccount[sugar] = 1
    return sugardiccount

def countsugar(sugarlist):
    sugardiccount = {}
    for sugar in sugarlist:
        if sugar in sugardiccount:
            sugardiccount[sugar] += 1
        else:
            sugardiccount[sugar] = 1
    return sugardiccount

def writesugar(sugar, singlelist, doublelist, triplelist, nodegrplist):
    #outputfile.write(sugar)
    # outputfile.write('single:'+ ';'.join(singlelist)+'\n')
    sigledic = countsugar(singlelist)
    key_value_single = [f"{key}:{value}" for key, value in sigledic.items()]
    outputfile.write('single;'+ sugar.strip("\n") + ';' + ','.join(key_value_single))
    outputfile.write('\n')
    # outputfile.write('double:'+ ';'.join(doublelist)+'\n')
    doubledic = countsugar(doublelist)
    key_value_double = [f"{key}:{value}" for key, value in doubledic.items()]
    outputfile.write('double;'+ sugar.strip("\n") + ';' + ','.join(key_value_double))
    outputfile.write("\n")
    # outputfile.write('trip:'+';'.join(triplelist)+'\n')
    tripledic = countsugar(triplelist)
    key_value_triple = [f"{key}:{value}" for key, value in tripledic.items()]
    outputfile.write('trip;'+ sugar.strip("\n") + ';' +','.join(key_value_triple))
    outputfile.write("\n")



    # for node_tem in nodegrplist:
    #     outputfile.write(node_tem)
    #     outputfile.write("\n")
if __name__ == '__main__':
    countsugar()
    writesugar()