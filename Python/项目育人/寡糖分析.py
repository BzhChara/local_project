import re

data=[]

with open("sugarlist.txt", "r") as f:
    for line in f.readlines():
        line = line.strip('\n')  #去掉列表中每一个元素的换行符
        data.append(line)

print(data[0:3],"\n")

num1=dict()
num2=dict()
num3=dict()

print("单糖如下:")
for item1 in data:
    if item1[0]=='s':
        tmp1=re.findall(r'\w+[\,]{0,1}\w+:[0-9]',item1)
        for item2 in tmp1:
            if item2[1]!=',':
                num1[item2[0:-2]]=1
print(num1)

print("二糖如下:")
for item3 in data:
    if item3[0]=='d':
        tmp2=re.findall(r'\w+[\,]{0,1}\w+[\-]\w+[\,]{0,1}\w+:[0-9]',item3)
        for item4 in tmp2:
            if item4[1]!=',':
                num2[item4[0:-2]]=1
print(num2)

print("三糖如下:")
for item5 in data:
    if item5[0]=='t':
        tmp3=re.findall(r'\w+[\,]{0,1}\w+[\-]\w+[\,]{0,1}\w+[\-]\w+[\,]{0,1}\w+:[0-9]',item5)
        for item6 in tmp3:
            if item6[1]!=',':
                num3[item6[0:-2]]=1
print(num3)