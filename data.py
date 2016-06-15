#!/usr/bin/env python
# coding=utf-8
list = []
list2 = []
average=[]
with open('./All_result','r') as file:
    for line in file:
        list.append(line.split(' ')[0]+' -'+line.split(' ')[2].replace('\n','')+'\n')
        list2.append(line.split(' ')[0]+' -'+line.split(' ')[3].replace('\n','')+'\n')

with open('RMSD_ROC','w') as file:
    for line in list:
        file.write(line)
with open('SAS_ROC','w') as file:
    for line in list2:
        file.write(line)

"""
list =[]
with open('./graduate1','r') as file:
    for line in file:
        list.append(line)

test_list = []
with open('./final_23C_4L_FSCOR_special_another','r') as file:
    for line in file:
        test_list.append(line)


for count in range(len(test_list)):
    for count2 in range(len(list)):
        if list[count2].split(' ')[1]==test_list[count].split(' ')[0]:
            print test_list[count].replace('\n','')#,list[count2]
"""
