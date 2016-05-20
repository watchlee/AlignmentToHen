#!/usr/bin/env python
# coding=utf-8
import subprocess

input_pdb_list = []
file_path = './test_data_from_lose_same'
#file_path = './test_data_for_correct_alignment'
with open(file_path,'r') as file:
    for line in file:
        input_pdb_list.append(line.replace('\n',''))


origin_path = "~/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/" 

#edit distance version
result_path_edit = "~/public_html/iPARTS2/result_edit_temp_data/"

#score version
#result_path_alignment = "~/public_html/iPARTS2/result_score_temp_data/"

number = 5
k = 2
l = 0.1
m = 1.5
result_path_alignment = "~/public_html/result_score_data_N"+str(number)+'_K'+str(k)+'_L'+str(l)+'_M'+str(m)+'/'
subprocess.call('mkdir '+result_path_alignment,shell=True)

for line in input_pdb_list:
    print origin_path+line
    outpath_align=  result_path_alignment+line
    outpath_edit=  result_path_edit+line
    cmd_align = "mkdir "+outpath_align
    cmd_edit = "mkdir "+outpath_edit
    exe_cmd_align = './align '+origin_path+line+'/semi_input.php'+' '+outpath_align+'/result.php'+' '+str(number)+' '+str(k)+' '+str(l)+' '+str(m)
#    exe_cmd_edit = './edit_iPARTS2_align '+origin_path+line+'/semi_input.php'+' '+outpath_edit+'/result.php'
    subprocess.call(cmd_align,shell=True) 
 #   subprocess.call(cmd_edit,shell=True) 
    subprocess.call(exe_cmd_align,shell=True) 
  #  subprocess.call(exe_cmd_edit,shell=True) 
    #print exe_cmd
    
