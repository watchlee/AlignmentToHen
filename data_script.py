#!/usr/bin/env python
# coding=utf-8
import subprocess

input_pdb_list = []
file_path = './test_data_from_lose_same'
with open(file_path,'r') as file:
    for line in file:
        input_pdb_list.append(line.replace('\n',''))


origin_path = "~/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/" 
result_path_edit = "~/public_html/iPARTS2/result_edit_temp_data/"
result_path_alignment = "~/public_html/iPARTS2/result_temp_data/"
for line in input_pdb_list:
    print origin_path+line
    outpath_align=  result_path_alignment+line
    outpath_edit=  result_path_edit+line
    cmd_align = "mkdir "+outpath_align
    cmd_edit = "mkdir "+outpath_edit
    exe_cmd_align = './align '+origin_path+line+'/semi_input.php'+' '+outpath_align+'/result.php'
    exe_cmd_edit = './edit_iPARTS2_align '+origin_path+line+'/semi_input.php'+' '+outpath_edit+'/result.php'
    subprocess.call(cmd_align,shell=True) 
    subprocess.call(cmd_edit,shell=True) 
    subprocess.call(exe_cmd_align,shell=True) 
    subprocess.call(exe_cmd_edit,shell=True) 
    #print exe_cmd
    
