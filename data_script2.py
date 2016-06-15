#!/usr/bin/env python
# coding=utf-8
import subprocess

input_pdb_list = []
#file_path = './test_Data_for_SARA_long'
#file_path = './test_Data_for_SARA'

#file_path = './debug_list'

#file_path = './True_under1K_list'
file_path = './SARA_FSCOR_under1K_alltoall'
#file_path = './ttttt'

#file_path = './debug_Data'
#file_path = './test_data_from_lose_same'
#file_path = './test_data_for_correct_alignment'
family_list = []
with open(file_path,'r') as file:
    for line in file:
        input_pdb_list.append(line.split(' ')[0].replace('\n',''))
        family_list.append(line.split(' ')[1].replace('\n',''))


origin_path = "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/" 

#edit distance version
result_path_edit = "~/public_html/iPARTS2/result_edit_temp_data/"
result_list= []


class PDB_DATA:
    def __init__(self,name,family,RMSD,SAS):
        self.name = name
        self.family = family
        self.RMSD = RMSD
        self.SAS = SAS
    def getName(self):
        return self.name
    def getFamily(self):
        return self.family
    def getRMSD(self):
        return self.RMSD
    def getSAS(self):
        return self.SAS

def loop_test():
    count=5
    count2=4.5
    count3=4
    count4=0.05
    count5=0.01
    open=9
    exten=1.0
    w_b=0.1
    pdb_name="1N66_A_to_4TRA_A"
    pdb_name="1FHK_A_to_1JBT_C"
    pdb_name="1FHK_A_to_1Q93_A"
    pdb_name="1AM0_A_to_1QVF_3"
    pdb_name="1N66_A_to_4TRA_A"

    list =["1AM0_A_to_1QVF_3","1AM0_A_to_1JO7_A","1AM0_A_to_1M82_A","1AM0_A_to_1BZU_A","1AM0_A_to_1BZT_A","1AM0_A_to_1BZ3_A","1AM0_A_to_1BZ2_A","1N66_A_to_4TRA_A"]
    #exe_cmd_align = './semi_affine_version_align '+origin_path+pdb_name+'/semi_input.php'+' ~/result.php'+' ~/error.php '+str(count)+' '+str(count2)+' '+str(count3)+' '+str(count4)+' '+str(count5)+' '+str(-w_b)+' '+str(open)+' '+str(exten)
    #subprocess.call(exe_cmd_align,shell=True) 
    #profit_process(pdb_name)
    for name in list: 
        exe_cmd_align = './semi_affine_version_align '+origin_path+name+'/semi_input.php'+' ~/result.php'+' ~/error.php '+str(count)+' '+str(count2)+' '+str(count3)+' '+str(count4)+' '+str(count5)+' '+str(-w_b)+' '+str(open)+' '+str(exten)
        subprocess.call(exe_cmd_align,shell=True) 
        profit_process(name)


    test = 0.0
    for x in range(10):
        print "loop=",x,test
        exe_cmd_align = './semi_affine_version_align '+origin_path+pdb_name+'/semi_input.php'+' ~/result.php'+' ~/error.php '+str(x)+' '+str(count2)+' '+str(count3)+' '+str(count4)+' '+str(count5)+' '+str(-w_b)+' '+str(open)+' '+str(exten)
        subprocess.call(exe_cmd_align,shell=True) 
        profit_process(pdb_name)
        test+=1

    '''

    with open('Loop_result','w') as file:
        for count in range(len(result_list)):
            file.write(result_list[count]+'\n')
    '''
###class list
pdb_list=[]
def profit_process(path,profit_result_path,result_path):
    pdb_list= path.split('_to_')
    seq1=''
    seq2=''
    with open(profit_result_path,'r') as file:
        seq1=file.readline().replace('\n','')
        seq2=file.readline().replace('\n','')
    new_align_number=0
    for count in range(len(seq1)):
        if seq1[count]!='-' and seq2[count]!='-':
            new_align_number+=1

    with open(result_path,'w') as file:
        file.write(">P1;"+pdb_list[0]+'\n')
        file.write("1234567890\n")
        file.write(seq1+'*\n')
        file.write(">P2;"+pdb_list[1]+'\n')
        file.write("1234567890\n")
        file.write(seq2+'*\n')
    reference = 'reference /home/watchlee/Research_Programming/RMSD/pdb/'+pdb_list[0]+'.pdb\n'
    mobilebl = 'mobile /home/watchlee/Research_Programming/RMSD/pdb/'+pdb_list[1]+'.pdb\n'
    readalignment= "readalignment "+result_path+'\n'
    atoms="atoms *\n"
    ignoremissing="ignoremissing\n"
    fit = "fit\n"
    quit ="quit\n"

    pf_script_path="/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/"+path+'/pf_script'
    profit_log_path="/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/"+path+'/profit_log0'
    print pf_script_path,profit_log_path
    with open(pf_script_path,'w') as file:
        file.write(reference)
        file.write(mobilebl)
        file.write(readalignment)
        file.write(atoms)
        file.write(ignoremissing)
        file.write(fit)
        file.write(quit)
    subprocess.call("./profit2.5.3 < "+pf_script_path+" > "+profit_log_path,shell=True)
    RMSD = "NAN" 
    with open(profit_log_path,'r') as file:
        for line in file:
            if(line.find('RMS:')!=-1):
                RMSD=line.split('RMS: ')[1].replace('\n','')
                print line.split('RMS: ')[1].replace('\n','')
    print 'RMSD',RMSD

    if(RMSD!="NAN" and float(new_align_number)!=0):
        SAS=float(RMSD)*100/float(new_align_number)
    else:
        SAS="NAN"
    print 'SAS',SAS
    result_list.append(pdb_list[0]+'_to_'+pdb_list[1]+' '+str(RMSD)+' '+str(SAS)+' '+str(new_align_number))
#score version
#result_path_alignment = "~/public_html/iPARTS2/result_score_temp_data/"













'''
((((---....)------)))
uuuZ---3rCCu------ZZ) 21
-uu-BXZu----9NIvu)ZZ- 21
-((-(((.----...)))))-
'''
number =1 
arc_match_multiple=4
arc_mismatch_multiple=0.5
arc_halfmatch_multiple=1.5
arc_mismatch_case1=3
arc_mismatch_case2=1
arc_mismatch_case3=1
arc_mismatch_case4=1
arc_breaking =0.5
arc_breaking2 =1.5
common_opp=9
common_exp = 1

#linaer gap penalty path
#result_path_alignment = "~/public_html/result_score_data_N"+str(number)+'_K'+str(arc_match_multiple)+'_L'+str(arc_mismatch_multiple)+'_M'+str(arc_halfmatch_multiple)+'/'
result_path_alignment = "~/public_html/result_score_data_B"+str(arc_breaking)+"O"+str(common_opp)+"E"+str(common_exp)+"_K"+str(arc_match_multiple)+'_L'+str(arc_mismatch_multiple)+'_M'+str(arc_halfmatch_multiple)+"/"
subprocess.call('mkdir '+result_path_alignment,shell=True)


def function():
    for line in input_pdb_list:
        print origin_path+line
        outpath_align=  result_path_alignment+line
        outpath_edit=  result_path_edit+line
        cmd_align = "mkdir "+outpath_align
        cmd_edit = "mkdir "+outpath_edit
        profit_result_path = origin_path+line+'/profit_result'
        result_path = origin_path+line+'/profit_data'
        #exe_cmd_align = './semi_affine_version_align '+origin_path+line+'/semi_input.php'+' '+outpath_align+'/result.php'+' '+outpath_align+'/error.php '+str(arc_match_multiple)+' '+str(arc_mismatch_case1)+' '+str(arc_mismatch_case2)+' '+str(arc_mismatch_case3)+' '+str(arc_mismatch_case4)+' '+str(arc_breaking)+' '+str(common_opp)+' '+str(common_exp)+' '+profit_result_path
        exe_cmd_align = './iPARTS3_cost '+origin_path+line+'/semi_input.php'+' '+outpath_align+'/result.php'+' '+outpath_align+'/error.php '+str(arc_match_multiple)+' '+str(arc_mismatch_case1)+' '+str(arc_mismatch_case2)+' '+str(arc_mismatch_case3)+' '+str(arc_mismatch_case4)+' '+str(arc_breaking)+' '+str(arc_breaking2)+' '+str(common_opp)+' '+str(common_exp)+' '+profit_result_path
    #    exe_cmd_align = './semi_global '+origin_path+line+'/semi_input.php'+' '+outpath_align+'/result.php'+' '+str(arc_match_multiple)+' '+str(arc_mismatch_case1)+' '+str(arc_mismatch_case2)+' '+str(arc_mismatch_case3)+' '+str(arc_mismatch_case4)+' '+str(arc_breaking)+' '+str(common_opp)+' '+str(common_exp)
        #exe_cmd_align = './affine_version_align '+origin_path+line+'/semi_input.php'+' '+outpath_align+'/result.php'+' '+outpath_align+'/error.php '+str(arc_match_multiple)+' '+str(arc_mismatch_case1)+' '+str(arc_mismatch_case2)+' '+str(arc_mismatch_case3)+' '+str(arc_mismatch_case4)+' '+str(arc_breaking)+' '+str(common_opp)+' '+str(common_exp)
        print exe_cmd_align
        #exe_cmd_align = './align '+origin_path+line+'/semi_input.php'+' '+outpath_align+'/result.php'+' '+str(number)+' '+str(k)+' '+str(l)+' '+str(m)+' '+str(base_open)+' '+str(base_exp)+' '+str(arc_open)+' '+str(arc_exp)
    #    exe_cmd_edit = './edit_iPARTS2_align '+origin_path+line+'/semi_input.php'+' '+outpath_edit+'/result.php'
        subprocess.call(cmd_align,shell=True) 
     #   subprocess.call(cmd_edit,shell=True) 
        subprocess.call(exe_cmd_align,shell=True) 
        profit_process(line,profit_result_path,result_path)
      #  subprocess.call(exe_cmd_edit,shell=True) 
        #print exe_cmd
    for count in range(len(result_list)):
        family = family_list[count]
        name = result_list[count].split(' ')[0]
        RMSD = result_list[count].split(' ')[1]
        SAS = result_list[count].split(' ')[2]
        pdb=PDB_DATA(name,family,RMSD,SAS)
        pdb_list.append(pdb)



    with open('True_result','w') as file:
        for count in range(len(result_list)):
            file.write(family_list[count]+' '+result_list[count]+'\n')

#loop_test()
function()
