#!/usr/bin/env python
# coding=utf-8
class PDB_DATA:
    def __init__(self,name,family,RMSD,SAS):
        self.name = name
        self.family = family
        self.RMSD = float(RMSD)
        self.SAS = float(SAS)
    def __repr__(self):
        return repr((self.name,self.family,self.RMSD,self.SAS))
    def setRMSD(self,RMSD):
        self.RMSD=RMSD
    def setSAS(self,SAS):
        self.SAS=SAS
    def getName(self):
        return self.name
    def getFamily(self):
        return self.family
    def getRMSD(self):
        return self.RMSD
    def getSAS(self):
        return self.SAS
    def getInfo(self):
        return self.name+'\t'+self.family,self.RMSD,self.SAS


def analysis_d2():
    
    sara_list=[]
    with open("./SARA_FSCOR_under1KK_2_SAS_another",'r') as file:
        for line in file:
            sara_list.append(line.split(' ')[1].replace('\n',''))
    family_list=[]
    iPARTS2_list= []
    with open("./final_23C_4L_FSCOR_2_SAS_another",'r') as file:
        for line in file:
            iPARTS2_list.append(line.split(' ')[1].replace('\n',''))
            family_list.append(line.split(' ')[0])    
    count=0
    number=0
    lost_list=[]
    with open('./graduate27','r') as file:
        for line in file:
            name= line.split(" ")[1]
            RMSD= line.split(" ")[2]
            SAS= line.split(" ")[3]
            lost_list.append(name+' '+family_list[count]+' '+RMSD+' '+SAS)
            count+=1
    sorted_sara_list=[]
    sorted_iPARTS2_list=[]
    sorted_iPARTS3_list=[]
    for count in range(len(lost_list)):
        if lost_list[count].split(' ')[3]!='NAN':
            #if lost_list[count].split(' ')[1]=='p,p':
            sorted_iPARTS3_list.append(PDB_DATA('',lost_list[count].split(' ')[1],0,lost_list[count].split(' ')[-1]))
            sorted_iPARTS2_list.append(PDB_DATA('',lost_list[count].split(' ')[1],0,iPARTS2_list[count].replace('-','')))
            sorted_sara_list.append(PDB_DATA('',lost_list[count].split(' ')[1],0,sara_list[count]))
    sorted_iPARTS3_list=sorted(sorted_iPARTS3_list,key=lambda PDB_DATA:PDB_DATA.SAS)
    sorted_iPARTS2_list=sorted(sorted_iPARTS2_list,key=lambda PDB_DATA:PDB_DATA.SAS)
    sorted_sara_list=sorted(sorted_sara_list,key=lambda PDB_DATA:PDB_DATA.SAS)

    with open('sorted_iPARTS2_under1K_d2_data','w') as file:
        for count in range(len(sorted_iPARTS2_list)): 
             file.write(sorted_iPARTS2_list[count].getFamily()+' '+str(sorted_iPARTS2_list[count].getSAS())+'\n')
    with open('sorted_iPARTS3_under1K_d2_data','w') as file:
        for count in range(len(sorted_iPARTS3_list)): 
             file.write(sorted_iPARTS3_list[count].getFamily()+' '+str(sorted_iPARTS3_list[count].getSAS())+'\n')
    with open('sorted_sara_under1K_d2_data','w') as file:
        for count in range(len(sorted_sara_list)): 
             file.write(sorted_sara_list[count].getFamily()+' '+str(sorted_sara_list[count].getSAS())+'\n')

     
            
            


        

def analysis():
    list = []

    lost_list = []
    sara_list = []
    sorted_sara_list= []
    sorted_iPARTS2_list=[]
    iPARTS2_origin_list=[]
    with open("./SARA_FSCOR_under1KK_0_SAS_another",'r') as file:
        for line in file:
            sorted_sara_list.append(PDB_DATA('',line.split(' ')[0],'0',line.split(' ')[1].replace('\n','')))
            sara_list.append(line.split(' ')[1].replace('\n',''))

    with open("./final_23C_4L_FSCOR_0_SAS_another",'r') as file:
        for line in file:
           iPARTS2_origin_list.append(line.split(' ')[1].replace('\n',''))
           sorted_iPARTS2_list.append(PDB_DATA('',line.split(' ')[0],'0',line.split(' ')[1].replace('\n','').replace('-','')))
    
    #sorted_sara_list=sorted(sorted_sara_list,key=lambda PDB_DATA:PDB_DATA.SAS)
    #sorted_iPARTS2_list=sorted(sorted_iPARTS2_list,key=lambda PDB_DATA:PDB_DATA.SAS)
    
    True_family_list =[]
    False_family_list=[]
    with open("./graduate30",'r') as file:
        for line in file:
            family = line.split(" ")[0]
            name= line.split(" ")[1]
            RMSD= line.split(" ")[2]
            SAS= line.split(" ")[3]

            lost_list.append(name+' '+family+' '+RMSD+' '+SAS)
            if family=='1':
                family="p,p"
                #True_family_list.append(name+' '+family)
            elif family=='0':
                family="n,p"
                #False_family_list.append(name+' '+family)
            if RMSD=="NAN":
                RMSD = 999.12345
                SAS = 999.123455
                
            #TP vs SAS使用
            list.append(PDB_DATA(name,family,RMSD,SAS))

    with open('sorted_iPARTS3_under1K_data','w') as file:
        for count in range(len(list)):
            #print list[count].getFamily()+' '+str(list[count].getSAS())
            file.write(list[count].getFamily()+' '+str(list[count].getSAS())+'\n')
    lost_analysis(lost_list,sara_list,iPARTS2_origin_list)
    '''
    write_file('True_under1K_list',True_family_list)
    write_file('False_under1K_list',False_family_list)
    '''

def lost_analysis(lost_list,sara_list,iPARTS2_origin_list):
    lost_sara_list=[]
    lost_iPARTS2_list=[]
    lost_my_list=[]
    lost_number=0
    average_SAS_SARA=0.0
    average_SAS_iPARTS2=0.0
    average_SAS_my=0.0
    average_true_SAS_SARA=0.0
    average_true_SAS_iPARTS2=0.0
    average_true_SAS_my=0.0
    average_false_SAS_SARA=0.0
    average_false_SAS_iPARTS2=0.0
    average_false_SAS_my=0.0
    lost_my_True_family_list =[]
    lost_my_False_family_list=[]
    lost_SARA_True_family_list =[]
    lost_SARA_False_family_list=[]
    lost_iPARTS2_True_family_list =[]
    lost_iPARTS2_False_family_list=[]

    analysis_lost_True_family_list=[]
    analysis_lost_False_family_list=[]
    
    count_win_number=0
    count_lose_number=0
    

    lost_SAS_list=[]
    
    true_count_win_number=0
    true_count_lose_number=0
    true_list=[]
    false_list=[]
    true_win_list=[]
    true_lose_list=[]
    false_count_win_number=0
    false_count_lose_number=0
    false_win_list=[]
    false_lose_list=[]
    lost_test_data=[]
    sorted_sara_list=[]
    sorted_iPARTS2_list=[]
    sorted_iPARTS3_list=[]
    
    lost_SAS_list.append('PDB vs PDB\tFamily\tmy SAS\tiPARTS2 SAS\tSARA SAS')
    for count in range(len(lost_list)):
        temp_string=lost_list[count].split(' ')[0]+' '+lost_list[count].split(' ')[1]+' '+lost_list[count].split(' ')[3].replace('\n','')
        if lost_list[count].find("NAN")!=-1:
            lost_number+=1
            lost_SAS_list.append(lost_list[count].split(' ')[0]+' '+lost_list[count].split(' ')[1]+' '+lost_list[count].split(' ')[3].replace('\n','')+' '+iPARTS2_origin_list[count].replace('-','')+' '+sara_list[count]) 
            lost_test_data.append(lost_list[count].split(' ')[0]) 
            #print lost_list[count]
        else:
            if lost_list[count].split(' ')[1]=='1':
                family='p,p'
            else:
                family='n,p'
            #print family,lost_list[count].split(' ')[3],iPARTS2_origin_list[count],sara_list[count]
            sorted_iPARTS3_list.append(PDB_DATA('',family,0,lost_list[count].split(' ')[3].replace('\n',''))) 
            sorted_iPARTS2_list.append(PDB_DATA('',family,0,iPARTS2_origin_list[count].replace('-',''))) 
            sorted_sara_list.append(PDB_DATA('',family,0,sara_list[count])) 
            
            iPARTS2_SAS= float(iPARTS2_origin_list[count].replace('-',''))
            my_SAS= float(lost_list[count].split(' ')[3].replace('\n',''))
            difference = iPARTS2_SAS-my_SAS
            
            if difference>=0.0:
                count_win_number+=1
            else:
                count_lose_number+=1
            #區分同family跟不同family
            if(lost_list[count].split(' ')[1]=='1'):
                if difference>=0.0:
                    true_count_win_number+=1
                    true_win_list.append(temp_string+' '+iPARTS2_origin_list[count].replace('-','')+' '+str(difference))
                else:
                    true_count_lose_number+=1
                    true_lose_list.append(temp_string+' '+iPARTS2_origin_list[count].replace('-','')+' '+str(difference))
                average_true_SAS_my+=my_SAS
                average_true_SAS_iPARTS2+=iPARTS2_SAS
                lost_my_True_family_list.append(lost_list[count])
                lost_SARA_True_family_list.append(sara_list[count])
                lost_iPARTS2_True_family_list.append(iPARTS2_origin_list[count])
                analysis_lost_True_family_list.append(temp_string+' '+sara_list[count].replace('\n','')+' '+iPARTS2_origin_list[count].replace('-','')+' '+str(difference))
                
                
            elif(lost_list[count].split(' ')[1]=='0'):
                if difference>=0.0:
                    false_count_win_number+=1
                    false_win_list.append(temp_string+' '+iPARTS2_origin_list[count].replace('-','')+' '+str(difference))
                else:
                    false_count_lose_number+=1
                    false_lose_list.append(temp_string+' '+iPARTS2_origin_list[count].replace('-','')+' '+str(difference))
                average_false_SAS_my+=my_SAS
                average_false_SAS_iPARTS2+=iPARTS2_SAS
                lost_my_False_family_list.append(lost_list[count])
                lost_SARA_False_family_list.append(sara_list[count])
                lost_iPARTS2_False_family_list.append(iPARTS2_origin_list[count])
                analysis_lost_False_family_list.append(temp_string+' '+sara_list[count].replace('\n','')+' '+iPARTS2_origin_list[count].split('-')[0]+' '+str(difference))


            average_SAS_iPARTS2+=float(iPARTS2_origin_list[count])
            average_SAS_my+=float(lost_list[count].split(' ')[3])
            average_SAS_SARA+=float(sara_list[count])

            

            lost_my_list.append(lost_list[count])
            lost_iPARTS2_list.append(iPARTS2_origin_list[count])
            lost_sara_list.append(lost_list[count].split(' ')[1]+' '+sara_list[count])
    sorted_iPARTS2_list=sorted(sorted_iPARTS2_list,key=lambda PDB_DATA:PDB_DATA.SAS) 
    sorted_iPARTS3_list=sorted(sorted_iPARTS3_list,key=lambda PDB_DATA:PDB_DATA.SAS) 
    sorted_sara_list=sorted(sorted_sara_list,key=lambda PDB_DATA:PDB_DATA.SAS) 
    with open('sorted_iPARTS2_under1K_data','w') as file:
        for count in range(len(sorted_iPARTS2_list)):
            file.write(sorted_iPARTS2_list[count].getFamily()+' '+str(sorted_iPARTS2_list[count].getSAS())+'\n')
    with open('sorted_iPARTS3_under1K_data','w') as file:
        for count in range(len(sorted_iPARTS3_list)):
            file.write(sorted_iPARTS3_list[count].getFamily()+' '+str(sorted_iPARTS3_list[count].getSAS())+'\n')
    with open('sorted_sara_under1K_data','w') as file:
        for count in range(len(sorted_sara_list)):
            file.write(sorted_sara_list[count].getFamily()+' '+str(sorted_sara_list[count].getSAS())+'\n')
    print lost_number
    print 'my total:',len(lost_my_list),' SARA total:',len(lost_sara_list),' iPARTS2 total:',len(lost_iPARTS2_list)
    print 'my average SAS=',average_SAS_my/len(lost_list),' SARA average SAS=',average_SAS_SARA/len(lost_sara_list),' iPARTS2 average SAS=',abs(average_SAS_iPARTS2/len(lost_iPARTS2_list))
    print count_win_number,count_lose_number
    print 'true win:',true_count_win_number,'true lose:',true_count_lose_number
    print 'false win:',false_count_win_number,'false lose:',false_count_lose_number
    print 'true average:',average_true_SAS_my/(true_count_win_number+true_count_lose_number),' true average iPARTS2',average_true_SAS_iPARTS2/(true_count_win_number+true_count_lose_number)
    print 'false average:',average_false_SAS_my/(false_count_win_number+false_count_lose_number),' true average iPARTS2',average_false_SAS_iPARTS2/(false_count_win_number+false_count_lose_number)
    write_file('lost_true_analysis',analysis_lost_True_family_list)
    write_file('lost_false_analysis',analysis_lost_False_family_list)
    write_file('true_win_analysis',true_win_list)
    write_file('true_lose_analysis',true_lose_list)
    write_file('false_win_analysis',false_win_list)
    write_file('false_lose_analysis',false_lose_list)
    write_file('lost_SAS_list',lost_SAS_list)
    write_file('lost_test_data',lost_test_data)

    with open('iPARTS2_roc','w') as file:
        for line in lost_my_list:
            file.write(line.split(' ')[1]+' -'+line.split(' ')[3].replace('\n','')+'\n')


def write_file(name,list):
    with open(name,'w') as file:
        for line in list:
            file.write(line+'\n')



analysis()
#analysis_d2()
