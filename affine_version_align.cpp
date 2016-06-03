/*************************************************************************
> File Name: align.cpp
> Author: 
> Mail: 
> Created Time: Thu Mar 24 11:35:41 2016
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <stack>
using namespace std;
#define TYPE_STRUCTURE
#define TYPE2_STRUCTURE
#define TYPE3_STRUCTURE
//#define _DEBUG 
//要debug的時候請把註解去掉
//#define _TRACE_DEBUG
/*alignment score*/
struct alignment{
    int p1,p2;
    double weight;
    struct alignment *next;
};
/*紀錄位置*/
struct four_tuple
{
    int l1;
    int r1;
    int l2;
    int r2;
};
/*Scoring matrix data structure*/
typedef struct index_matrix{
    char alphabet;
    int pos;
}IndexMatrix;
/*Scoring matrix merge sort*/
class Merge_Sort
{
    public:
    static void Sort(IndexMatrix* array,int length);
    private:
    static void Sort(IndexMatrix* array, IndexMatrix* temp,int length,int start,int count);
    static void Merge(IndexMatrix* array, IndexMatrix* temp,int length,int leftstart,int leftcount,int rightstart,int rightcount);
};
void Merge_Sort::Sort(IndexMatrix* array,int length)
{
    IndexMatrix *temp = new IndexMatrix[length];
    Sort(array,temp,length,0,length);
}
void Merge_Sort::Sort(IndexMatrix* array,IndexMatrix* temp,int length,int start,int count)
{
    /*太小不必比較*/
    if(count<2)
    return ;
    Sort(array,temp,length,start,count/2);
    Sort(array,temp,length,start+count/2,count-count/2);
    Merge(array,temp,length,start,count/2,start+count/2,count-count/2);
}
void Merge_Sort::Merge(IndexMatrix* array, IndexMatrix* temp,int length,int leftstart,int leftcount,int rightstart,int rightcount)
{
    int i = leftstart, j = rightstart, leftbound = leftstart+leftcount, rightbound = rightstart+rightcount, index=leftstart;
    while(i<leftbound||j<rightbound)
    {
        if(i< leftbound&& j < rightbound)
        {
            int val,val2;
            val = array[j].alphabet;
            val2 = array[i].alphabet;
            if(val < val2)
            {
                temp[index].alphabet=array[j].alphabet;
                temp[index].pos=array[j++].pos;
            }       
            else
            {
                temp[index].alphabet=array[i].alphabet;
                temp[index].pos=array[i++].pos;
            }
        }
        else if(i<leftbound)
        {
            temp[index].alphabet=array[i].alphabet;
            temp[index].pos=array[i++].pos;
        }
        else
        {
            temp[index].alphabet=array[j].alphabet;
            temp[index].pos=array[j++].pos;
        }
        ++index;
    }
    for(i = leftstart;i<index;++i)
    {
        array[i].alphabet=temp[i].alphabet;
        array[i].pos=temp[i].pos;
    }
}
/*-------------------------------------------------------------------------------------*/
/*------------------------------------global variable-----------------------------------*/
IndexMatrix *alphabet_index=NULL;
static int size=0;
static int gap_opp,gap_exp;//Input gap penalty arguments
static string seq1,arc1,seq2,arc2,matrixpath;
static string aseq1,aseq2,astr1,astr2;
vector<double> weights;
vector<int>    L1,R1,I1,L2,R2,I2;
stack<int>     str_stack;
vector<vector<double> > M;
vector<vector<double> > lower;
vector<vector<double> > upper;
vector<vector<double> > middle;
/*記錄矩陣各個entry來的方向*/
vector<vector<double> > direct;
vector<vector<double> > D;
const double eps=0.0000001;
static int **scoring_matrix;
double s_d,s_aa,s_am,s_breaking,s_altering;
//similarity
//arc match = 4 * w_a;
double MSP_alignment_value = 0.0;
int not_free1    (int pos)          { return (arc1[pos]=='.' ? 0:1)      ;  }
int not_free2    (int pos)          { return (arc2[pos]=='.' ? 0:1)      ;  }
/*沒這麼簡單了*/
int arc_mismatch(int pos1,int pos2){ return (seq1[pos1]!=seq2[pos2]?1:0);  }
int base_matching(int,int);
double arc_match_weight = 4;//arc match 
double arc_mismatch_weight= 1;//arc mismatch
double arc_halfmatch_weight= 3;//when one base score >=0 and other base score < 0, m 不能比k大
double number=1;
double w_d =-1*number;  // base deletion
double w_r =-2*number;  // arc  removing
//double w_b =-1.5*number;  // arc  breaking
double w_b =-0;  // arc  breaking
double common_opp = 9;
double common_exp = 1;
double deletion_cost = 0.0;
double match_cost= 0.0;
/*
double w_am=-1.5;  // arc  mismatch
double w_aa=0;     // arc  match
*/
double w_m=1; //base mismatch
double arc_operation(int p1,int p2,int p3,int p4)
{
    //當alignment score > = 0 views as arc match
    //if(base_matching(p1,p2)>=0&&base_matching(p3,p4)>=0)
    //when both base score > 0
    if(base_matching(p1,p2)>=0&&base_matching(p3,p4)>=0)
    {
        //return 0;
        if(base_matching(p1,p2)==0&&base_matching(p3,p4)==0)
            return arc_match_weight;
        else
            return arc_match_weight*(base_matching(p1,p2)+base_matching(p3,p4));

    }
    //當alignment score <0 views as arc-mismatch
    else
    {
        if((base_matching(p1,p2)<0 && base_matching(p3,p4)>=0)||(base_matching(p1,p2)>=0&&base_matching(p3,p4)<0))
        {
            if((base_matching(p1,p2)+base_matching(p3,p4))>=0)
            {
                return (base_matching(p1,p2)+base_matching(p3,p4))*arc_halfmatch_weight;
            } 
            else
            {
                return (base_matching(p1,p2)+base_matching(p3,p4))*(arc_mismatch_weight/2);
            }
        }
        else if(base_matching(p1,p2)<=0 && base_matching(p3,p4)<0)
        //return (base_matching(p1,p2)+base_matching(p3,p4))*0.5*w_am;
        return (base_matching(p1,p2)+base_matching(p3,p4))*arc_mismatch_weight;
    }
}
int BinarySearch(char ,IndexMatrix *,int );
void traceback();
void insert(alignment*,int,int,double);
int** readmat(char *);
void read_data(const char *);
void write_data(double,const char *);
void write_error(double,double,const char *);
double computation();
double test_alignment();
string special_character_processing(string);
//做一般的global alignment
/*需要寫一個python 去處理input file data*/
/* input format
*$seq1=""
*$arc1=""
*$seq2=""
*$arc2=""
*$opp=""
*$exp=""
*
*/
double max(double a,double b)
{
    return ((a<b)? b:a);
}
double max(double a,double b,double c)
{
    double value=(a<b)?b:a;
    return ((value<c)? c:value);
}
double max4(double a,double b,double c,double d)
{
    double value=(a<b)? b:a;
    double value2=(c<d)? d:c;
    return ((value<value2)? value2:value);
    /*
    if (a>=b && a>=c && a>=d)
    return a;
    else if (b>=a && b>=c && b>=d)
    return b;
    else if (c>=a && c>=b && c>=d)
    return c;
    return d;
    */
}
double min4(double a,double b,double c,double d)
{
    if (a<=b && a<=c && a<=d)
    return a;
    else if (b<=a && b<=c && b<=d)
    return b;
    else if (c<=a && c<=b && c<=d)
    return c;
    return d;
}
int main(int argc,char* argv[])
{
    const char *pdb_compare_path ;
    const char *pdb_result;
    const char *error_result;
    cout<<argc<<endl;
    if(argc!=10)
    {
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/test_data/1A9N_Q_to_1E7K_C/semi_input.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1M90_B_to_1NKW_9/semi_input.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1IBK_A_to_1N33_A/semi_input.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1FG0_A_to_1BZ3_A/semi_input.php";
       //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1D0T_A_to_1ZDK_R/semi_input.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1BAU_A_to_1S9S_A/semi_input.php";
       //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1FHK_A_to_1S9S_A/semi_input.php";
       //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1MT4_A_to_1I3Y_A/semi_input.php";
       //特殊例子
      // pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1AM0_A_to_1FMN_A/semi_input.php";
       //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1FHK_A_to_1ZIF_A/semi_input.php";
       //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1FHK_A_to_1BYJ_A/semi_input.php";
       //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1BN0_A_to_1AQ3_R/semi_input.php";
       //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1Q9A_A_to_1QA6_C/semi_input.php";
       pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1J4Y_A_to_1Q2S_E/semi_input.php";
       //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1Q9A_A_to_1HC8_C/semi_input.php";
       //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1G70_A_to_1M5L_A/semi_input.php";
      // pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1UN6_E_to_1M90_B/semi_input.php";
       //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/5MSF_S_to_7MSF_R/semi_input.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/AlignmentToHen/Test_for_GLOBAL/1FG0_A_to_1BZ3_A.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/AlignmentToHen/Test_for_GLOBAL/1D0T_A_to_1ZDK_R.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/AlignmentToHen/Test_for_GLOBAL/1FHK_A_to_1FKA_A.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/AlignmentToHen/Test_for_GLOBAL/1FHK_A_to_1KC8_A.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/AlignmentToHen/Test_for_GLOBAL/1FHK_A_to_1KC8_A.php";
       //pdb_compare_path= "/home/watchlee/Research_Programming/AlignmentToHen/Test_for_GLOBAL/1MT4_A_to_1HC8_C.php";
       //pdb_compare_path= "/home/watchlee/Research_Programming/AlignmentToHen/Test_for_GLOBAL/1G70_A_to_1M5L_A.php";
        pdb_result= "/home/watchlee/result.php";
        error_result= "/home/watchlee/error.php";
        number = 1;
        arc_match_weight = 3;
        arc_mismatch_weight = 1;
        arc_halfmatch_weight = 1.5;
    }
    else
    {
        pdb_compare_path=argv[1] ;
        pdb_result= argv[2];
        error_result= argv[3];

        number=1;
        arc_match_weight = atof(argv[4]);
        arc_mismatch_weight = atof(argv[5]);
        arc_halfmatch_weight = atof(argv[6]); 
        w_b =atof(argv[7]);
        

        //k = atof(argv[4]);
        //l = atof(argv[5]);
        //m = atof(argv[6]);
        common_opp=atof(argv[8]);
        common_exp=atof(argv[9]);
        /*
        base_opp = atoi(argv[7]);
        base_exp = atoi(argv[8]);
        arc_opp = atoi(argv[9]);
        arc_exp = atoi(argv[10]);
        */
    }
    char path[100];
    //-------------Test Scoring Matrix
    //sprintf(path,"./SM/BLOSUM-like_scoring_matrix");
    //sprintf(path,"./SM/iPARTS2_new_23C_4L_matrix");
    sprintf(path,"./SM/23-4L_matrix");
    //sprintf(path,"./SM/SARA_23C_4L_matrix");
    scoring_matrix=readmat(path);
    //cout<<"arc max = "<<arc_max<<endl;
    //cout<<"seq max = "<<seq_max<<endl;
    int count = 0;
    #ifdef _DEBUG 
    char test_char[100];
    for(count = 0;count<size;count++)
    {
        test_char[count]=alphabet_index[count].alphabet;
        printf("%d %c %d\n",alphabet_index[count].alphabet,alphabet_index[count].alphabet,alphabet_index[count].pos);
    }
    for(count = 0;count<size;count++)
    printf("%c ",test_char[count]);
    printf("\n");
    #endif
    Merge_Sort sorting;
    sorting.Sort(alphabet_index,size);
    /*取得處理資料需要的data*/
    //read_data("./test_file_2","./result"); 
    //read_data("./test_file3","./result 
    //const char *pdb_compare= "/home/watchlee/Research_Programming/X3DNA/test_data/2FMT_C_to_1J2B_D/semi_input.php";
    read_data(pdb_compare_path); 
    //read_data("./test_file4","./result"); 
    //read_data("./test_file5","./result"); 
    //read_data("./test_file","./result"); 
    /*Computation*/
    cout<<arc1<<endl;
    cout<<seq1<<endl;
    cout<<seq2<<endl;
    cout<<arc2<<endl;
    double score = computation();
    traceback();
    cout<<"Score="<<score<<endl;
    double total = 0.0;
    for(int count = 0;count<weights.size();count++)
    {
        //cout<<weights[count]<<" ";
        total+=weights[count];
    }

    cout<<astr1<<endl;
    cout<<aseq1<<" "<<aseq1.size()<<endl;
    cout<<aseq2<<" "<<aseq2.size()<<endl;
    cout<<astr2<<endl;
    cout<<total<<endl;
    cout<<"path="<<pdb_compare_path<<endl;
    cout<<"check out here! "<<score<<" "<<total<<endl;
    //cout<<"score="<<test_alignment()<<endl;
    if(abs(double(score-total))<0.01)
    {
        cout<<"correct!"<<endl;
        write_data(score,pdb_result);
    }
    else
    {
        cout<<"incorrect!"<<endl;
        write_error(score,total,error_result);
    }

    /*釋放記憶體*/
    free(alphabet_index);//來源 line:33
return 0;
}
/*判斷方向*/
string determingDIRECT(double temp_middle,double temp_lower,double temp_upper,double current_score)
{
    string determing ;
    if(temp_upper==current_score) 
        determing="7";
    if(temp_lower==current_score)
        determing="6";
    if(temp_middle==current_score)
        determing="5";
    if(temp_lower==current_score&&temp_upper==current_score)
        determing="4";
    if(temp_middle==current_score&&temp_upper==current_score)
        determing="3";
    if(temp_middle==current_score&&temp_lower==current_score)
        determing="2";
    if(temp_middle==current_score&&temp_lower==current_score&&temp_upper==current_score)
        determing="1";
    return determing; 
}
/*------------affine gap peanlty alignemnt-----------------------*/
double test_alignment()
{
    int gap_opp = 9;
    int gap_exp = 1;
    double LARGE_NUMBER = 999999;
    vector<vector<double> > temp_D;
    vector<vector<double> > temp_middle;
    vector<vector<double> > temp_lower;
    vector<vector<double> > temp_upper;
    vector<vector<string> >direct_middle;
    vector<vector<string> >direct_lower;
    vector<vector<string> >direct_upper;
    vector<vector<string> >direct;
    temp_D.resize(arc1.size()+1);
    temp_middle.resize(arc1.size()+1);
    temp_lower.resize(arc1.size()+1);
    temp_upper.resize(arc1.size()+1);
    direct_middle.resize(arc1.size()+1);
    direct_lower.resize(arc1.size()+1);
    direct_upper.resize(arc1.size()+1);
    direct.resize(arc1.size()+1);
    for(int count = 0;count<=arc1.size();count++)
    {
        temp_D[count].resize(arc2.size()+1);
        temp_middle[count].resize(arc2.size()+1);
        temp_lower[count].resize(arc2.size()+1);
        temp_upper[count].resize(arc2.size()+1);
        direct_middle[count].resize(arc2.size()+1);
        direct_lower[count].resize(arc2.size()+1);
        direct_upper[count].resize(arc2.size()+1);
        direct[count].resize(arc2.size()+1);
    }
    temp_middle[0][0]=temp_lower[0][0]=temp_upper[0][0]=temp_D[0][0]=0.0;
    direct[0][0]="o";
    /*某些部分在affine gap penalty中entry不會使用到設定為負無限大當作不會參考*/
    /*對Y軸進行初始*/
    for(int count = 1; count<=arc1.size();count++)
    {
        temp_lower[count][0]=-LARGE_NUMBER;
        temp_middle[count][0]=-LARGE_NUMBER;
        temp_D[count][0]=temp_upper[count][0]=-gap_opp-(count*gap_exp);
        direct_middle[count][0]=direct_lower[count][0]=direct_upper[count][0]=direct[count][0]="B";
    }
    /*對X軸進行初始*/
    for(int count = 1; count<=arc2.size();count++)
    {
        temp_middle[0][count]=-LARGE_NUMBER;
        temp_upper[0][count]=-LARGE_NUMBER;
        temp_D[0][count]=temp_lower[0][count]=-gap_opp-(count*gap_exp);
        direct_middle[0][count]=direct_lower[0][count]=direct_upper[0][count]=direct[0][count]="C";
    }
    for(int count = 1;count<=arc1.size();count++)
    {
        for(int count2 = 1;count2<=arc2.size();count2++)
        {
            double temp_middles,temp_lowers,temp_uppers;
            double temp_value=max(temp_middle[count-1][count2-1],temp_lower[count-1][count2-1],temp_upper[count-1][count2-1]);
            temp_middle[count][count2]=temp_value+base_matching(count-1,count2-1);
            temp_middles=temp_middle[count-1][count2-1];
            temp_lowers=temp_lower[count-1][count2-1];
            temp_uppers=temp_upper[count-1][count2-1];
            string temp_direct_middle = determingDIRECT(temp_middles,temp_lowers,temp_uppers,temp_value);
            direct_middle[count][count2]=temp_direct_middle;
            temp_value=max(temp_lower[count-1][count2]-gap_exp,temp_middle[count-1][count2]-gap_opp-gap_exp,temp_upper[count-1][count2]-gap_opp-gap_exp);
            temp_lower[count][count2]=temp_value;
            temp_middles=temp_middle[count-1][count2]-gap_opp-gap_exp;
            temp_lowers=temp_lower[count-1][count2]-gap_exp;
            temp_uppers=temp_upper[count-1][count2]-gap_opp-gap_exp;
            string temp_direct_lower=determingDIRECT(temp_middles,temp_lowers,temp_uppers,temp_value);
            direct_lower[count][count2]=temp_direct_lower;
            temp_value=max(temp_upper[count][count2-1]-gap_exp,temp_middle[count][count2-1]-gap_opp-gap_exp,temp_lower[count][count2-1]-gap_opp-gap_exp);
            temp_upper[count][count2]=temp_value;
            temp_middles=temp_middle[count][count2-1]-gap_opp-gap_exp;
            temp_lowers=temp_lower[count][count2-1]-gap_opp-gap_exp;
            temp_uppers=temp_upper[count][count2-1]-gap_exp;
            string temp_direct_upper=determingDIRECT(temp_middles,temp_lowers,temp_uppers,temp_value);
            direct_upper[count][count2]=temp_direct_upper;
            //arc annotated sequence must consider four tuples lower,upper,middle,special
            temp_D[count][count2]=max(temp_middle[count][count2],temp_lower[count][count2],temp_upper[count][count2]);
            string temp_direct;
            if(temp_upper[count][count2]==temp_D[count][count2])
                temp_direct="C"; 
            if(temp_lower[count][count2]==temp_D[count][count2])
                temp_direct="B";
            if(temp_middle[count][count2]==temp_D[count][count2])
                temp_direct="A";
            direct[count][count2]=temp_direct;
        }
    }
    /*清除記憶體 http://stackoverflow.com/questions/10464992/c-delete-vector-objects-free-memory*/
    vector<vector<double> >().swap(temp_middle);
    vector<vector<double> >().swap(temp_lower);
    vector<vector<double> >().swap(temp_upper);
    vector<vector<string> > D = direct;
    cout<<"為跑之前"<<endl;
    for(int count = 0;count<=arc1.size();count++)
    {
        for(int count2 = 0;count2<=arc2.size();count2++)
        {
            cout<<D[count][count2]<<" ";
        }
        cout<<endl;
    }
    int i = arc1.size(),j = arc2.size(); 
    string s1="",s2="";
    while(i+j!=0)
    {
        if(D[i][j]=="A")
        {
            if(direct_middle[i][j]=="5"||direct_middle[i][j]=="1"||direct_middle[i][j]=="2"||direct_middle[i][j]=="3")
                D[i-1][j-1]="A"; 
            if(direct_middle[i][j]=="6"||direct_middle[i][j]=="4"&&j>1&&i>1)
                D[i-1][j-1]="B";
            if(direct_middle[i][j]=="7"&&j>1&&i>1)
                D[i-1][j-1]="C";
            s1.insert(s1.begin(),1,seq1[i-1]);
            s2.insert(s2.begin(),1,seq2[j-1]);
            i--;
            j--;
        }
        if(D[i][j]=="B")
        {
            if(direct_lower[i][j]=="5"||direct_lower[i][j]=="1"||direct_lower[i][j]=="2"||direct_lower[i][j]=="3")
                D[i-1][j]="A";
            if(direct_lower[i][j]=="6"||direct_lower[i][j]=="4"&&i>1)
                D[i-1][j]="B";
            if(direct[i][j]=="7"&&j>1)
                D[i-1][j]="C";
            s1.insert(s1.begin(),1,seq1[i-1]);
            i--;
            s2.insert(s2.begin(),1,'-');
        }
        if(D[i][j]=="C")
        {
            if(direct_upper[i][j]=="5"||direct_upper[i][j]=="1"||direct_upper[i][j]=="2"||direct_upper[i][j]=="3")
                D[i][j-1]="A";
            if(direct_upper[i][j]=="6"||direct_upper[i][j]=="4"&&j>1)
                D[i][j-1]="B";
            if(direct_upper[i][j]=="7"&&j>1)
                D[i][j-1]="C";
            s2.insert(s2.begin(),1,seq2[j-1]);
            s1.insert(s1.begin(),1,'-');
            j--;
        }
    }
    /*已確認computation是正確的*/
    for(int count = 0;count<=arc1.size();count++)
    {
        for(int count2 = 0;count2<=arc2.size();count2++)
        {
            cout<<M[count][count2]<<" ";
        }
        cout<<endl;
    }
    cout<<"跑之後"<<endl;
    for(int count = 0;count<=arc1.size();count++)
    {
        for(int count2 = 0;count2<=arc2.size();count2++)
        {
            cout<<D[count][count2]<<" ";
        }
        cout<<endl;
    }
    cout<<s1<<endl;
    cout<<s2<<endl;
    return temp_D[arc1.size()][arc2.size()];
}
string special_character_processing(string arc_string)
{
    string processed_string="";
    for(int i = 0;i<arc_string.size();i++)
    {
        if(arc_string[i]=='\"'||arc_string[i]=='\\'||arc_string[i]=='$')
        {
            processed_string.insert(processed_string.end(),1,'\\');
            processed_string.insert(processed_string.end(),1,arc_string[i]);
        }   
        else
        processed_string.insert(processed_string.end(),1,arc_string[i]);
    }
    return processed_string;
}

void write_error(double score,double seq_score ,const char *path)
{
    fstream file;
    file.open(path,ios::out);
    file<<"score is not correct!"<<endl;
    file<<astr1<<endl;
    file<<aseq1<<endl;
    file<<aseq2<<endl;
    file<<astr2<<endl;
    file<<"match score"<<match_cost<<endl;
    file<<"deletion score"<<deletion_cost<<endl;
    file<<"score"<<score<<endl;
    file<<"seq score="<<seq_score<<endl;
    file.close();
    
}
void write_data(double score,const char *path)
{
    /*處理特殊字元*/
    string processed_aseq1 = special_character_processing(aseq1);
    string processed_aseq2 = special_character_processing(aseq2);
    int Nmat= 0;
    for(int count = 0;count<aseq1.size();count++)
    {
        if(aseq1[count]!='-' && aseq2[count]!='-')
        ++Nmat;
    }
    fstream file;
    file.open(path,ios::out);
    file<<"<?"<<endl;
    file<<"$result_list=array();"<<endl;
    //file<<"$result_list[0][\"seq1\"] = \""<<astr1<<"\";"<<endl;
    //cout<<astr1<<endl;
    //file<<"$result_list[0][\"seq1\"] = \""<<aseq1<<"\";"<<endl;
    file<<"$result_list[0][\"seq1\"] = \""<<processed_aseq1<<"\";"<<endl;
    file<<"$result_list[0][\"arc1\"] = \""<<astr1<<"\";"<<endl;
    //file<<"$result_list[0][\"seq2\"] = \""<<aseq2<<"\";"<<endl;
    file<<"$result_list[0][\"seq2\"] = \""<<processed_aseq2<<"\";"<<endl;
    file<<"$result_list[0][\"arc2\"] = \""<<astr2<<"\";"<<endl;
    file<<"$result_list[0][\"start1\"] = \"0\";"<<endl;
    file<<"$result_list[0][\"start2\"] = \"0\";"<<endl;
    file<<"$result_list[0][\"end1\"] = \""<<seq1.size()-1<<"\";"<<endl;
    file<<"$result_list[0][\"end2\"] = \""<<seq2.size()-1<<"\";"<<endl;
    file<<"$result_list[0][\"Nmat\"] = \""<<Nmat<<"\";"<<endl;
    file<<"$result_list[0][\"score\"] = \""<<score<<"\";"<<endl;
    file<<"?>"<<endl;
    //file<<astr2<<endl;
    //file<<score<<endl;
    file.close();
}
void insert(alignment* start,int p1,int p2,double weight)
{
    //Construction of the alignment
    alignment* insertion = new alignment;
    insertion->p1     = p1;
    insertion->p2     = p2;
    insertion->weight = weight;
    insertion->next   = NULL; 
    alignment* iter   = start;
    if (p1==-1 && p2!=-1)
    while (iter->next!=NULL && p2>iter->next->p2)
    iter=iter->next;
    else if (p2==-1 && p1!=-1)
    while (iter->next!=NULL && p1>iter->next->p1)
    iter=iter->next;
    else
    while (iter->next!=NULL && p2>iter->next->p2 && p1>iter->next->p1 )
    iter=iter->next;
    insertion->next=iter->next;
    iter->next=insertion;
}
void traceback()
{
    vector<vector<string> >direct_lower;
    vector<vector<string> >direct_upper;
    vector<vector<double> >().swap(middle);
    vector<vector<double> >().swap(lower);
    vector<vector<double> >().swap(upper);
    vector<vector<double> >().swap(M);
    static double total=0.0;
    static double gap_insert = 0.0;
    static double gap_delete = 0.0;
    static double match_score = 0.0;
    double LARGE_NUMBER=999999;
    // stores aligned sequences and weights in aseq1,aseq2,astr1,astr2
    alignment* ali=new alignment;
    ali->p1=-1;  
    ali->p2=-1; 
    ali->next=NULL;
    stack<double> weight;
    double v1,v2,v3,v4;
    double gappenalty2; 
    double gappenalty1;
    double open_gappenalty2;
    double open_gappenalty1;
    double arc_cost = 0.0;
    // range is the currently computed sequence range
    four_tuple range;
    range.l1=0;
    range.l2=0;
    range.r1=seq1.size()-1;
    range.r2=seq2.size()-1;
    stack<four_tuple> ranges;
    ranges.push(range);
    vector<vector<string> >direct_array;
    while(!ranges.empty())
    {
        int l1=ranges.top().l1;
        int r1=ranges.top().r1;
        int l2=ranges.top().l2;
        int r2=ranges.top().r2;
          //cout<<"input r1="<<r1<<" l1="<<l1<<" r2="<<r2<<" l2="<<l2<<endl;
        ranges.pop();
        //幹你娘bug在這
        bool open_flag=true;
        if (l1>r1 && l2<=r2)
        {
            for (int s=r2;s>=l2;s--)
            {
                #ifdef _TRACE_DEBUG
                cout<<"seq2 start inserting gap"<<endl;
                #endif 
                if(open_flag==true)
                {
                    #ifdef _TRACE_DEBUG
                    cout<<"- "<<seq2[s]<<" open "<<(-common_exp-common_opp)<<endl;
                    #endif 
                    deletion_cost+=(-common_exp-common_opp);
                    insert(ali,-1,s,-common_exp-common_opp);
                    open_flag=false;
                }
                else
                {
                    #ifdef _TRACE_DEBUG
                    cout<<"- "<<seq2[s]<<" exp "<<(-common_exp)<<endl;
                    #endif 
                    insert(ali,-1,s,-common_exp);
                    deletion_cost+=(-common_exp);
                }
            }
            
        }
        else if (l1<=r1 && l2>r2)
        {
            for (int s=r1;s>=l1;s--)
            {
                #ifdef _TRACE_DEBUG
                cout<<"seq1 start inserting gap"<<endl;
                #endif 
                if(open_flag==true)
                {
                    #ifdef _TRACE_DEBUG
                    cout<<seq1[s]<<" - open "<<(-common_exp-common_opp)<<endl;
                    #endif 
                    insert(ali,s,-1,-common_exp-common_opp);
                    open_flag=false;
                    deletion_cost+=(-common_exp-common_opp);
                }
                else
                {
                    #ifdef _TRACE_DEBUG
                    cout<<seq1[s]<<" - exp "<<(-common_exp)<<endl;
                    #endif 
                    insert(ali,s,-1,-common_exp);
                    deletion_cost+=(-common_exp);
                }
            }
        }
        else if (l1<=r1 && l2<=r2)
        {
            // init and compute M
            M.resize(r1-l1+2);
            lower.resize(r1-l1+2);
            upper.resize(r1-l1+2);
            direct_array.resize(r1-l1+2);
            direct_lower.resize(r1-l1+2);
            direct_upper.resize(r1-l1+2);
            for(int s=0;s<M.size();s++)
            {
                M[s].resize(r2-l2+2);
                lower[s].resize(r2-l2+2);
                upper[s].resize(r2-l2+2);
                direct_array[s].resize(r2-l2+2);
                direct_lower[s].resize(r2-l2+2);
                direct_upper[s].resize(r2-l2+2);
            }
            /*---------affine gap penalty ------2016/5/20*/ 
            M[0][0]=lower[0][0]=upper[0][0]=0;
            direct_array[0][0]="o";
            direct_lower[0][0]="o";
            direct_upper[0][0]="o";
            for (int k=1;k<r1-l1+2;k++)
            {
                lower[k][0]=-LARGE_NUMBER; 
                //M[k][0]=upper[k][0]=-base_opp-(k*base_exp)+not_free1(l1+k-1)*(0.5*(-arc_opp-k*arc_exp)+(base_opp+(k*base_exp)));
                M[k][0]=upper[k][0]=-common_opp-(k*common_exp);
                direct_array[k][0]="B";
                direct_lower[k][0]="C";
                direct_upper[k][0]="B";
            }
            for (int l=1;l<r2-l2+2;l++)
            {
                upper[0][l]=-LARGE_NUMBER;
                //M[0][l]=lower[0][l]=-base_opp-(l*base_exp)+not_free2(l2+l-1)*(0.5*(-arc_opp-l*arc_exp)+(base_opp+(l*base_exp)));
                M[0][l]=lower[0][l]=-common_opp-(l*common_exp);
                direct_array[0][l]="C";
                direct_upper[0][l]="B";
                direct_lower[0][l]="C";
            }
            direct_lower[0][1]="A";
            direct_upper[1][0]="A";
            
            for (int k=1;k<r1-l1+2;k++)
                for (int l=1;l<r2-l2+2;l++)
                {
                    //max
                    v1=v2=v3=v4=-LARGE_NUMBER;
                    //min         v1=v2=v3=v4=LARGE_NUMBER;
                    int a1=l1+k-1;                      // a1,a2 sequence positions 
                    int a2=l2+l-1;
                    v1=lower[k][l]=max(lower[k][l-1]-common_exp,M[k][l-1]-common_exp-common_opp);
                    if(lower[k][l]==lower[k][l-1]-common_exp)
                        direct_lower[k][l]="C";
                    else if(lower[k][l]==M[k][l-1]-common_exp-common_opp)
                        direct_lower[k][l]="A";
                    
                    v2=upper[k][l]=max(upper[k-1][l]-common_exp,M[k-1][l]-common_exp-common_opp);
                    if(upper[k][l]==upper[k-1][l]-common_exp)
                        direct_upper[k][l]="B";
                    else if(upper[k][l]==M[k-1][l]-common_exp-common_opp)
                        direct_upper[k][l]="A";


                    //v1=M[k-1][l]+w_d+not_free1(a1)*(0.5*w_r-w_d);
                    //v2=M[k][l-1]+w_d+not_free2(a2)*(0.5*w_r-w_d);
                    v3=M[k-1][l-1]+base_matching(a1,a2)*w_m+(not_free1(a1)+not_free2(a2))*0.5*w_b;
                    if (arc1[a1]==')' && arc2[a2]==')') 
                    {
                        int i1=L1[I1[a1]];
                        int j1=L2[I2[a2]];
                        v4=M[i1-l1][j1-l2]+D[I1[a1]][I2[a2]]+arc_operation(i1,j1,a1,a2);
                        //(base_matching(i1,j1)+base_matching(a1,a2))*0.5*w_am;
                    }
                    //M[k][l]=min4(v1,v2,v3,v4);
                    M[k][l]=max4(v1,v2,v3,v4);

                    string temp_direct;
                    if(M[k][l]==v1)
                        temp_direct="C";
                    else if(M[k][l]==v2)
                        temp_direct="B";
                    else if(M[k][l]==v3)
                        temp_direct="A";
                    else if(M[k][l]==v4)
                        temp_direct="S";
                    else
                    {
                        //cout<<"error!"<<endl;
                        //cout<<"M["<<k<<"]["<<l<<"]="<<M[k][l]<<" v1="<<v1<<" v2="<<v2<<" v3="<<v3<<" v4="<<v4<<endl;
                        temp_direct="M";
                    }
                    direct_array[k][l]=temp_direct;
                }
            /* debug用*/
#ifdef _TRACE_DEBUG
            cout<<"score"<<endl;
            for(int count=0;count<=r1-l1+1;count++)
            {
                for(int count2=0;count2<=r2-l2+1;count2++)
                    cout<<M[count][count2]<<"\t";
                cout<<endl;
            }
            
            cout<<"lower"<<endl;
            for(int count=0;count<=r1-l1+1;count++)
            {
                for(int count2=0;count2<=r2-l2+1;count2++)
                    cout<<direct_lower[count][count2]<<"\t";
                cout<<endl;
            }
            cout<<"M"<<endl;
            for(int count=0;count<=r1-l1+1;count++)
            {
                for(int count2=0;count2<=r2-l2+1;count2++)
                    cout<<direct_array[count][count2]<<"\t";
                cout<<endl;
            }

            cout<<"upper"<<endl;
            for(int count=0;count<=r1-l1+1;count++)
            {
                for(int count2=0;count2<=r2-l2+1;count2++)
                    cout<<direct_upper[count][count2]<<"\t";
                cout<<endl;
            }
            /**/
#endif
            int k=r1-l1+1;
            int l=r2-l2+1;
            #ifdef _TRACE_DEBUG
            cout<<"origin k and l"<<endl;
            cout<<"k="<<k<<" l="<<l<<endl;
            #endif
          //cout<<"input r1="<<r1<<" l1="<<l1<<" r2="<<r2<<" l2="<<l2<<" k="<<k<<" l="<<l<<endl;
            bool seqaln=true;
            // sequence alignment
            while (seqaln)
            {

                int a1=l1+k-1;                      // a1,a2 sequence positions 
                int a2=l2+l-1;
                #ifdef _TRACE_DEBUG
                cout<<"a1="<<a1<<" a2="<<a2<<endl;
                cout<<"k="<<k<<" l="<<l<<endl;
                #endif
                if (k==0 && l==0)
                {
                    seqaln=false;
                }
                else if (k>0 && direct_array[k][l]=="B" )
                {
                    /*跳至direct_upper矩陣進行處理*/ 
                    bool search_flag = true;
                    while(search_flag)
                    {
                        a1=l1+k-1;
                        if(k!=0)
                        {
                            #ifdef _TRACE_DEBUG
                            cout<<"insert operation "<<k<<" "<<direct_upper[k][l]<<endl;
                            cout<<seq1[a1]<<" - index "<<a1<<endl;
                            #endif
                            if(direct_upper[k][l]=="B")
                            {
                                insert(ali,a1,-1,-common_exp);
                                deletion_cost+=-common_exp;
                                #ifdef _TRACE_DEBUG
                                cout<<"exp score="<<(-common_exp)<<endl;
                                #endif 
                            }
                            else if(direct_upper[k][l]=="A")
                            {
                                insert(ali,a1,-1,-common_exp-common_opp);
                                deletion_cost+=-common_exp-common_opp;
                                search_flag=false;
                                #ifdef _TRACE_DEBUG
                                cout<<"open score="<<(-common_exp-common_opp)<<endl;
                                #endif 
                            }

                        k--;
                        }
                        else
                        {
                            search_flag=false;
                        }
                    }
                    
                }
                //else if (l>0 && fabs(M[k][l]-(M[k][l-1]+w_d+not_free2(a2)*(0.5*w_r-w_d)))<eps )
                //else if (l>0 && (M[k][l]-lower[k][l])==0)
                else if (l>0 && direct_array[k][l]=="C")
                {
                    //cout<<k<<" "<<l<<" "<<direct_array[k][l]<<endl;
                    bool search_flag= true;
                    while(search_flag)
                    {
                        a2=l2+l-1;
                        
                        if(l!=0)
                        {
                            #ifdef _TRACE_DEBUG
                            cout<<"insert operation "<<l<<" "<<direct_lower[k][l]<<endl;
                            cout<<"- "<<seq2[a2]<<" index "<<a2<<endl;
                            #endif
                            if(direct_lower[k][l]=="C")
                            {
                                insert(ali,-1,a2,-common_exp);
                                #ifdef _TRACE_DEBUG
                                cout<<"exp score="<<(-common_exp)<<endl;
                                #endif 
                                deletion_cost+=-common_exp;
                            }
                            else if(direct_lower[k][l]=="A")
                            {
                                insert(ali,-1,a2,-common_exp-common_opp);
                                #ifdef _TRACE_DEBUG
                                cout<<"open score="<<(-common_exp-common_opp)<<endl;
                                #endif 
                                deletion_cost+=-common_exp-common_opp;
                                search_flag=false;
                            }

                        l--;
                        }
                        else
                            search_flag=false;
                    }
                }
                //else if (k>0 && l>0 && (M[k][l]-M[k-1][l-1]+base_matching(a1,a2)+(not_free1(a1)+not_free2(a2))*0.5*w_b)==0)
                else if (k>0 && l>0 && direct_array[k][l]=="A") 
                {
                    //cout<<k<<" "<<l<<" "<<direct_array[k][l]<<endl;
                    //cout<<seq1[a1]<<" "<<seq2[a2]<<base_matching(a1,a2)+(not_free1(a1)+not_free2(a2))*0.5*w_b<<endl;
                    match_score+=base_matching(a1,a2)+(not_free1(a1)+not_free2(a2))*0.5*w_b;
                    #ifdef _TRACE_DEBUG
                    cout<<"Base match"<<endl;
                    cout<<arc1[a1]<<" "<<arc2[a2]<<endl;
                    cout<<seq1[a1]<<" "<<seq2[a2]<<" index "<<a1<<" "<<a2<<endl;
                    cout<<"score="<<base_matching(a1,a2)+(not_free1(a1)+not_free2(a2))*0.5*w_b<<endl;

                    #endif
                    insert(ali,a1,a2,base_matching(a1,a2)+(not_free1(a1)+not_free2(a2))*0.5*w_b);
                    total+=base_matching(a1,a2)+(not_free1(a1)+not_free2(a2))*0.5*w_b;
                    k--;
                    l--;
                }
                else
                {
                    //cout<<"move to S"<<endl;
                    seqaln=false;
                }
            }
            int a1=l1+k-1;                      // a1,a2 sequence positions 
            int a2=l2+l-1;                      // right arc ends
            // base-pair alignment
            //cout<<"a1="<<a1<<" a2="<<a2<<" arc "<<arc1[a1]<<" "<<arc2[a2]<<endl;
            if (arc1[a1]==')' && arc2[a2]==')')
            {
                double w=M[L1[I1[a1]]-l1][L2[I2[a2]]-l2]+D[I1[a1]][I2[a2]]+arc_operation(L1[I1[a1]],L2[I2[a2]],a1,a2);//(base_matching(L1[I1[a1]],L2[I2[a2]])+base_matching(a1,a2))*0.5*w_am;
                if (fabs(M[k][l]-w)<eps)
                {

                    int i1=L1[I1[a1]];              // left arc ends
                    int j1=L2[I2[a2]];
                    double edge_weight=0.5*arc_operation(L1[I1[a1]],L2[I2[a2]],a1,a2); //(base_matching(L1[I1[a1]],L2[I2[a2]])+base_matching(a1,a2))*0.5*w_am;
                    arc_cost+=edge_weight*2;
                    #ifdef _TRACE_DEBUG
                    cout<<"confimed the path come from arc match"<<endl;
                    cout<<"arc match"<<endl;
                    cout<<arc1[i1]<<" "<<arc1[a1]<<endl;
                    cout<<seq1[i1]<<" "<<seq1[a1]<<" index "<<i1<<" "<<a1<<endl;
                    cout<<seq2[j1]<<" "<<seq2[a2]<<" index "<<j1<<" "<<a2<<endl;
                    cout<<arc2[j1]<<" "<<arc2[a2]<<endl; 
                    cout<<"score="<<(edge_weight*2)<<endl;
                    #endif
                    total+=edge_weight*2;
                    match_cost+=edge_weight*2;
                    insert(ali,i1,j1,edge_weight);
                    insert(ali,a1,a2,edge_weight);
                    four_tuple CR1,CR2;
                    CR1.l1=l1   ; CR1.r1=i1-1 ; CR1.l2=l2   ; CR1.r2=j1-1 ;
                    CR2.l1=i1+1 ; CR2.r1=a1-1 ; CR2.l2=j1+1 ; CR2.r2=a2-1 ;
                    #ifdef _TRACE_DEBUG
                    cout<<"之後處理 左邊CR1.r1="<<i1-1<<" CR1.l1="<<l1<<" CR1.r2="<<j1-1<<" CR1.l2="<<l2<<endl;
                    cout<<"優先處理 右邊CR2.r1="<<a1-1<<" CR2.l1="<<i1+1<<" CR2.r2="<<a2-1<<" CR2.l2="<<j1+1<<endl;
                    #endif
                    ranges.push(CR1);
                    ranges.push(CR2);
                }
            }
        }
    }
    cout<<"match score:"<<arc_cost<<"+"<<match_score<<endl;
    // write aligned sequences
    weights.resize(0);
    for(alignment* iter=ali->next;iter!=NULL;iter=iter->next)
    {
        aseq1.push_back((iter->p1==-1?'-':seq1[iter->p1]));
        astr1.push_back((iter->p1==-1?'-':arc1[iter->p1]));
        aseq2.push_back((iter->p2==-1?'-':seq2[iter->p2]));
        astr2.push_back((iter->p2==-1?'-':arc2[iter->p2]));
        weights.push_back(iter->weight);
    }
}
          /*------------------利用二元搜尋快速索引位置-----------------------*/
      int BinarySearch(char character,IndexMatrix *array,int length)
      {
          int low = 0;
          int high =length-1;
          while(low<=high)
          {
              int mid = (low+high)/2;
              if(array[mid].alphabet==character)
              {
                  //printf("%c %d",array[mid].alphabet,array[mid].pos);
                  return array[mid].pos;
              }
              else if(array[mid].alphabet>character)
              high = mid-1;
              else if(array[mid].alphabet<character)
              low = mid+1;
          }
          return -1;
      }
      int base_matching(int pos1,int pos2)
      {
          int index_x,index_y;
          //time complexity = O(logN) * 2
          index_x = BinarySearch(seq1[pos1],alphabet_index,size);
          index_y = BinarySearch(seq2[pos2],alphabet_index,size);
          //cout<<seq1[pos1]<<" vs "<<seq2[pos2]<<" "<<scoring_matrix[index_x][index_y]<<endl;
          /*
          if(scoring_matrix[index_x][index_y]>=0)
          {
          return 0;
      }
          else
          return 1;
          */
          return scoring_matrix[index_x][index_y];
      }
      /*讀取序列資訊*/
      void read_data(const char *path)
      {
          //FILE* file = fopen(path,"r");
          fstream file;
          file.open(path,ios::in);
          if(!file)
          {
              printf("Can't find your input file %s\n",path);
              exit(-1);
          }
          /*暫存buffer*/
          string sub;
          int temp_size;
          int option=0;
          /*read <?php*/
          getline(file,sub,'\n');
          /*read seq1*/    
          getline(file,sub,'\n');
          temp_size = sub.size();
          seq1= sub.substr(7,temp_size-9);
          //cout<<sub<<" length="<<sub.size()<<" seq length="<<sub.size()-9<<endl;
          /*read arc1*/    
          getline(file,sub,'\n');
          temp_size = sub.size();
          arc1 = sub.substr(7,temp_size-9);
          //cout<<sub<<" length="<<sub.size()<<" arc length="<<sub.size()-9<<endl;
          /*read seq2*/    
          getline(file,sub,'\n');
          temp_size = sub.size();
          seq2 = sub.substr(7,temp_size-9);
          //cout<<sub<<" length="<<sub.size()<<" seq length="<<sub.size()-9<<endl;
          /*read arc2*/    
          getline(file,sub,'\n');
          temp_size = sub.size();
          arc2 = sub.substr(7,temp_size-9);
          //cout<<sub<<" length="<<sub.size()<<" arc length="<<sub.size()-9<<endl;
          /*read matrixpath*/    
          getline(file,sub,'\n');
          temp_size = sub.size();
          matrixpath = sub.substr(10,temp_size-12);
          //cout<<sub<<endl;
          /*read gap_opp*/    
          getline(file,sub,'\n');
          temp_size = sub.size();
          sub = sub.substr(6,temp_size-8);
          gap_opp = atoi(sub.c_str());
          //cout<<sub<<endl;
          /*read gap_exp*/    
          getline(file,sub,'\n');
          temp_size = sub.size();
          sub = sub.substr(6,temp_size-8);
          gap_exp = atoi(sub.c_str());
          //cout<<sub<<endl;
          /*Prevent the errors happend! If the sequence's length is not equal to the arc's length */
          int length_seq1 = seq1.size(),length_arc1=arc1.size(),length_seq2=seq2.size(),length_arc2=arc2.size();
          if(length_seq1!=length_arc1 || length_seq2!=length_arc2)
          {
              cout<<seq1<<"\n"<<arc1<<"\n"<<seq2<<"\n"<<arc2<<"\n"<<gap_opp<<"\n"<<gap_exp<<endl;
              cout<<"Fatel error!  the size of the sequence is not equal to the size of the arc"<<endl;
              cout<<"Maybe is the seq1 or seq2..."<<endl;
              cout<<seq1.size()<<endl;
              cout<<arc1.size()<<endl;
              cout<<seq2.size()<<endl;
              cout<<arc2.size()<<endl;
              exit(1);
          }
      }
      /*讀取矩陣資訊*/
      int** readmat(char *file_name)
      {
          /*開啟Scoring Matrix*/   
          FILE *fptr = fopen(file_name,"r");
          if(fptr==NULL)
          {
              printf("File open failed!\n");
              exit(EXIT_FAILURE);
          }
          else
          {
              printf("Successful open file!\n");
              puts(file_name);
          }
          /*預設每一段長度最大為400*/
          char line[400];
          /*first reading data is comment, so lets ingnore it!*/
          fgets(line,400,fptr);
          /*Start reading data from scoring matrix*/
          fgets(line,400,fptr);
          /*配置暫存矩陣進行內容的讀取*/
          char* row_array = (char*)malloc(400*sizeof(char));
          char* column_array = (char*)malloc(400*sizeof(char));
          /*---------------counting------------------------*/
          /*
          *
          *First at all, I used test_index to caculate how many data that I should store into array. And also, I used very unsmartly way to store alphabet characters. 
          *Comment Time : 2015/7/10
          *
          *
          */
          /*吃掉空白部分只留下字元*/ 
          char *temp = strtok(line," ");
          char str[92]="";
          row_array[0] = *temp;
          column_array[0]=*temp;
          int test_index=1;
          int loop = 1;
          do
          {
              strcat(str,temp); 
              temp = strtok(NULL," ");
              if(temp!=NULL)
              {
                  column_array[loop]=*temp;
                  row_array[loop++]=*temp;
                  test_index++;
              }
          }while(temp!=NULL);
          free(row_array);
          free(column_array);
          /*----------------------------------------------*/
      /*initialize*/ 
      size = test_index;
      int **array = (int**)malloc(size*sizeof(void *));
      int count;
      alphabet_index = (IndexMatrix*)malloc(size*sizeof(IndexMatrix));
      for( count = 0;count<size;count++  )
      {
          array[count] = (int*)malloc(size*sizeof(int));
          alphabet_index[count].alphabet=str[count];
          alphabet_index[count].pos=count;
      }
      /*Array index value*/
      int array_row ,array_col = -1;
      loop =0;
      /*讀取每一個字並存入*/
      while(!feof(fptr))
      {
          array_row = 0; 
          char *temp = strtok(line," ");
          //printf("%s ",temp);
          temp = strtok(NULL," ");
          while(temp!=NULL)
          {
              array[array_col][array_row] = atoi(temp);
              //debug function
              //printf("%d ",atoi(temp));
              array_row++;
              temp = strtok(NULL," ");
          }
          array_col++;
          fgets(line,400,fptr);
          //debug function
          //printf("\n");
      }
      fclose(fptr);
      /*----------------Testing-------------------*/
      return array;
    }
/*計算best alignment*/
double computation()
{
    vector<vector<string> > B;
    vector<vector<string> > C;

    vector<vector<string> > direct_array;
    double LARGE_NUMBER=999999;
    /*對I1 I2兩個vector配置arc1 and arc2大小的記憶體空間*/
    I1.resize(arc1.size());
    I2.resize(arc2.size());
    /*依照由內而外的方式取得arc1*/
    int index = 0;
    for(int i = 0;i<arc1.size();i++)
    {
        I1[i]=-1;
        if(arc1[i]=='(')
           str_stack.push(i);
           else if(arc1[i]==')')
        {
            int last = str_stack.top();
            str_stack.pop();
            I1[i]=I1[last]=index++;
            L1.push_back(last);
            R1.push_back(i);
        }
    }
    //cout<<"total base pair = "<<L1.size()<<endl;
    /*依照由內而外的方式取得arc2*/
    index=0;
    for(int i=0;i<arc2.size();i++)
    {
        I2[i]=-1;
        if(arc2[i]=='(')
           str_stack.push(i);
           else if(arc2[i]==')')
        {
            int last = str_stack.top();
            str_stack.pop();
            I2[i]=I2[last]=index++;
            L2.push_back(last);
            R2.push_back(i);
        }
    }
    #ifdef _TRACE_DEBUG
    cout<<"base pair 1"<<endl;
    for(int i = 0;i<L1.size();i++)
    cout<<L1[i]<<" "<<R1[i]<<" "<<I1[i]<<endl;
    cout<<"base pair 2"<<endl;
    for(int i = 0;i<L2.size();i++)
    cout<<L2[i]<<" "<<R2[i]<<" "<<I2[i]<<endl;
    #endif
    /*initialize*/
    D.resize(L1.size());
    //cout<<"D matrix szie = "<<D.size()<<endl;
    for(int i = 0;i<D.size();i++)
    {
        D[i].resize(L2.size());
    }
/*First processing arc align arc */
for(int i = 0;i<L1.size();i++)
{
    for(int j = 0;j<L2.size();j++)
    {
        //added affine gap penaly matrix : lower, upper, middle
        lower.resize(R1[i]-L1[i]);
        upper.resize(R1[i]-L1[i]);
        M.resize(R1[i]-L1[i]);

        #ifdef _TRACE_DEBUG
        B.resize(R1[i]-L1[i]);
        C.resize(R1[i]-L1[i]);
        direct_array.resize(R1[i]-L1[i]);
        #endif


        for(int s = 0;s<M.size();s++)
        {
            M[s].resize(R2[j]-L2[j]);
            //affine gap penalty 2016/5/20
            lower[s].resize(R2[j]-L2[j]);
            upper[s].resize(R2[j]-L2[j]);
            #ifdef _TRACE_DEBUG
            B[s].resize(R2[j]-L2[j]);
            C[s].resize(R2[j]-L2[j]);
            direct_array[s].resize(R2[j]-L2[j]);
            #endif
        }
        M[0][0]=0;
        //affine gap penalty 2016/5/20
        lower[0][0]=upper[0][0]=0;
        #ifdef _TRACE_DEBUG
        B[0][0]=C[0][0]=direct_array[0][0]="o";
        #endif

        for(int k = 1;k<R1[i]-L1[i];k++)
        {
            lower[k][0]=-LARGE_NUMBER; 
            M[k][0]=upper[k][0]=-common_opp-k*common_exp;
            //M[k][0]=M[k-1][0]+w_d+not_free1(L1[i]+k)*(0.5*w_r-w_d);
            #ifdef _TRACE_DEBUG
            B[k][0]="B";
            C[k][0]="C";
            direct_array[k][0]="B";
            #endif
        }
        for (int l=1;l<R2[j]-L2[j];l++)
        {
            upper[0][l]=-LARGE_NUMBER;
            M[0][l]=lower[0][l]=-common_opp-l*common_exp;
            //M[0][l]=M[0][l-1]+w_d+not_free2(L2[j]+l)*(0.5*w_r-w_d);
            #ifdef _TRACE_DEBUG
            B[0][l]="B";
            C[0][l]="C";
            direct_array[0][l]="C";
            #endif
        }
        //compute M
        double v1,v2,v3,v4;
        for(int k=1;k<R1[i]-L1[i];k++)
        for(int l=1;l<R2[j]-L2[j];l++)
        {
            //max 
            v1 = v2 = v3 = v4 = -LARGE_NUMBER;
            //min                 v1 = v2 = v3 = v4 = LARGE_NUMBER;
            int a1 = L1[i]+k;
            int a2 = L2[j]+l;
            //edit operation
            //v1 = lower gap penalty , v2 = upper gap penalty, v3 = mismatch or match gap penalty, v4=special case
            v1=lower[k][l]=max(lower[k][l-1]-common_exp,M[k][l-1]-common_opp-common_exp);
            #ifdef _TRACE_DEBUG
            if(lower[k][l]==lower[k][l-1]-common_exp)
                C[k][l]="C";
            else if(lower[k][l]==M[k][l-1]-common_exp-common_opp)
                C[k][l]="A";
            #endif
            v2=upper[k][l]=max(upper[k-1][l]-common_exp,M[k-1][l]-common_opp-common_exp);
            #ifdef _TRACE_DEBUG
            if(upper[k][l]==upper[k-1][l]-common_exp)
                B[k][l]="B";
            else if(upper[k][l]==M[k-1][l]-common_opp-common_exp)
                B[k][l]="A";
            #endif
            //v3=middle[k][l]=max(middle[k-1][l-1]+base_matching(a1,a2)+((not_free1(a1)+not_free2(a2))*0.5*w_b),lower[k-1][l-1],upper[k-1][l-1]);
            //v1=M[k-1][l]+w_d+not_free1(a1)*(0.5*w_r-w_d);
            //v2=M[k][l-1]+w_d+not_free2(a2)*(0.5*w_r-w_d);
            v3=M[k-1][l-1]+base_matching(a1,a2)*w_m+(not_free1(a1)+not_free2(a2))*0.5*w_b;
            if(arc1[a1]==')' && arc2[a2]==')')
            {
                int leftpoint = L1[I1[a1]];
                int leftpoint2 = L2[I2[a2]];
                v4 = M[leftpoint - L1[i]-1][leftpoint2 - L2[j]-1]+D[I1[a1]][I2[a2]]+arc_operation(L1[I1[a1]],L2[I2[a2]],R1[I1[a1]],R2[I2[a2]]);
                //(base_matching(L1[I1[a1]],L2[I2[a2]])+base_matching(R1[I1[a1]],R2[I2[a2]]))*0.5*w_am;
            }
            //M[k][l]=min4(v1,v2,v3,v4);
            //
            M[k][l]=max4(v1,v2,v3,v4);
            #ifdef _TRACE_DEBUG
            if(M[k][l]==v1)
                direct_array[k][l]="C";
            else if(M[k][l]==v2)
                direct_array[k][l]="B";
            else if(M[k][l]==v3)
                direct_array[k][l]="A";
            else if(M[k][l]==v4)
                direct_array[k][l]="S";
            #endif
        }
        D[i][j]=M[R1[i]-L1[i]-1][R2[j]-L2[j]-1];
        #ifdef _TRACE_DEBUG
        cout<<"base pair vs base pair"<<endl;
        cout<<"arc1: "<<R1[i]<<" "<<L1[i]<<endl;
        cout<<"arc2: "<<R2[j]<<" "<<L2[j]<<endl;
        cout<<"score"<<endl;
        for(int x = 0;x<R1[i]-L1[i];x++)
        {
            for(int y=0;y<R2[j]-L2[j];y++)
                cout<<M[x][y]<<"\t";
            cout<<endl;
        }
        cout<<"lower"<<endl;
        for(int x = 0;x<R1[i]-L1[i];x++)
        {
            for(int y=0;y<R2[j]-L2[j];y++)
                cout<<C[x][y]<<"\t";
            cout<<endl;
        }
        cout<<"M"<<endl;
        for(int x = 0;x<R1[i]-L1[i];x++)
        {
            for(int y=0;y<R2[j]-L2[j];y++)
                cout<<direct_array[x][y]<<"\t";
            cout<<endl;
        }
        cout<<"upper"<<endl;
        for(int x = 0;x<R1[i]-L1[i];x++)
        {
            for(int y=0;y<R2[j]-L2[j];y++)
                cout<<B[x][y]<<"\t";
            cout<<endl;
        }
        cout<<endl;
        #endif
    }
}
/*affine gap penalty 2016/5/20*/
lower.resize(arc1.size()+1);
upper.resize(arc1.size()+1);
M.resize(arc1.size()+1);
direct_array.resize(arc1.size()+1);
for(int i =0; i<M.size();i++)
{
    M[i].resize(arc2.size()+1);
    lower[i].resize(arc2.size()+1);
    upper[i].resize(arc2.size()+1);
    direct_array[i].resize(arc2.size()+1);
}
M[0][0]=lower[0][0]=upper[0][0]=0;
direct_array[0][0]="o";
/*initalize row*/
for(int k=1;k<=arc1.size();k++)
{
    lower[k][0]=-LARGE_NUMBER; 
    M[k][0]=upper[k][0]=-common_opp-k*common_exp;
    direct_array[k][0]="B";
    //M[k][0]=M[k-1][0]+w_d+not_free1(k-1)*(0.5*w_r-w_d);
}
/*initalize column*/
for(int l = 1;l<=arc2.size();l++)
{
    upper[0][l]=-LARGE_NUMBER;
    M[0][l]=lower[0][l]=-common_opp-l*common_exp;
    direct_array[0][l]="C";
    //M[0][l]=M[0][l-1]+w_d+not_free2(l-1)*(0.5*w_r-w_d);
}
//compute M
double v1,v2,v3,v4;
for (int k=1;k<=arc1.size();k++)
for (int l=1;l<=arc2.size();l++)
{
    //min v1 = v2 = v3 = v4 = LARGE_NUMBER;
    //max
    v1 = v2 = v3 = v4 = -LARGE_NUMBER;
    //origin version
    //v1=M[k-1][l]+w_d+not_free1(k-1)*(0.5*w_r-w_d);
    //v2=M[k][l-1]+w_d+not_free2(l-1)*(0.5*w_r-w_d);
    //v3=M[k-1][l-1]+base_matching(k-1,l-1)*w_m+(not_free1(k-1)+not_free2(l-1))*0.5*w_b;
    //affine gap penalty version 2016/5/20 
    v1=lower[k][l]=max(lower[k][l-1]-common_exp,M[k][l-1]-common_exp-common_opp);
    v2=upper[k][l]=max(upper[k-1][l]-common_exp,M[k-1][l]-common_exp-common_opp);
    //v1=lower[k][l]=max(lower[k][l-1]-common_exp,middle[k][l-1]-common_opp-common_exp,upper[k][l-1]-common_opp-common_exp);
    //v2=upper[k][l]=max(upper[k-1][l]-common_exp,middle[k-1][l]-common_opp-common_exp,lower[k-1][l]-common_opp-common_exp);
    //v3=middle[k][l]=max(middle[k-1][l-1]+base_matching(k-1,l-1)+((not_free1(k-1)+not_free2(l-1))*0.5*w_b),lower[k-1][l-1],upper[k-1][l-1]);
    v3=M[k-1][l-1]+base_matching(k-1,l-1)*w_m+(not_free1(k-1)+not_free2(l-1))*0.5*w_b;
    if(arc1[k-1]==')'&& arc2[l-1]==')')
    {
        int leftpoint = L1[I1[k-1]];
        int leftpoint2 = L2[I2[l-1]];
        //v4=middle[leftpoint][leftpoint2]+D[I1[k-1]][I2[l-1]]+arc_operation(L1[I1[k-1]],L2[I2[l-1]],R1[I1[k-1]],R2[I2[l-1]]);
        v4=M[leftpoint][leftpoint2]+D[I1[k-1]][I2[l-1]]+arc_operation(L1[I1[k-1]],L2[I2[l-1]],R1[I1[k-1]],R2[I2[l-1]]);
        //(base_matching(L1[I1[k-1]],L2[I2[l-1]])+base_matching(R1[I1[k-1]],R2[I2[l-1]]))*0.5*w_am;
    }
    //M[k][l]=min4(v1,v2,v3,v4);
    M[k][l]=max4(v1,v2,v3,v4);
    if(M[k][l]==v1)
        direct_array[k][l]="C";
    else if(M[k][l]==v2)
        direct_array[k][l]="B";
    else if(M[k][l]==v3)
        direct_array[k][l]="A";
    else if(M[k][l]==v4)
        direct_array[k][l]="S";
    //M[k][l]=max4(lower[k][l],upper[k][l],middle[k][l],v4);
}
#ifdef _TRACE_DEBUG
cout<<"Computation model results"<<endl;
for(int count = 0;count<=arc1.size();count++)
{
    for(int count2= 0;count2<=arc2.size();count2++)
        cout<<M[count][count2]<<"\t";
    cout<<endl;
}
cout<<endl;
for(int count=0;count<=arc1.size();count++)
{
    for(int count2=0;count2<=arc2.size();count2++)
        cout<<direct_array[count][count2]<<"\t";
    cout<<endl;
}
cout<<"\n"<<endl;
#endif

return M[arc1.size()][arc2.size()];
}
