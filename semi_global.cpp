/*************************************************************************
	> File Name: semi_global.cpp
	> Author: 
	> Mail: 
	> Created Time: Thu Mar 24 11:35:41 2016
    > Description: semi序列比對且未考慮affine gap penalty
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

#define _TRACE_DEBUG
#define _DISPLAY
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
static double max_row_value=0.0,max_column_value=0.0;
static int current_column,current_row;
static double match_score =0.0;
static double arc_score = 0.0;
static double del_score = 0.0;
vector<double> weights;
vector<int>    L1,R1,I1,L2,R2,I2;
stack<int>     str_stack;
vector<vector<double> > M;
vector<vector<double> > lower;
vector<vector<double> > upper;

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


double arc_match_case1 = 0.0; //if both > 0
double arc_match_case2 = 0.0; // if both =0
double arc_mismatch_case1=0.0;// if sum >0
double arc_mismatch_case2=0.0;// if sum = 0
double arc_mismatch_case3=0.0;// if sum < 0
double arc_mismatch_case4=0.0;// if both < 0
double w_d =0;  // base deletion
double w_r =0;  // arc  removing
double w_b =0;  // arc  breaking

double number=1;



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
        if((base_matching(p1,p2)+base_matching(p3,p4))==0)
            return arc_mismatch_case2;
        else
            return arc_match_case1*(base_matching(p1,p2)+base_matching(p3,p4));
    }
    //當alignment score <0 views as arc-mismatch
    else
    {
        if((base_matching(p1,p2)<0 && base_matching(p3,p4)>=0)||(base_matching(p1,p2)>=0&&base_matching(p3,p4)<0))
        {
            if((base_matching(p1,p2)+base_matching(p3,p4))>0)
            {

                return (base_matching(p1,p2)+base_matching(p3,p4))*arc_mismatch_case1;
            } 
            else if((base_matching(p1,p2)+base_matching(p3,p4))==0)
                return arc_mismatch_case2;
            else
            {
                return (base_matching(p1,p2)+base_matching(p3,p4))*(arc_mismatch_case3);
            }
        }
        else if(base_matching(p1,p2)<=0 && base_matching(p3,p4)<0)
        //return (base_matching(p1,p2)+base_matching(p3,p4))*0.5*w_am;
            return (base_matching(p1,p2)+base_matching(p3,p4))*arc_mismatch_case4;
    }
}



int BinarySearch(char ,IndexMatrix *,int );
void traceback();
void insert(alignment*,int,int,double);
int** readmat(char *);
void read_data(const char *);
void write_data(double,const char *);
void write_profit(const char *);
double computation();
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
    const char *profit_result;
    if(argc!=11)
    {
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/test_data/1A9N_Q_to_1E7K_C/semi_input.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1M90_B_to_1NKW_9/semi_input.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1IBK_A_to_1N33_A/semi_input.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1BAU_A_to_1S9S_A/semi_input.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1FQZ_A_to_1KP7_A/semi_input.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1F84_A_to_1P5N_A/semi_input.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1BYJ_A_to_1I6U_C/semi_input.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1FHK_A_to_1S9S_A/semi_input.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1JJ2_9_to_1NKW_9/semi_input.php";
        pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1UN6_F_to_1JUR_A/semi_input.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1E8O_E_to_1JID_B/semi_input.php";
       //pdb_compare_path= "/home/watchlee/Research_Programming/AlignmentToHen/Test_for_GLOBAL/1NJI_B_to_1NMY_9.php";
       //pdb_compare_path= "/home/watchlee/Research_Programming/AlignmentToHen/Special_for_debug/semi_input.php";
    //   pdb_compare_path= "/home/watchlee/Research_Programming/AlignmentToHen/Special_for_debug/1E8O_Evs1JID_B/semi_input.php";
       //pdb_compare_path= "/home/watchlee/Research_Programming/AlignmentToHen/Special_for_debug/input2.php";
        pdb_result= "/home/watchlee/result.php";
        profit_result= "/home/watchlee/profit_result";
        number = 1;
         arc_match_case1 = 1; //if both > 0
         arc_mismatch_case1=1;// if sum >0
         arc_mismatch_case2=0;// if sum = 0
         arc_mismatch_case3=1;// if sum < 0
         arc_mismatch_case4=1;// if both < 0
         w_d =-100;  // base deletion
         w_r =-1;  // arc  removing
         w_b =-1;  // arc  breaking
    }
    else
    {
       pdb_compare_path=argv[1] ;
       pdb_result= argv[2];
         arc_match_case1 = atof(argv[3]); //if both > 0
         arc_mismatch_case1=atof(argv[4]);// if sum >0
         arc_mismatch_case2=atof(argv[5]);// if sum = 0
         arc_mismatch_case3=atof(argv[6]);// if sum < 0
         arc_mismatch_case4=atof(argv[7]);// if both < 0
         w_d =atof(argv[8]);  // base deletion
         w_r =atof(argv[9]);  // arc  removing
         w_b =atof(argv[10]);  // arc  breaking
    }

    char path[100];
    //-------------Test Scoring Matrix
    //sprintf(path,"./SM/BLOSUM-like_scoring_matrix");
    //sprintf(path,"./SM/iPARTS2_new_23C_4L_matrix");
    sprintf(path,"./SM/23-4L_matrix");
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
    double score = computation();
    traceback();
    #ifdef _TRACE_DEBUG
    cout<<"Match Score="<<match_score<<endl;
    cout<<"Arc Score="<<arc_score<<endl;
    cout<<"del Score="<<del_score<<endl;
    cout<<"Score="<<score<<endl;
    #endif 
    double total = 0.0;
    for(int count = 0;count<weights.size();count++)
    {
        cout<<weights[count]<<"\t";

        total+=weights[count];
    }
    cout<<endl;
    #ifdef _DISPLAY
    cout<<astr1<<endl;
    cout<<aseq1<<" "<<aseq1.size()<<endl;
    cout<<aseq2<<" "<<aseq2.size()<<endl;
    cout<<astr2<<endl;
    #endif
    #ifdef _DISPLAY
    cout<<"result = "<<score<<" "<<total<<endl;
    #endif
    if(abs(total-score)>0.01)
    {
        cout<<"incorrect"<<endl;
        cout<<pdb_compare_path<<endl;
    }
    else
        cout<<"correct"<<endl;

    write_data(score,pdb_result);
    
    write_profit(profit_result);

    /*釋放記憶體*/
    free(alphabet_index);//來源 line:33
    return 0;
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
void write_profit(const char *path)
{ 
    /*處理特殊字元*/

    fstream file;
    file.open(path,ios::out);
    file<<aseq1<<endl;
    file<<aseq2<<endl;


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
  int times=0;
  // stores aligned sequences and weights in aseq1,aseq2,astr1,astr2
  alignment* ali=new alignment;
  ali->p1=-1;  
  ali->p2=-1; 
  ali->next=NULL;
  stack<double> weight;
  double v1,v2,v3,v4;
  vector<vector <string> > direct;
  // range is the currently computed sequence range
  four_tuple range;
  range.l1=0;
  range.l2=0;
  range.r1=seq1.size()-1;
  range.r2=seq2.size()-1;
  /*做第一次*/
  stack<four_tuple> ranges;
  ranges.push(range);
  bool arc_flag= false; 
  bool first_time = true;
  while(!ranges.empty())
    {
      int l1=ranges.top().l1;
      int r1=ranges.top().r1;
      int l2=ranges.top().l2;
      int r2=ranges.top().r2;
      ranges.pop();
        /*seq2插入gap*/
      #ifdef _TRACE_DEBUG
      cout<<"current r1="<<r1<<" l1="<<l1<<" r2="<<r2<<" l2="<<l2<<endl;
      #endif
      if (l1>r1 && l2<=r2)
        {

              #ifdef _TRACE_DEBUG
            cout<<"seq2 insert gap r1="<<r1<<" l1="<<l1<<" r2="<<r2<<" l2="<<l2<<endl;
              #endif
            for (int s=r2;s>=l2;s--)
            {
                if(l1!=0)
                {
                    insert(ali,-1,s,w_d+not_free2(s)*(0.5*w_r-w_d));
                    del_score+=w_d;
                }
                else
                {
                    insert(ali,-1,s,0);

                }
            }
        }
      else if (l1<=r1 && l2>r2)
        {
              #ifdef _TRACE_DEBUG
          cout<<"seq1 insert gap r1="<<r1<<" l1="<<l1<<" r2="<<r2<<" l2="<<l2<<endl;
              #endif
            /*seq1插入gap*/
            for (int s=r1;s>=l1;s--)
            {
                if(l2!=0)
                {
                    
                     insert(ali,s,-1,w_d+not_free1(s)*(0.5*w_r-w_d));
                    del_score+=w_d;
                }
                else
                     insert(ali,s,-1,0);
                    
            }

        }
        /*做正常的alignment*/
      else if (l1<=r1 && l2<=r2)
	{
	  // init and compute M
	  M.resize(r1-l1+2);
      direct.resize(r1-l1+2);
	  for(int s=0;s<M.size();s++)
	    {
	      M[s].resize(r2-l2+2);
	      direct[s].resize(r2-l2+2);
	      for(int t=0;t<M[s].size();t++)
            M[s][t]=0;
	    }
	#ifdef _TRACE_DEBUG
        cout<<"ARC_FLAG="<<arc_flag<<endl;
    #endif
	  M[0][0]=0;
	  direct[0][0]="o";
	  for (int k=1;k<r1-l1+2;k++)
        {
            if(!arc_flag&&l1==0)
                M[k][0]=0; 
            else
                M[k][0]=M[k-1][0]+w_d+not_free1(l1+k-1)*(0.5*w_r-w_d);
        direct[k][0]="B";
        }
	  
	  for (int l=1;l<r2-l2+2;l++)
        {
            if(!arc_flag&&l2==0)
                M[0][l]=0;
            else
                M[0][l]=M[0][l-1]+w_d+not_free2(l2+l-1)*(0.5*w_r-w_d);
            direct[0][l]="C";
        } 
      if(arc_flag)
        arc_flag=false;

	  for (int k=1;k<r1-l1+2;k++)
	    for (int l=1;l<r2-l2+2;l++)
	      {
		//max
        v1=v2=v3=v4=-10000;
        //min 		v1=v2=v3=v4=10000;
		int a1=l1+k-1;                      // a1,a2 sequence positions 
		int a2=l2+l-1;
		
		v1=M[k-1][l]+w_d+not_free1(a1)*(0.5*w_r-w_d);
		v2=M[k][l-1]+w_d+not_free2(a2)*(0.5*w_r-w_d);
		v3=M[k-1][l-1]+base_matching(a1,a2)*w_m+(not_free1(a1)+not_free2(a2))*w_b;
		
		if (arc1[a1]==')' && arc2[a2]==')') 
		  {
		    int i1=L1[I1[a1]];
		    int j1=L2[I2[a2]];
		    v4=M[i1-l1][j1-l2]+D[I1[a1]][I2[a2]]+arc_operation(i1,j1,a1,a2);
		      //(base_matching(i1,j1)+base_matching(a1,a2))*0.5*w_am;
		  }
		//M[k][l]=min4(v1,v2,v3,v4);
		M[k][l]=max4(v1,v2,v3,v4);
             if(M[k][l]==v1)
                direct[k][l]="B";
             else if(M[k][l]==v2)
                direct[k][l]="C";
             else if(M[k][l]==v3)
                direct[k][l]="A";
             else if(M[k][l]==v4)
                direct[k][l]="S";
	      }
	#ifdef _TRACE_DEBUG
    cout<<"\t";
    for(int count =0;count<r2-l2+2;count++)
    {
        cout<<"a2:"<<count<<"\t";
    }
    cout<<"\t";
    cout<<endl;
    for(int count =0; count<r1-l1+2;count++)
    {
        cout<<"a1:"<<count<<"\t";
        for(int count2=0;count2<r2-l2+2;count2++)
            cout<<M[count][count2]<<"\t";
        cout<<endl;
    }       
    cout<<"Direct array"<<endl;
    cout<<"\t";
    for(int count =0;count<r2-l2+2;count++)
    {
        cout<<"a2:"<<count<<"\t";
    }
    cout<<"\t";
    cout<<endl;
    for(int count =0; count<r1-l1+2;count++)
    {
        cout<<"a1:"<<count<<"\t";
        for(int count2=0;count2<r2-l2+2;count2++)
            cout<<direct[count][count2]<<"\t";
        cout<<endl;
    }
    #endif
	  bool seqaln=true;
      /*由外而內*/
	  int k=r1-l1+1;
	  int l=r2-l2+1;
#ifdef _TRACE_DEBUG
cout<<"Origin r1 and l1"<<endl;
cout<<"r1="<<r1<<" l1="<<l1<<endl;
cout<<"Origin r2 and l2"<<endl;
cout<<"r2="<<r2<<" l2="<<l2<<endl;
cout<<"Origin k and l"<<endl;
cout<<"k="<<k<<" l="<<l<<endl;
cout<<"Last column and Last row"<<endl;
cout<<max_column_value<<" "<<max_row_value<<endl;
cout<<current_column<<" "<<current_row<<endl;
#endif

      if(first_time)
      {
          if(max_column_value>max_row_value)
            k=l1+current_column;
          else
              l=l2+current_row;
        for(int count = r1;count>k-1;count--)
          {
              cout<<"insert!!"<<count<<endl;
            insert(ali,count,-1,0);
          }
        for(int count = r2;count>l-1;count--)
          {
              cout<<"deletion!!"<<count<<endl;
            insert(ali,-1,count,0);
          }

      }
      first_time=false;
        while(seqaln)
        {
            int a1=l1+k-1;
           int a2=l2+l-1;
            #ifdef _TRACE_DEBUG
                cout<<"a1="<<a1<<"\ta2="<<a2<<endl;
            #endif
            if(k==0&&l==0)
                seqaln=false;
            else if(direct[k][l]=="B")
            {
                #ifdef _TRACE_DEBUG
                cout<<"發生insert"<<endl;
                cout<<arc1[a1]<<endl;
                cout<<seq1[a1]<<endl;
                cout<<"-"<<endl;
                cout<<"-"<<endl;
                cout<<"l = "<<l<<" l2="<<l2<<" insert cost= "<<w_d+not_free1(a1)*(0.5*w_r-w_d)<<endl;
                #endif
                if(l!=0||l2!=0)
                  insert(ali,a1,-1,w_d+not_free1(a1)*(0.5*w_r-w_d));
                else
                  insert(ali,a1,-1,0);
                k--;
            }
             else if(direct[k][l]=="C")
            {
                #ifdef _TRACE_DEBUG
                cout<<"發生deletion"<<endl;
                cout<<"-"<<endl;
                cout<<"-"<<endl;
                cout<<seq2[a2]<<endl;
                cout<<arc2[a2]<<endl;
                cout<<"k = "<<k<<" l1="<<l1<<" delete cost= "<<w_d+not_free1(a1)*(0.5*w_r-w_d)<<endl;
                #endif
               
                if(k!=0||l1!=0)
                      insert(ali,-1,a2,w_d+not_free2(a2)*(0.5*w_r-w_d));
                else
                      insert(ali,-1,a2,0);
                l--;
            }
            else if(direct[k][l]=="A")
            {
              insert(ali,a1,a2,base_matching(a1,a2)*w_m+(not_free1(a1)+not_free2(a2))*w_b);
              #ifdef _TRACE_DEBUG
              cout<<"發生ＭＡＴＣＨ"<<endl;
                cout<<arc1[a1]<<endl;
                cout<<seq1[a1]<<endl;
                cout<<seq2[a2]<<endl;
                cout<<arc2[a2]<<endl;
                cout<<"match cost ="<<base_matching(a1,a2)*w_m+(not_free1(a1)+not_free2(a2))*w_b<<endl;
              #endif
                k--;
                l--;
            }
            else
                seqaln=false;

        }
	  // sequence alignment
      /*
	  while (seqaln)
	    {
	      int a1=l1+k-1;                      // a1,a2 sequence positions 
	      int a2=l2+l-1;
        #ifdef _TRACE_DEBUG
        cout<<"TREACEBACKing alignemnt"<<endl;
        cout<<"a1="<<a1<<" a2="<<a2<<endl;
        cout<<"k="<<k<<" l="<<l<<endl;
        cout<<"direct: "<<direct[k][l]<<endl;
        #endif

	      if (k==0 && l==0)
		seqaln=false;
	      else if (k>0 && fabs(M[k][l]-(M[k-1][l]+w_d+not_free1(a1)*(0.5*w_r-w_d)))<eps )
		{
            if(l!=0||l2!=0)
            {

              insert(ali,a1,-1,w_d+not_free1(a1)*(0.5*w_r-w_d));
             del_score=w_d+not_free1(a1)*(0.5*w_r-w_d);
            }
            else
              insert(ali,a1,-1,0);
		  k--;
		}
	      else if (l>0 && fabs(M[k][l]-(M[k][l-1]+w_d+not_free2(a2)*(0.5*w_r-w_d)))<eps )
		{
            if(k!=0||l1!=0)
            {
                  insert(ali,-1,a2,w_d+not_free2(a2)*(0.5*w_r-w_d));
                    del_score+=w_d+not_free2(a2)*(0.5*w_r-w_d);    
            }
            else
                  insert(ali,-1,a2,0);

		  l--;
		}
	      else if (k>0 && l>0 && fabs(M[k][l]-(M[k-1][l-1]+base_matching(a1,a2)*w_m+(not_free1(a1)+not_free2(a2))*w_b))<eps)
		{
		  insert(ali,a1,a2,base_matching(a1,a2)*w_m+(not_free1(a1)+not_free2(a2))*w_b);
            match_score+=base_matching(a1,a2)+(not_free1(a1)+not_free2(a2)*w_b);
		  k--;
		  l--;
		
        }
            //發現若來源並非上述三種
	      else
		seqaln=false;
	    }
        */
	  int a1=l1+k-1;                      // a1,a2 sequence positions 
	  int a2=l2+l-1;                      // right arc ends

	  // base-pair alignment
	  if (arc1[a1]==')' && arc2[a2]==')')
	    {
          arc_flag = true;
	      double w=M[L1[I1[a1]]-l1][L2[I2[a2]]-l2]+D[I1[a1]][I2[a2]]+arc_operation(L1[I1[a1]],L2[I2[a2]],a1,a2);//(base_matching(L1[I1[a1]],L2[I2[a2]])+base_matching(a1,a2))*0.5*w_am;
	      if (fabs(M[k][l]-w)<eps)
		{
		  int i1=L1[I1[a1]];              // left arc ends
		  int j1=L2[I2[a2]];
		  
		  double edge_weight=0.5*arc_operation(L1[I1[a1]],L2[I2[a2]],a1,a2); //(base_matching(L1[I1[a1]],L2[I2[a2]])+base_matching(a1,a2))*0.5*w_am;
          arc_score+=2*edge_weight;
		  insert(ali,i1,j1,edge_weight);
		  insert(ali,a1,a2,edge_weight);
		   
		  four_tuple CR1,CR2;
		  CR1.l1=l1   ; CR1.r1=i1-1 ; CR1.l2=l2   ; CR1.r2=j1-1 ;
		  CR2.l1=i1+1 ; CR2.r1=a1-1 ; CR2.l2=j1+1 ; CR2.r2=a2-1 ;
          #ifdef _TRACE_DEBUG
          cout<<"發生ARC operation"<<endl;
          cout<<arc1[i1]<<"\t"<<arc1[a1]<<endl;
          cout<<seq1[i1]<<"\t"<<seq1[a1]<<endl;
          cout<<seq2[j1]<<"\t"<<seq2[a2]<<endl;
          cout<<arc2[j1]<<"\t"<<arc2[a2]<<endl;
          cout<<"arc cost="<<2*edge_weight<<endl;
          cout<<"之後處理 CR1.r1="<<i1-1<<" CR1.l1="<<l1<<" CR1.r2="<<j1-1<<" CR1.l2="<<l2<<endl;
          cout<<"優先處理 CR2.r1="<<a1-1<<" CR2.l1="<<i1+1<<" CR2.r2="<<a2-1<<" CR2.l2="<<j1+1<<endl;
         #endif
		  ranges.push(CR1);
		  ranges.push(CR2);
		}
	    }
	}
    }

    
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
        #ifdef _TRACE_DEBUG
        printf("Successful open file!\n");
        puts(file_name);
        #endif
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
    /*對I1 I2兩個vector配置arc1 and arc2大小的記憶體空間*/
    vector<vector <string> > direct; 
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
    cout<<"seq1's arc"<<endl;
    for(int i = 0;i<L1.size();i++)
        cout<<L1[i]<<"\t"<<R1[i]<<endl;
    cout<<"total base pair 1= "<<L1.size()<<endl;
    cout<<"seq2's arc"<<endl;
    for(int i = 0;i<L2.size();i++)
        cout<<L2[i]<<"\t"<<R2[i]<<endl;
    cout<<"total base pair 2= "<<L2.size()<<endl;
    for(int i = 0;i<L1.size();i++)
        cout<<L1[i]<<"\t"<<R1[i]<<"\t"<<I1[i]<<endl;
    cout<<endl;
    for(int i = 0;i<L2.size();i++)
        cout<<L2[i]<<"\t"<<R2[i]<<"\t"<<I2[i]<<endl;
#endif
    /*initialize*/
    D.resize(L1.size());
    //lower.resize(L1.size());
    //upper.resize(L1.size());
    //cout<<"D matrix szie = "<<D.size()<<endl;

    
    for(int i = 0;i<D.size();i++)
    {
        D[i].resize(L2.size());
        //lower.resize(L2.size());
        //upper.resize(L2.size());
    }

    for(int i = 0;i<L1.size();i++)
    {
        for(int j = 0;j<L2.size();j++)
        {
            //added affine gap penaly matrix : lower and upper
            M.resize(R1[i]-L1[i]);
            direct.resize(R1[i]-L1[i]);
            for(int s = 0;s<M.size();s++)
            {
                M[s].resize(R2[j]-L2[j]);
                direct[s].resize(R2[j]-L2[j]);
                
            }

            M[0][0]=0;
            direct[0][0]="O";
            for(int k = 1;k<R1[i]-L1[i];k++)
            {
                M[k][0]=M[k-1][0]+w_d+not_free1(L1[i]+k)*(0.5*w_r-w_d);
               direct[k][0]="B";

            }
            for (int l=1;l<R2[j]-L2[j];l++)
            {
                M[0][l]=M[0][l-1]+w_d+not_free2(L2[j]+l)*(0.5*w_r-w_d);
               direct[0][l]="C";
            }
        //compute M
        double v1,v2,v3,v4;
        for(int k=1;k<R1[i]-L1[i];k++)
            for(int l=1;l<R2[j]-L2[j];l++)
            {
                //max 
                v1 = v2 = v3 = v4 = -10000;
                //min                 v1 = v2 = v3 = v4 = 10000;
                int a1 = L1[i]+k;
                int a2 = L2[j]+l;
                //edit operation
                //v1=M[k-1][l]+w_d+not_free1(a1)*(0.5*w_r-w_d);
                //
                v1=M[k-1][l]+w_d+not_free1(a1)*(0.5*w_r-w_d);
                v2=M[k][l-1]+w_d+not_free2(a2)*(0.5*w_r-w_d);
                //
                v3=M[k-1][l-1]+base_matching(a1,a2)*w_m+(not_free1(a1)+not_free2(a2))*w_b;
                if(arc1[a1]==')' && arc2[a2]==')')
                {
                    int leftpoint = L1[I1[a1]];
                    int leftpoint2 = L2[I2[a2]];
                    v4 = M[leftpoint - L1[i]-1][leftpoint2 - L2[j]-1]+D[I1[a1]][I2[a2]]+arc_operation(L1[I1[a1]],L2[I2[a2]],R1[I1[a1]],R2[I2[a2]]);
                    //(base_matching(L1[I1[a1]],L2[I2[a2]])+base_matching(R1[I1[a1]],R2[I2[a2]]))*0.5*w_am;
                }
                //M[k][l]=min4(v1,v2,v3,v4);
                M[k][l]=max4(v1,v2,v3,v4);
             if(M[k][l]==v1)
                direct[k][l]="B";
             else if(M[k][l]==v2)
                direct[k][l]="C";
             else if(M[k][l]==v3)
                direct[k][l]="A";
             else if(M[k][l]==v4)
                direct[k][l]="S";

                }
        D[i][j]=M[R1[i]-L1[i]-1][R2[j]-L2[j]-1];

        
        #ifdef _TRACE_DEBUG
        cout<<"ARC1 base pair"<<endl;
        cout<<R1[i]<<"\t"<<L1[i]<<endl;
        cout<<"ARC2 base pair"<<endl;
        cout<<R2[j]<<"\t"<<L2[j]<<endl;
         
        cout<<"COMPUTATION\n\t";
        for(int count =0;count<R2[j]-L2[j];count++)
        {
            cout<<"a2:"<<count<<"\t";
        }
        cout<<"\t";
        cout<<endl;
        for(int count =0; count<R1[i]-L1[i];count++)
        {
            cout<<"a1:"<<count<<"\t";
            for(int count2=0;count2<R2[j]-L2[j];count2++)
                cout<<M[count][count2]<<"\t";
            cout<<endl;
        }       
        cout<<"Direct array"<<endl;
        cout<<"\t";
        for(int count =0;count<R2[j]-L2[j];count++)
        {
            cout<<"a2:"<<count<<"\t";
        }
        cout<<"\t";
        cout<<endl;
        for(int count =0; count<R1[i]-L1[i];count++)
        {
            cout<<"a1:"<<count<<"\t";
            for(int count2=0;count2<R2[j]-L2[j];count2++)
                cout<<direct[count][count2]<<"\t";
            cout<<endl;
        }
        #endif
        
        }
    }
    
    direct.resize(arc1.size()+1);
    M.resize(arc1.size()+1);
    for(int i =0; i<M.size();i++)
    {
        M[i].resize(arc2.size()+1);
        direct[i].resize(arc2.size()+1);
    }

    M[0][0]=0;
    direct[0][0]="o";
    //這裏應該初始值設定成0
    for(int k=1;k<=arc1.size();k++)
    {
        //M[k][0]=M[k-1][0]+w_d+not_free1(k-1)*(0.5*w_r-w_d);
        M[k][0]=0;
        direct[k][0]="B";

    }
    for(int l = 1;l<=arc2.size();l++)
    {
        //M[0][l]=M[0][l-1]+w_d+not_free2(l-1)*(0.5*w_r-w_d);
        M[0][l]=0;
        direct[0][l]="C";
        
    }

    //compute M
    double v1,v2,v3,v4;
    for (int k=1;k<=arc1.size();k++)
        for (int l=1;l<=arc2.size();l++)
        {
                //min v1 = v2 = v3 = v4 = 10000;
            //max
            v1 = v2 = v3 = v4 = -10000;
            v1=M[k-1][l]+w_d+not_free1(k-1)*(0.5*w_r-w_d);
            v2=M[k][l-1]+w_d+not_free2(l-1)*(0.5*w_r-w_d);
            v3=M[k-1][l-1]+base_matching(k-1,l-1)*w_m+(not_free1(k-1)+not_free2(l-1))*w_b;
            if(arc1[k-1]==')'&& arc2[l-1]==')')
            {
                int leftpoint = L1[I1[k-1]];
                int leftpoint2 = L2[I2[l-1]];
                v4=M[leftpoint][leftpoint2]+D[I1[k-1]][I2[l-1]]+arc_operation(L1[I1[k-1]],L2[I2[l-1]],R1[I1[k-1]],R2[I2[l-1]]);
                //(base_matching(L1[I1[k-1]],L2[I2[l-1]])+base_matching(R1[I1[k-1]],R2[I2[l-1]]))*0.5*w_am;

            }
            //M[k][l]=min4(v1,v2,v3,v4);
            M[k][l]=max4(v1,v2,v3,v4);
             if(M[k][l]==v1)
                direct[k][l]="B";
             else if(M[k][l]==v2)
                direct[k][l]="C";
             else if(M[k][l]==v3)
                direct[k][l]="A";
             else if(M[k][l]==v4)
                direct[k][l]="S";
                
        }
    #ifdef _TRACE_DEBUG
    cout<<"COMPUTATION"<<endl;
    cout<<"\t";
    for(int count = 0;count<=arc2.size();count++)
    {
        cout<<"a2:"<<count<<"\t";
    }
    cout<<endl;
    for(int count = 0;count<=arc1.size();count++)
    {
        cout<<"a1:"<<count<<"\t";
        for(int count2=0;count2<=arc2.size();count2++)
        {
            cout<<left<<M[count][count2]<<"\t";
        }
        cout<<endl;
    }
    cout<<"Direct"<<endl;
    cout<<"\t";
    for(int count = 0;count<=arc2.size();count++)
    {
        cout<<"a2:"<<count<<"\t";
    }
    cout<<endl;
    for(int count = 0;count<=arc1.size();count++)
    {
        cout<<"a1:"<<count<<"\t";
        for(int count2=0;count2<=arc2.size();count2++)
        {
            cout<<left<<direct[count][count2]<<"\t";
        }
        cout<<endl;
    }
    cout<<"done!"<<endl;

    #endif
    //找最大值
    max_column_value=M[0][arc2.size()];
    #ifdef _TRACE_DEBUG
    for(int count = 0;count<=arc1.size();count++)
        cout<<M[count][arc2.size()]<<"\t";
    cout<<endl;
    #endif
    for(int count = 0;count<=arc1.size();count++)
    {
        if(max_column_value<M[count][arc2.size()])
        {
            max_column_value=M[count][arc2.size()];
            current_column=count;
        }
    }
    max_row_value=M[arc1.size()][0];
    #ifdef _TRACE_DEBUG
    for(int count = 0;count<=arc2.size();count++)
        cout<<M[arc1.size()][count]<<"\t";
    cout<<endl;
    #endif
    for(int count =0;count<=arc2.size();count++)
    {
        if(max_row_value<M[arc1.size()][count])
        {
            max_row_value=M[arc1.size()][count];
            current_row=count;
        }
    }
    
    //回傳last column or last row's max value
    //return M[arc1.size()][arc2.size()];
    #ifdef _TRACE_DEBUG
        cout<<"max_column_value : "<<max_column_value<<endl;
        cout<<"max_row_value : "<<max_row_value<<endl;
    #endif
    if(max_row_value>max_column_value)
        return M[arc1.size()][current_row];
    else
        return M[current_column][arc2.size()];

}



