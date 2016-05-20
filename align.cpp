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
//#define _COM_DEBUG

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
double arc_max = 0.0,seq_max = 0.0;
vector<double> weights;
vector<int>    L1,R1,I1,L2,R2,I2;
stack<int>     str_stack;
vector<vector<double> > M;
vector<vector<double> > lower;
vector<vector<double> > upper;
vector<vector<double> > middle;

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

double k = 2;//arc match 
double l= 0.1;//arc mismatch
double m= 1.5;//when one base score >=0 and other base score < 0, m 不能比k大

double number=1;

double w_d =-1*number;  // base deletion
double w_r =-2*number;  // arc  removing
double w_b =-1.5*number;  // arc  breaking

double common_opp = 9;
double common_exp = 1;
double base_opp = 9;  // base opening gap
double base_epp = 1;  // base extension gap
double arc_opp = 18;   // arc opening gap
double arc_exp = 2;   // arc extension gap_opp
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
        return k*(base_matching(p1,p2)+base_matching(p3,p4));
    }
    //當alignment score <0 views as arc-mismatch
    else
    {
        if((base_matching(p1,p2)<0 && base_matching(p3,p4)>=0)||(base_matching(p1,p2)>=0&&base_matching(p3,p4)<0))
        {
            if((base_matching(p1,p2)+base_matching(p3,p4))>=0)
            {

                return (base_matching(p1,p2)+base_matching(p3,p4))*m;
            } 
            else
            {
                return (base_matching(p1,p2)+base_matching(p3,p4))*(l/2);
            }
        }
        else if(base_matching(p1,p2)<=0 && base_matching(p3,p4)<0)
        //return (base_matching(p1,p2)+base_matching(p3,p4))*0.5*w_am;
            return (base_matching(p1,p2)+base_matching(p3,p4))*l;
    }
}



int BinarySearch(char ,IndexMatrix *,int );
void traceback();
void insert(alignment*,int,int,double);
int** readmat(char *);
void read_data(const char *);
void write_data(double,const char *);
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
    cout<<argc<<endl;
    if(argc!=7)
    {
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/test_data/1A9N_Q_to_1E7K_C/semi_input.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1M90_B_to_1NKW_9/semi_input.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1IBK_A_to_1N33_A/semi_input.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1BAU_A_to_1S9S_A/semi_input.php";
        pdb_compare_path= "/home/watchlee/Research_Programming/X3DNA/23-4L_SARA_FSCOR_structure/1BAU_A_to_1S9S_A/semi_input.php";
        //pdb_compare_path= "/home/watchlee/Research_Programming/AlignmentToHen/Test_for_GLOBAL/1IBK_A_to_1N32_A.php";
        pdb_result= "/home/watchlee/result.php";
        number = 1;
        k = 2;
        l = 0.1;
        m = 1.5;
    }
    else
    {
       pdb_compare_path=argv[1] ;
       pdb_result= argv[2];
       number = atoi(argv[3]);
       k = atoi(argv[4]);
       l = atoi(argv[5]);
       m = atoi(argv[6]);
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

    cout<<"score="<<test_alignment()<<endl;
    
    write_data(score,pdb_result);

    /*釋放記憶體*/
    free(alphabet_index);//來源 line:33
    return 0;
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
    temp_D.resize(arc1.size()+1);
    temp_middle.resize(arc1.size()+1);
    temp_lower.resize(arc1.size()+1);
    temp_upper.resize(arc1.size()+1);
    for(int count = 0;count<=arc1.size();count++)
    {
        temp_D[count].resize(arc2.size()+1);
        temp_middle[count].resize(arc2.size()+1);
        temp_lower[count].resize(arc2.size()+1);
        temp_upper[count].resize(arc2.size()+1);
    }
    temp_middle[0][0]=temp_lower[0][0]=temp_upper[0][0]=0.0;
    /*某些部分在affine gap penalty中entry不會使用到設定為負無限大當作不會參考*/
    /*對Y軸進行初始*/
    for(int count = 1; count<=arc1.size();count++)
    {
        temp_lower[count][0]=-LARGE_NUMBER;
        temp_middle[count][0]=-LARGE_NUMBER;
        temp_D[count][0]=temp_upper[count][0]=-gap_opp-(count*gap_exp);

    }
    /*對X軸進行初始*/
    for(int count = 1; count<=arc2.size();count++)
    {
        temp_middle[0][count]=-LARGE_NUMBER;
        temp_upper[0][count]=-LARGE_NUMBER;
        temp_D[0][count]=temp_lower[0][count]=-gap_opp-(count*gap_exp);
    }
    for(int count = 1;count<=arc1.size();count++)
    {
        for(int count2 = 1;count2<=arc2.size();count2++)
        {
            temp_middle[count][count2]=max(temp_middle[count-1][count2-1],temp_lower[count-1][count2-1],temp_upper[count-1][count2-1])+base_matching(count-1,count2-1);

            temp_lower[count][count2]=max(temp_lower[count][count2-1]-gap_exp,temp_middle[count][count2-1]-gap_opp-gap_exp,temp_upper[count][count2-1]-gap_opp-gap_exp);
            //
            temp_upper[count][count2]=max(temp_upper[count-1][count2]-gap_exp,temp_middle[count-1][count2]-gap_opp-gap_exp,temp_lower[count-1][count2]-gap_opp-gap_exp);
            //arc annotated sequence must consider four tuples lower,upper,middle,special
            temp_D[count][count2]=max(temp_middle[count][count2],temp_lower[count][count2],temp_upper[count][count2]);
            
        }
    }

    int i = arc1.size(),j = arc2.size(); 
    string s1=seq1,s2=seq2;
    while(i>0||j>0)
    {
        if(temp_D[i][j]==temp_lower[i][j])
        {
            s1.insert(i,1,'-');
            j--;
        }
        else if(temp_D[i][j]==temp_upper[i][j])
        {
            s2.insert(j,1,'-');
            i--;
        }
        else
        {

            i--;
            j--;
        }
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

  // range is the currently computed sequence range
  four_tuple range;
  range.l1=0;
  range.l2=0;
  range.r1=seq1.size()-1;
  range.r2=seq2.size()-1;

  stack<four_tuple> ranges;
  ranges.push(range);
  
  while(!ranges.empty())
    {
      int l1=ranges.top().l1;
      int r1=ranges.top().r1;
      int l2=ranges.top().l2;
      int r2=ranges.top().r2;
      ranges.pop();

      if (l1>r1 && l2<=r2)
	for (int s=r2;s>=l2;s--)
	  insert(ali,-1,s,w_d);
      else if (l1<=r1 && l2>r2)
	for (int s=r1;s>=l1;s--)
	  insert(ali,s,-1,w_d);
      else if (l1<=r1 && l2<=r2)
	{
	  // init and compute M
	  M.resize(r1-l1+2);
	  lower.resize(r1-l1+2);
	  middle.resize(r1-l1+2);
	  upper.resize(r1-l1+2);
	  for(int s=0;s<M.size();s++)
        {
	      M[s].resize(r2-l2+2);
	      middle[s].resize(r2-l2+2);
	      lower[s].resize(r2-l2+2);
	      upper[s].resize(r2-l2+2);
	      for(int t=0;t<M[s].size();t++)
            M[s][t]=0;
	    }
	  /*---------affine gap penalty ------2016/5/20*/ 
	  M[0][0]=middle[0][0]=lower[0][0]=upper[0][0]=0;
	  for (int k=1;k<r1-l1+2;k++)
        {

        lower[k][0]=-LARGE_NUMBER; 
        middle[k][0]=-LARGE_NUMBER; 
	    M[k][0]=M[k-1][0]+w_d+not_free1(l1+k-1)*(0.5*w_r-w_d);
        }
	  
	  for (int l=1;l<r2-l2+2;l++)
        {
        upper[0][l]=-LARGE_NUMBER;
        middle[0][l]=-LARGE_NUMBER;
	    M[0][l]=M[0][l-1]+w_d+not_free2(l2+l-1)*(0.5*w_r-w_d);
        }
	  
	  for (int k=1;k<r1-l1+2;k++)
	    for (int l=1;l<r2-l2+2;l++)
	      {
		//max
        v1=v2=v3=v4=-LARGE_NUMBER;
        //min 		v1=v2=v3=v4=LARGE_NUMBER;
		int a1=l1+k-1;                      // a1,a2 sequence positions 
		int a2=l2+l-1;
		
        gappenalty2 = -base_epp+(not_free2(a2)*(-arc_exp/2+base_epp));
        gappenalty1 = -base_epp+(not_free1(a1)*(-arc_exp/2+base_epp));
        open_gappenalty2 = -base_opp-base_epp+(not_free2(a2)*(-(arc_opp/2)-(arc_exp/2)+base_opp+base_epp));
        open_gappenalty1 = -base_opp-base_epp+(not_free1(a1)*(-(arc_opp/2)-(arc_exp/2)+base_opp+base_epp));
        v1=lower[k][l]=max(lower[k][l-1]+gappenalty2,middle[k][l-1]+open_gappenalty2,upper[k][l-1]+gappenalty2);
        v2=upper[k][l]=max(upper[k-1][l]+gappenalty1,middle[k-1][l]+open_gappenalty1,lower[k-1][l]+gappenalty1);
        v3=middle[k][l]=max(middle[k-1][l-1],lower[k-1][l-1],upper[k-1][l-1])+base_matching(a1,a2)+(not_free1(a1)+not_free2(a2))*0.5*w_b; 
		//v1=M[k-1][l]+w_d+not_free1(a1)*(0.5*w_r-w_d);
		//v2=M[k][l-1]+w_d+not_free2(a2)*(0.5*w_r-w_d);
		//v3=M[k-1][l-1]+base_matching(a1,a2)*w_m+(not_free1(a1)+not_free2(a2))*0.5*w_b;
		
		if (arc1[a1]==')' && arc2[a2]==')') 
		  {
		    int i1=L1[I1[a1]];
		    int j1=L2[I2[a2]];
		    //v4=middle[i1-l1][j1-l2]+D[I1[a1]][I2[a2]]+arc_operation(i1,j1,a1,a2);
		    v4=M[i1-l1][j1-l2]+D[I1[a1]][I2[a2]]+arc_operation(i1,j1,a1,a2);
		      //(base_matching(i1,j1)+base_matching(a1,a2))*0.5*w_am;
		  }
		//M[k][l]=min4(v1,v2,v3,v4);
		 M[k][l]=max4(v1,v2,v3,v4);
	      }
	  
	  bool seqaln=true;
	  int k=r1-l1+1;
	  int l=r2-l2+1;
	  // sequence alignment
	  while (seqaln)
	    {
	      int a1=l1+k-1;                      // a1,a2 sequence positions 
	      int a2=l2+l-1;
        gappenalty2 = -base_epp+(not_free2(a2)*(-arc_exp/2+base_epp));
        gappenalty1 = -base_epp+(not_free1(a1)*(-arc_exp/2+base_epp));
        open_gappenalty2 = -base_opp-base_epp+(not_free2(a2)*(-(arc_opp/2)-(arc_exp/2)+base_opp+base_epp));
        open_gappenalty1 = -base_opp-base_epp+(not_free1(a1)*(-(arc_opp/2)-(arc_exp/2)+base_opp+base_epp));

	      if (k==0 && l==0)
            seqaln=false;
	      //else if (k>0 && fabs(M[k][l]-(M[k-1][l]+w_d+not_free1(a1)*(0.5*w_r-w_d)))<eps )
	      else if (k>0 && fabs(M[k][l]-upper[k][l])<eps )
		  {
		  //insert(ali,a1,-1,w_d+not_free1(a1)*(0.5*w_r-w_d));
		  insert(ali,a1,-1,open_gappenalty1);
		  k--;
		}
	      else if (l>0 && fabs(M[k][l]-lower[k][l])<eps )
		{
		  //insert(ali,-1,a2,w_d+not_free2(a2)*(0.5*w_r-w_d));
		  insert(ali,-1,a2,open_gappenalty2);
		  l--;
		}
	      //else if (k>0 && l>0 && fabs(M[k][l]-(middle[k-1][l-1]+base_matching(a1,a2)*w_m+(not_free1(a1)+not_free2(a2))*0.5*w_b))<eps)
	      else if (k>0 && l>0 && fabs(M[k][l]-(M[k-1][l-1]+base_matching(a1,a2)*w_m+(not_free1(a1)+not_free2(a2))*0.5*w_b))<eps)
		{
		  insert(ali,a1,a2,base_matching(a1,a2)*w_m+(not_free1(a1)+not_free2(a2))*0.5*w_b);
		  k--;
		  l--;
		}
	      else
            seqaln=false;
	    }
	  
	  int a1=l1+k-1;                      // a1,a2 sequence positions 
	  int a2=l2+l-1;                      // right arc ends

	  // base-pair alignment
	  if (arc1[a1]==')' && arc2[a2]==')')
	    {
	      double w=M[L1[I1[a1]]-l1][L2[I2[a2]]-l2]+D[I1[a1]][I2[a2]]+arc_operation(L1[I1[a1]],L2[I2[a2]],a1,a2);//(base_matching(L1[I1[a1]],L2[I2[a2]])+base_matching(a1,a2))*0.5*w_am;
	      if (fabs(M[k][l]-w)<eps)
		{
		  int i1=L1[I1[a1]];              // left arc ends
		  int j1=L2[I2[a2]];
		  
		  double edge_weight=0.5*arc_operation(L1[I1[a1]],L2[I2[a2]],a1,a2); //(base_matching(L1[I1[a1]],L2[I2[a2]])+base_matching(a1,a2))*0.5*w_am;

		  insert(ali,i1,j1,edge_weight);
		  insert(ali,a1,a2,edge_weight);
		  
		  four_tuple CR1,CR2;
		  CR1.l1=l1   ; CR1.r1=i1-1 ; CR1.l2=l2   ; CR1.r2=j1-1 ;
		  CR2.l1=i1+1 ; CR2.r1=a1-1 ; CR2.l2=j1+1 ; CR2.r2=a2-1 ;
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
    cout<<"total base pair = "<<L1.size()<<endl;
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
    cout<<"total base pair 2= "<<L2.size()<<endl;
#ifdef _COM_DEBUG
    for(int i = 0;i<L1.size();i++)
        cout<<L1[i]<<" "<<R1[i]<<" "<<I1[i]<<endl;
    cout<<endl;
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
            middle.resize(R1[i]-L1[i]);

            M.resize(R1[i]-L1[i]);
            for(int s = 0;s<M.size();s++)
            {
                M[s].resize(R2[j]-L2[j]);
                //affine gap penalty 2016/5/20
                lower[s].resize(R2[j]-L2[j]);
                upper[s].resize(R2[j]-L2[j]);
                middle[s].resize(R2[j]-L2[j]);

            }

            M[0][0]=0;
            //affine gap penalty 2016/5/20
            lower[0][0]=upper[0][0]=0;

            for(int k = 1;k<R1[i]-L1[i];k++)
            {
                lower[k][0]=-LARGE_NUMBER; 
                middle[k][0]=-LARGE_NUMBER; 
                M[k][0]=upper[k][0]=-base_opp-(k*base_epp)+not_free1(L1[i]+k)*(0.5*(-arc_opp-k*arc_exp)+(base_opp+(k*base_epp)));
                //M[k][0]=M[k-1][0]+w_d+not_free1(L1[i]+k)*(0.5*w_r-w_d);

            }
            for (int l=1;l<R2[j]-L2[j];l++)
            {
                upper[0][l]=-LARGE_NUMBER;
                middle[0][l]=-LARGE_NUMBER;
                M[0][l]=lower[0][l]=-base_opp-(l*base_epp)+not_free1(L1[i]+l)*(0.5*(-arc_opp-l*arc_exp)+(base_opp+(l*base_epp)));
                //M[0][l]=M[0][l-1]+w_d+not_free2(L2[j]+l)*(0.5*w_r-w_d);
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
                double gappenalty2 = -base_epp+(not_free2(a2)*(-arc_exp/2+base_epp));
                double gappenalty1 = -base_epp+(not_free1(a1)*(-arc_exp/2+base_epp));
                double open_gappenalty2 = -base_opp-base_epp+(not_free2(a2)*(-(arc_opp/2)-(arc_exp/2)+base_opp+base_epp));
                double open_gappenalty1 = -base_opp-base_epp+(not_free1(a1)*(-(arc_opp/2)-(arc_exp/2)+base_opp+base_epp));
                v1=lower[k][l]=max(lower[k][l-1]+gappenalty2,middle[k][l-1]+open_gappenalty2,upper[k][l-1]+gappenalty2);
                v2=upper[k][l]=max(upper[k-1][l]+gappenalty1,middle[k-1][l]+open_gappenalty1,lower[k-1][l]+gappenalty1);
                v3=middle[k][l]=max(middle[k-1][l-1],lower[k-1][l-1],upper[k-1][l-1])+base_matching(a1,a2)+(not_free1(a1)+not_free2(a2))*0.5*w_b; 


                //v1=M[k-1][l]+w_d+not_free1(a1)*(0.5*w_r-w_d);
                //v2=M[k][l-1]+w_d+not_free2(a2)*(0.5*w_r-w_d);
                //v3=M[k-1][l-1]+base_matching(a1,a2)*w_m+(not_free1(a1)+not_free2(a2))*0.5*w_b;
                if(arc1[a1]==')' && arc2[a2]==')')
                {
                    int leftpoint = L1[I1[a1]];
                    int leftpoint2 = L2[I2[a2]];
                    v4 = M[leftpoint - L1[i]-1][leftpoint2 - L2[j]-1]+D[I1[a1]][I2[a2]]+arc_operation(L1[I1[a1]],L2[I2[a2]],R1[I1[a1]],R2[I2[a2]]);
                    //v4 = middle[leftpoint - L1[i]-1][leftpoint2 - L2[j]-1]+D[I1[a1]][I2[a2]]+arc_operation(L1[I1[a1]],L2[I2[a2]],R1[I1[a1]],R2[I2[a2]]);
                    //(base_matching(L1[I1[a1]],L2[I2[a2]])+base_matching(R1[I1[a1]],R2[I2[a2]]))*0.5*w_am;
                }
                //M[k][l]=min4(v1,v2,v3,v4);
                M[k][l]=max4(v1,v2,v3,v4);
            }
        D[i][j]=M[R1[i]-L1[i]-1][R2[j]-L2[j]-1];

        }
    }
    
    /*affine gap penalty 2016/5/20*/
    lower.resize(arc1.size()+1);
    middle.resize(arc1.size()+1);
    upper.resize(arc1.size()+1);

    M.resize(arc1.size()+1);
    for(int i =0; i<M.size();i++)
    {
        M[i].resize(arc2.size()+1);
        lower[i].resize(arc2.size()+1);
        middle[i].resize(arc2.size()+1);
        upper[i].resize(arc2.size()+1);

    }

    M[0][0]=lower[0][0]=middle[0][0]=upper[0][0]=0;
    /*initalize row*/
    for(int k=1;k<=arc1.size();k++)
    {
        lower[k][0]=-LARGE_NUMBER; 
        middle[k][0]=-LARGE_NUMBER; 
        M[k][0]=upper[k][0]=-base_opp-(k*base_epp)+not_free1(k-1)*(0.5*(-arc_opp-k*arc_exp)+(base_opp+(k*base_epp)));
        //M[k][0]=M[k-1][0]+w_d+not_free1(k-1)*(0.5*w_r-w_d);
    }
    /*initalize column*/
    for(int k=1;k<=arc1.size();k++)
    for(int l = 1;l<=arc2.size();l++)
    {
        upper[0][l]=-LARGE_NUMBER;
        middle[0][l]=-LARGE_NUMBER;
        M[0][l]=lower[0][l]=-base_opp-(l*base_epp)+not_free2(l-1)*(0.5*(-arc_opp-l*arc_exp)+(base_opp+(l*base_epp)));
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
            double gappenalty1 = -base_epp-(not_free1(k-1)*((arc_exp/2)-base_epp));
            double open_gappenalty1 = -base_opp-base_epp-(not_free1(k-1)*((arc_opp/2+arc_exp/2)-base_opp-base_epp));
            double gappenalty2 = -base_epp-(not_free2(l-1)*((arc_exp/2)-base_epp));
            double open_gappenalty2 = -base_opp-base_epp-(not_free2(l-1)*((arc_opp/2+arc_exp/2)-base_opp-base_epp));
            

            //origin version
            //v1=M[k-1][l]+w_d+not_free1(k-1)*(0.5*w_r-w_d);
            //v2=M[k][l-1]+w_d+not_free2(l-1)*(0.5*w_r-w_d);
            //v3=M[k-1][l-1]+base_matching(k-1,l-1)*w_m+(not_free1(k-1)+not_free2(l-1))*0.5*w_b;
            
            //affine gap penalty version 2016/5/20 
            v1=lower[k][l]=max(lower[k][l-1]+gappenalty2,middle[k][l-1]+open_gappenalty2,upper[k][l-1]+gappenalty2);
            v2=upper[k][l]=max(upper[k-1][l]+gappenalty1,middle[k-1][l]+open_gappenalty1,lower[k-1][l]+gappenalty1);
            //v1=lower[k][l]=max(lower[k][l-1]-common_exp,middle[k][l-1]-common_opp-common_exp,upper[k][l-1]-common_opp-common_exp);
            //v2=upper[k][l]=max(upper[k-1][l]-common_exp,middle[k-1][l]-common_opp-common_exp,lower[k-1][l]-common_opp-common_exp);
            v3=middle[k][l]=max(middle[k-1][l-1],lower[k-1][l-1],upper[k-1][l-1])+base_matching(k-1,l-1)+(not_free1(k-1)+not_free2(l-1))*0.5*w_b; 

            if(arc1[k-1]==')'&& arc2[l-1]==')')
            {
                int leftpoint = L1[I1[k-1]];
                int leftpoint2 = L2[I2[l-1]];
                //v4=middle[leftpoint][leftpoint2]+D[I1[k-1]][I2[l-1]]+arc_operation(L1[I1[k-1]],L2[I2[l-1]],R1[I1[k-1]],R2[I2[l-1]]);
                v4=M[leftpoint][leftpoint2]+D[I1[k-1]][I2[l-1]]+arc_operation(L1[I1[k-1]],L2[I2[l-1]],R1[I1[k-1]],R2[I2[l-1]]);
                //(base_matching(L1[I1[k-1]],L2[I2[l-1]])+base_matching(R1[I1[k-1]],R2[I2[l-1]]))*0.5*w_am;

            }
            //M[k][l]=min4(v1,v2,v3,v4);
            //M[k][l]=max4(v1,v2,v3,v4);
            M[k][l]=max4(lower[k][l],upper[k][l],middle[k][l],v4);
        }

    return M[arc1.size()][arc2.size()];
}



