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
#include <math.h>
#include <stack>
using namespace std;

//#define _DEBUG 
#define _COM_DEBUG

/*alignment score*/
typedef struct alignment{
    int p1,p2;
    double weight;
    alignment *next;
}Alignment;
/*紀錄位置*/
typedef struct fourtuple
{
    int leftpoint1;
    int rightpoint1;
    int leftpoint2;
    int rightpoint2;
}FourTuple;



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
vector<double> weights;
vector<int>    L1,R1,I1,L2,R2,I2;
stack<int>     str_stack;
vector<vector<double> > M;
vector<vector<double> > D;
const double eps=0.0000001;
static int **scoring_matrix;

double w_d =-1.0;  // base deletion
double w_m =-1.0;  // base mismatch
double w_r =-2.0;  // arc  removing
double w_b =-1.5;  // arc  breaking
double w_am=-1.8;  // arc  mismatc

int not_free1    (int pos)          { return (arc1[pos]=='.' ? 0:1)      ;  }
int not_free2    (int pos)          { return (arc2[pos]=='.' ? 0:1)      ;  }
/*沒這麼簡單了*/
int scoring_base_mismatch(int pos1,int pos2){ return (seq1[pos1]!=seq2[pos2]?1:0);  }
int base_mismatch(int,int);

int BinarySearch(char ,IndexMatrix *,int );
void traceback();
void insert(Alignment*,int,int,double);
int** readmat(char *);
void read_data(const char *,const char *);
double computation();
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
double max4(double a,double b,double c,double d)
{

  if (a>=b && a>=c && a>=d)
    return a;
  else if (b>=a && b>=c && b>=d)
    return b;
  else if (c>=a && c>=b && c>=d)
    return c;
  return d;
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
    

    char path[100];
    //sprintf(path,"./SM/BLOSUM-like_scoring_matrix");
    sprintf(path,"./SM/iPARTS2_new_23C_4L_matrix");
    scoring_matrix=readmat(path);
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
    read_data("./test_file","./result"); 

    /*Computation*/
    double score = computation();
    cout<<score<<endl;
    traceback();
    /*釋放記憶體*/
    free(alphabet_index);//來源 line:33
    return 0;
}

void insert(Alignment* start,int pos1,int pos2,double weight)
{
    //Construction of the alignment
    Alignment* insertion = new Alignment;
    insertion->p1     = pos1;
    insertion->p2     = pos2;
    insertion->weight = weight;
    insertion->next   = NULL; 
    Alignment* iter   = start;
    if (pos1==-1 && pos2!=-1)
        while (iter->next!=NULL && pos2>iter->next->p2)
            iter=iter->next;
    else if (pos2==-1 && pos1!=-1)
        while (iter->next!=NULL && pos1>iter->next->p1)
            iter=iter->next;
    else
        while (iter->next!=NULL && pos2>iter->next->p2 && pos1>iter->next->p1 )
            iter=iter->next;
      
    insertion->next=iter->next;
    iter->next=insertion;
}

void traceback()
{
    Alignment* ali = new Alignment;
    ali->p1=-1;
    ali->p2=-1;
    ali->next=NULL;
    stack<double> weight;
    double v1,v2,v3,v4;

    //range is the currently computed sequence range
    FourTuple range;
    range.leftpoint1=0;
    range.rightpoint1=arc1.size()-1;
    range.leftpoint2=0;
    range.rightpoint2=arc2.size()-1;

    stack<FourTuple> range_stack;
    range_stack.push(range);
    int round = 0;
    while(!range_stack.empty())
    {
        int leftpoint1 = range_stack.top().leftpoint1;
        int rightpoint1 = range_stack.top().rightpoint1;
        int leftpoint2 = range_stack.top().leftpoint2;
        int rightpoint2 = range_stack.top().rightpoint2;
        range_stack.pop();

        if(leftpoint1 > rightpoint1 && leftpoint2<=rightpoint2)
            for(int count = rightpoint2;count>=leftpoint2;--count)
                insert(ali,-1,count,w_d);
        else if(leftpoint1 <= rightpoint1 && leftpoint2> rightpoint2)
            for(int count = rightpoint1;count>=leftpoint1;--count)
                insert(ali,count,-1,w_d);
        else if(leftpoint1<=rightpoint1 && leftpoint2<=rightpoint2)
        {
            M.resize(rightpoint1-leftpoint1+2);
            for(int count = 0;count<M.size();--count)
            {
                M[count].resize(rightpoint2-leftpoint2+2);
                for(int inner_count=0;inner_count<M[count].size();inner_count++)
                    M[count][inner_count]=0;
            }
            /*initalize*/
            M[0][0]=0;
            for(int k=1;k<rightpoint1-leftpoint1+2;k++)
                M[k][0]=M[k-1][0]+w_d+not_free1(leftpoint1+k-1)*(0.5*w_r-w_d);
            for(int l=1;l<rightpoint2-leftpoint2+2;l++)
                M[0][l]=M[0][l-1]+w_d+not_free2(leftpoint2+l-1)*(0.5*w_r-w_d);
            //compute M
            for(int k = 1;k<rightpoint1-leftpoint1+2;++k)
                for(int l=1;l<rightpoint2-leftpoint2+2;++l)
                {
                    v1=v2=v3=v4=-99999;
                    int offset_index1=leftpoint1+k-1;
                    int offset_index2=leftpoint2+l-1;
                    //case1
                    v1=M[k-1][l]+w_d+not_free1(offset_index1)*(0.5*w_r-w_d);
                    //case2
                    v2=M[k][l-1]+w_d+not_free2(offset_index1)*(0.5*w_r-w_d);
                    //case3
                    v3=M[k-1][l-1]+base_mismatch(offset_index1,offset_index2)+(not_free1(offset_index1)+not_free2(offset_index2))*0.5*w_b;
                    //case4
                    if(arc1[offset_index1]==')' && arc2[offset_index2]==')')
                    {
                        int current_basepair_leftpoint1=L1[I1[offset_index1]];
                        int current_basepair_leftpoint2=L2[I2[offset_index2]];
                        v4 = M[current_basepair_leftpoint1-leftpoint1][current_basepair_leftpoint2-leftpoint2]+
                        D[I1[offset_index1]][I2[offset_index2]]+(base_mismatch(current_basepair_leftpoint1,current_basepair_leftpoint2)+base_mismatch(offset_index1,offset_index2))*0.5*w_am;
                    }
                    M[k][l]=max4(v1,v2,v3,v4);

                }
            bool seq_flag = true;

            
        }
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

int base_mismatch(int pos1,int pos2)
{
    int index_x,index_y;
    index_x = BinarySearch(seq1[pos1],alphabet_index,size);
    index_y = BinarySearch(seq2[pos2],alphabet_index,size);
    //cout<<seq1[pos1]<<" vs "<<seq2[pos2]<<" "<<scoring_matrix[index_x][index_y]<<endl;
    return scoring_matrix[index_x][index_y];  
}
/*讀取序列資訊*/
void read_data(const char *path,const char *outpath)
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
    cout<<sub<<" length="<<sub.size()<<" seq length="<<sub.size()-9<<endl;

    /*read arc1*/    
    getline(file,sub,'\n');
    temp_size = sub.size();
    arc1 = sub.substr(7,temp_size-9);
    cout<<sub<<" length="<<sub.size()<<" arc length="<<sub.size()-9<<endl;

    /*read seq2*/    
    getline(file,sub,'\n');
    temp_size = sub.size();
    seq2 = sub.substr(7,temp_size-9);
    cout<<sub<<" length="<<sub.size()<<" seq length="<<sub.size()-9<<endl;

    /*read arc2*/    
    getline(file,sub,'\n');
    temp_size = sub.size();
    arc2 = sub.substr(7,temp_size-9);
    cout<<sub<<" length="<<sub.size()<<" arc length="<<sub.size()-9<<endl;

    /*read matrixpath*/    
    getline(file,sub,'\n');
    temp_size = sub.size();
    matrixpath = sub.substr(10,temp_size-12);
    cout<<sub<<endl;

    /*read gap_opp*/    
    getline(file,sub,'\n');
    temp_size = sub.size();
    sub = sub.substr(6,temp_size-8);
    gap_opp = atoi(sub.c_str());
    cout<<sub<<endl;

    /*read gap_exp*/    
    getline(file,sub,'\n');
    temp_size = sub.size();
    sub = sub.substr(6,temp_size-8);
    gap_exp = atoi(sub.c_str());
    cout<<sub<<endl;


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
#ifdef _COM_DEBUG
    for(int i = 0;i<L1.size();i++)
        cout<<L1[i]<<" "<<R1[i]<<" "<<I1[i]<<endl;
    cout<<endl;
    for(int i = 0;i<L2.size();i++)
        cout<<L2[i]<<" "<<R2[i]<<" "<<I2[i]<<endl;
#endif
    /*initialize*/
    D.resize(L1.size());
    for(int i = 0;i<D.size();i++)
        D[i].resize(L2.size());
    for(int i = 0;i<L1.size();i++)
    {
        for(int j = 0;j<L2.size();j++)
        {
            M.resize(R1[i]-L1[i]);
            for(int s = 0;s<M.size();s++)
                M[s].resize(R2[j]-L2[j]);
            M[0][0]=0;
            for(int k = 1;k<R1[i]-L1[i];k++)
            {
                M[k][0]=M[k-1][0]+w_d+not_free1(L1[i]+k)*(0.5*w_r-w_d);
            }
            for (int l=1;l<R2[j]-L2[j];l++)
            {
                M[0][l]=M[0][l-1]+w_d+not_free2(L2[j]+l)*(0.5*w_r-w_d);
            }
        //compute M
        double v1,v2,v3,v4;
        for(int k=1;k<R1[i]-L1[i];k++)
            for(int l=1;l<R2[j]-L2[j];l++)
            {
                v1 = v2 = v3 = v4 = -100000;
                int a1 = L1[i]+k;
                int a2 = L2[j]+l;

                v1=M[k-1][l]+w_d+not_free1(a1)*(0.5*w_r-w_d);
                v2=M[k][l-1]+w_d+not_free2(a2)*(0.5*w_r-w_d);
                v3=M[k-1][l-1]+base_mismatch(a1,a2)+(not_free1(a1)+not_free2(a2))*0.5*w_b;
                if(arc1[a1]==')' && arc2[a2]==')')
                {
                    int leftpoint = L1[I1[a1]];
                    int leftpoint2 = L2[I2[a2]];
                    v4 = M[leftpoint - L1[i]-1][leftpoint2 - L2[j]-1]+D[I1[a1]][I2[a2]]+(base_mismatch(L1[I1[a1]],L2[I2[a2]])+base_mismatch(R1[I1[a1]],R2[I2[a2]]))*0.5*w_am;
                }
                //M[k][l]=min4(v1,v2,v3,v4);
                M[k][l]=max4(v1,v2,v3,v4);
            }
        D[i][j]=M[R1[i]-L1[i]-1][R2[j]-L2[j]-1];

        }
    }

    M.resize(arc1.size()+1);
    for(int i =0; i<M.size();i++)
        M[i].resize(arc2.size()+1);

    M[0][0]=0;

    for(int k=1;k<=arc1.size();k++)
        M[k][0]=M[k-1][0]+w_d+not_free1(k-1)*(0.5*w_r-w_d);
    for(int l = 1;l<=arc2.size();l++)
        M[0][l]=M[0][l-1]+w_d+not_free2(l-1)*(0.5*w_r-w_d);

    //compute M
    double v1,v2,v3,v4;
    for (int k=1;k<=arc1.size();k++)
        for (int l=1;l<=arc2.size();l++)
        {
            v1=v2=v3=v4=-10000;
            v1=M[k-1][l]+w_d+not_free1(k-1)*(0.5*w_r-w_d);
            v2=M[k][l-1]+w_d+not_free2(l-1)*(0.5*w_r-w_d);
            v3=M[k-1][l-1]+base_mismatch(k-1,l-1)+(not_free1(k-1)+not_free2(l-1))*0.5*w_b;
            if(arc1[k-1]==')'&& arc2[l-1]==')')
            {
                int leftpoint = L1[I1[k-1]];
                int leftpoint2 = L2[I2[l-1]];
                v4=M[leftpoint][leftpoint2]+D[I1[k-1]][I2[l-1]]
                +(base_mismatch(L1[I1[k-1]],L2[I2[l-1]])+base_mismatch(R1[I1[k-1]],R2[I2[l-1]]))*0.5*w_am;

            }
            //M[k][l]=min4(v1,v2,v3,v4);
            M[k][l]=max4(v1,v2,v3,v4);
        }

    return M[arc1.size()][arc2.size()];
}
