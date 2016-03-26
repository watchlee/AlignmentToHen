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
#include<iostream>
using namespace std;

//#define _DEBUG 

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
            printf("%c %d",array[mid].alphabet,array[mid].pos);
            return array[mid].pos;
        }
        else if(array[mid].alphabet>character)
            high = mid-1;
        else if(array[mid].alphabet<character)
            low = mid+1;

    }
    return -1;
}

/*------------------------------------global variable-----------------------------------*/
IndexMatrix *temp_node=NULL;
static int size=0;

int** readmat(char *);
void read_data(const char *,const char *);
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
int main(int argc,char* argv[])
{
    

    char path[100];
    //sprintf(path,"./SM/BLOSUM-like_scoring_matrix");
    sprintf(path,"./SM/iPARTS2_new_23C_4L_matrix");
    int **scoring_matrix=readmat(path);
    int count = 0;
#ifdef _DEBUG 
    char test_char[100];
    for(count = 0;count<size;count++)
    {
        test_char[count]=temp_node[count].alphabet;
        printf("%d %c %d\n",temp_node[count].alphabet,temp_node[count].alphabet,temp_node[count].pos);
    }
    for(count = 0;count<size;count++)
        printf("%c ",test_char[count]);
    printf("\n");
#endif


    Merge_Sort sorting;
    sorting.Sort(temp_node,size);

    read_data("./test_file","./result"); 

    /*釋放記憶體*/
    free(temp_node);//來源 line:33
    return 0;
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
    char *seq1,*seq2,*arc1,*arc2,*matrixpath;
    int gap_opp,gap_exp;
    string sub;
    int temp_size;
    int option=0;
    /*read <?php*/
    getline(file,sub,'\n');
    cout<<sub<<endl;

    getline(file,sub,'\n');
    temp_size = sub.size();
    sub = sub.substr(7,temp_size-9);
    cout<<sub<<endl;

    getline(file,sub,'\n');
    temp_size = sub.size();
    sub = sub.substr(7,temp_size-9);
    cout<<sub<<endl;

    getline(file,sub,'\n');
    temp_size = sub.size();
    sub = sub.substr(7,temp_size-9);
    cout<<sub<<endl;

    getline(file,sub,'\n');
    temp_size = sub.size();
    sub = sub.substr(7,temp_size-9);
    cout<<sub<<endl;

    getline(file,sub,'\n');
    temp_size = sub.size();
    sub = sub.substr(10,temp_size-12);
    cout<<sub<<endl;

    getline(file,sub,'\n');
    temp_size = sub.size();
    sub = sub.substr(6,temp_size-8);
    cout<<sub<<endl;

    getline(file,sub,'\n');
    temp_size = sub.size();
    sub = sub.substr(6,temp_size-8);
    cout<<sub<<endl;
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
    temp_node = (IndexMatrix*)malloc(size*sizeof(IndexMatrix));

    for( count = 0;count<size;count++  )
    {
        array[count] = (int*)malloc(size*sizeof(int));
        temp_node[count].alphabet=str[count];
        temp_node[count].pos=count;
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
