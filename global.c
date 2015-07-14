/*************************************************************************
	> File Name: global.c
	> Author: LeePoHan
	> Mail: windvergil@gmail.com
	> Created Time: Wed Jul  1 13:34:18 2015
 ************************************************************************/

#include<stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>



/*
 * function: semiglobal alignment
 * Input : Working Path
 * Output: result.php
 *
*/

/*
 * Input:
 * $seq1 $seq2: sequence with structural alphabet
 * $matfile: Substitution matrix file
 * $opp: open gap penelty
 * $exp: extended gap penalty
 * $suboptimal: # of output solution
*/

/*-------------------------DEBUG setting-------------------------------*/
#define DEBUG 1


/************************************************************************/
/*-------------------------Setting Variable-------------------------------*/
static int input_length;
FILE *input_file;
char* seq1;
char* seq2;
char* matfile;
int gap_opp;
int gap_exp;
int gap_suboptimal;
int **scoring_matrix;
/*-------------------------Setting Array----------------------------------*/
int **A;
int **B;
int **C;
int **score;
char **dir;
char **dirA;
char **dirB;
char **dirC;
char *prefix;
char *column_array;
char *row_array;
/***************************************************************************/

/*-------------------------Setting function-------------------------------*/
void input_function(char*);
void print_result(char*,int);
void initialize_setting();
int max_compare_function(int,int);
int** readmat(char*);
void setting_scoring_function();
int file_exists(const char*);
bool array_key_exits();
/************************************************************************/
int main(int argc,char* argv[])
{
    /*-----------Here we must add prefix process function----------------2015/7/8--*/
    prefix = argv[1];
    if(strcmp(&prefix[strlen(prefix)-1],"/")!=0)
        strcat(prefix,"/");

    char *file = (char*)malloc(strlen(prefix)*sizeof(char));
    strcpy(file,prefix);
    strcat(file,"input.php");

    input_function(file);
    fclose(input_file);
    
    return 0;
}

void input_function(char *path)
{
    input_file = fopen(path,"r");
    if(input_file==NULL)
    {
        printf("where is input file?\n");
        exit(-1);
    }

    char *temp=(char*)malloc(100*sizeof(char));
//    char *temp;
    /*In any case, I decide to choice dynamic allocate memory,but not done yet! */
    int data_size = 0;
    int maximun_length = 0;
    int temporal_length = 0;
    while(!feof(input_file))
    {
        fscanf(input_file,"%s",temp);
        if(strncmp(temp,"=",1)!=0&&strncmp(temp,"$",1)!=0&&strncmp(temp,"?>",1)!=0&&strncmp(temp,"<?",1)!=0)
        {
            /*first input*/
            if(data_size==0)
            {
                maximun_length = strlen(temp);
            }
            else
            /*trying to find maximun length to allocate memory*/
            {
                temporal_length = strlen(temp);
                maximun_length = max_compare_function(maximun_length,temporal_length);
            }
            data_size++;
            
            
        }
    }
    fclose(input_file);

    input_file = fopen(path,"r");
    
    /*We have to get data from php */
    /*Input 
     *1 line = $seq1
     *2 line = $seq2
     *3 line = $matfile
     *4 line = $opp
     *5 line = $exp
     *6 line = $suboptimal
     */
    printf("Getting input data...\n");
    int loop  = 0;

    /*created temporal array to store data*/

/*
 *      In current statuation, I don't use it! 
 *      2015/7/2
 *
    char **temp_array,*temp_length;
    temp_array = (char**)malloc(data_size*sizeof(char*)+data_size*maximun_length*sizeof(char));
    for(loop = 0,temp_length=(char*)(temp_array+data_size);loop<data_size;loop++,temp_length+=maximun_length)
        temp_array[loop] = temp_length;
 *
 *       
*/
    loop = 0;
    
    /* It is very stupid way------2015/7/2*/
    while(!feof(input_file))
    {
        fscanf(input_file,"%s",temp);
        temp = strtok(temp,"\"");
        printf("%s\n",temp);
        if(strncmp(temp,"=",1)!=0&&strncmp(temp,"$",1)!=0&&strncmp(temp,"?>",1)!=0&&strncmp(temp,"<?",1)!=0)
        {
            /*When you copy data ,you must inialize your char pointer first, nor it will happen
             * segmentation fault 11 */
            
            switch(loop)
            {
                case 0:
                    seq1 = (char*)malloc(strlen(temp)*sizeof(char));
                    strcpy(seq1,temp);
                break;
                case 1:
                    seq2 = (char*)malloc(strlen(temp)*sizeof(char));
                    strcpy(seq2,temp);
                break;
                case 2:
                    matfile=(char*)malloc(strlen(temp)*sizeof(char));
                    strcpy(matfile,temp);
                break;
                case 3:
                    gap_opp = atoi(temp);
                break;
                case 4:
                    gap_exp = atoi(temp);
                break;
                case 5:
                    gap_suboptimal = atoi(temp);
                break;
                default:
                break;
            }
            loop++;
        }
    }

    /*Get scoring_matrix*/
    scoring_matrix= readmat(matfile);
#if DEBUG
    printf("seq1 = %s\n",seq1);
    printf("seq2 = %s\n",seq2);
    printf("matfile = %s\n",matfile);
    printf("opp = %d\n",gap_opp);
    printf("exp = %d\n",gap_exp);
    printf("suboptimal = %d\n",gap_suboptimal);
    int i,j;
    /*
    for(i =0;i<24;i++)
    {
        for(j = 0;j<24;j++)
            printf("%d ",scoring_matrix[i][j]);
        printf("\n");
    }
    */
#endif

    initialize_setting();
    int total_len = strlen(seq1)+strlen(seq2);
    char* temp_string = (char*)malloc(strlen(prefix)*sizeof(char));
    strcpy(temp_string,prefix);
    strcat(temp_string,"dirA_glo");
    /*------------------------Due to we don't need to consider memory problem, just skip this function
    if(total_len>4000&&file_exists(temp_string))
    {
        remove(temp_string);

        strcpy(temp_string,prefix);
        strcat(temp_string,"dirB_glo");
        remove(temp_string);

        strcpy(temp_string,prefix);
        strcat(temp_string,"dirC_glo");
        remove(temp_string);

        strcpy(temp_string,prefix);
        strcat(temp_string,"dir_glo");
        remove(temp_string);

        strcpy(temp_string,prefix);
        strcat(temp_string,"score_glo");
        remove(temp_string);
    }
    */
    
    /*1.找到Seq的X Y位置
     *2.利用此XY找到Matrix的對應位置
     *Sol:利用seq1 X and Seq2 Y取得的字元到matrix找 看是否能找到相對應的XY並回傳int的X Y
     */
     int inner_loop;
     for(loop = 0;loop<strlen(seq1);loop++)
     {
        for(inner_loop = 0;inner_loop<strlen(seq2);inner_loop++)
        {

        }
     }
    
/*------------------free memory-----------------*/ 
    free(scoring_matrix);
    free(A);
    free(B);
    free(C);
    free(dir);
    free(dirA);
    free(dirB);
    free(dirC);
    free(seq1);
    free(seq2);
    free(matfile);
    free(temp_string);
    free(row_array);
    free(column_array);
}

int max_compare_function(int first_data,int second_data)
{
    if(first_data>second_data)
        return first_data;
    else
        return second_data;
}

void print_result(char* array,int length)
{
    int loop;
    printf("print result\n");
    for(loop =0;loop<length;loop++)
        printf("%s\n",&array[loop]);
    
}

int** readmat(char *file_name)
{

    
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

    char line[100];
    /*first reading data is comment, so lets ingnore it!*/
    fgets(line,100,fptr);
    /*Start reading data from scoring matrix*/
    fgets(line,100,fptr);
    row_array = (char*)malloc(24*sizeof(char));
    column_array = (char*)malloc(24*sizeof(char));
    
    int loop = 1;
/*---------------Testing------------------------*/
    /*
     *
     *First at all, I used test_index to caculate how many data that I should store into array. And also, I used very unsmartly way to store alphabet characters. 
     *Comment Time : 2015/7/10
     *
     *
     */
    char *temp = strtok(line," ");
    row_array[0] = *temp;
    column_array[0]=*temp;
    
    int test_index=1;
    do
    {
        
        temp = strtok(NULL," ");
        if(temp!=NULL)
        {
            column_array[loop]=*temp;
            row_array[loop++]=*temp;
            test_index++;
        }
    }while(temp!=NULL);

#if DEBUG
    printf("Row array = \n");
    for(loop = 0;loop<24;loop++)
        printf("%c ",row_array[loop]);
    printf("\n");
    printf("Column array = \n");
    for(loop = 0;loop<24;loop++)
        printf("%c ",column_array[loop]);
    printf("\ntest_index = %d\n\n",test_index);
#endif

/*----------------------------------------------*/
    
    const int size = test_index;
    int **array = (int**)malloc(size*sizeof(void *));
    int count;
    for( count = 0;count<size;count++  )
    {
        array[count] = (int*)malloc(size*sizeof(int));
    }
    /*Array index value*/
    int array_row ,array_col = -1;
    loop =0;

    

    while(!feof(fptr))
    {
        array_row = 0; 
        char *temp = strtok(line," ");
        printf("%s ",temp);
        temp = strtok(NULL," ");
        while(temp!=NULL)
        {
            array[array_col][array_row] = atoi(temp);
            printf("%d ",atoi(temp));
            array_row++;
            temp = strtok(NULL," ");
        }
        
        array_col++;
        fgets(line,100,fptr);
        printf("\n");
    }
    

    fclose(fptr);
    printf("Done!\n");
/*----------------Testing-------------------*/
#if DEBUG

/*
 *Should add output_test_file 2015/7/7
 *
 *
 */
    int i,j;
    for(i =0;i<24;i++)
    {
        printf("i = %d    ",i);
        for(j = 0;j<24;j++)
        {
//            printf("%d ",data_array[i][j]);
            printf("%d ",array[i][j]);
            
        }
        printf("\n");
    }
    fptr = fopen(file_name,"r");
    fgets(line,100,fptr);
    /*Start reading data from scoring matrix*/
    fgets(line,100,fptr);
    fgets(line,100,fptr); 
    i = 0;
    while(!feof(fptr))
    {
        j = 0;
        char *temp = strtok(line," ");
        temp = strtok(NULL," ");
        
        while(temp!=NULL)
        {
            if(atoi(temp)==array[i][j])
            {
                temp = strtok(NULL," ");
                j++;
                
            }
            else
            {
                printf("比對失敗!!!!!!!!!!\n");
                printf("%d and array[%d][%d] = %d was compared failed!",atoi(temp),i,j,array[i][j]);
                exit(EXIT_FAILURE);
            }
            
        }
        if(j==24&&i==23)
        {
            printf("完全符合!\t\n");
        }
        i++;
        fgets(line,100,fptr);
    }
    fclose(fptr);
#endif
    return array;
}

/*-----------------------initializing all setting variable--------------------------*/
void initialize_setting()
{
    int i;
    A = (int**)malloc(strlen(seq1)*sizeof(int*)+strlen(seq1)*strlen(seq2)*sizeof(int));
    B = (int**)malloc(strlen(seq1)*sizeof(int*)+strlen(seq1)*strlen(seq2)*sizeof(int));
    C = (int**)malloc(strlen(seq1)*sizeof(int*)+strlen(seq1)*strlen(seq2)*sizeof(int));
    score=(int**)malloc(strlen(seq1)*sizeof(int*)+strlen(seq1)*strlen(seq2)*sizeof(int));
    dir = (char**)malloc(strlen(seq1)*sizeof(char*)+strlen(seq1)*strlen(seq2)*sizeof(char));
    dirA = (char**)malloc(strlen(seq1)*sizeof(char*)+strlen(seq1)*strlen(seq2)*sizeof(char));
    dirB = (char**)malloc(strlen(seq1)*sizeof(char*)+strlen(seq1)*strlen(seq2)*sizeof(char));
    dirC = (char**)malloc(strlen(seq1)*sizeof(char*)+strlen(seq1)*strlen(seq2)*sizeof(char));

    int *temp;
    for(i = 0,temp = (int*)(A+strlen(seq1));i<strlen(seq1);temp+=strlen(seq2),i++)
        A[i] = temp;

    int *tempB;
    for(i = 0,tempB = (int*)(B+strlen(seq1));i<strlen(seq1);tempB+=strlen(seq2),i++)
        B[i] = tempB;

    int *tempC;
    for(i = 0,tempC = (int*)(C+strlen(seq1));i<strlen(seq1);tempC+=strlen(seq2),i++)
        C[i] = tempC;

    int *tempScore;
    for(i =0,tempScore = (int*)(score+strlen(seq1));i<strlen(seq1);i++,tempScore+=strlen(seq2))
        score[i] = tempScore;

    char* temp_dir;
    for(i = 0,temp_dir = (char*)(dir+strlen(seq1));i<strlen(seq1);temp_dir+=strlen(seq2),i++)
        dir[i] = temp_dir;

    char* temp_dirA;
    for(i = 0,temp_dirA = (char*)(dirA+strlen(seq1));i<strlen(seq1);temp_dirA+=strlen(seq2),i++)
        dirA[i]= temp_dirA;

    char* temp_dirB;
    for(i = 0,temp_dirB=(char*)(dirC+strlen(seq1));i<strlen(seq1);temp_dirB+=strlen(seq2),i++)
        dirB[i]=temp_dirB;

    char* temp_dirC;
    for(i = 0,temp_dirC = (char*)(dirC+strlen(seq1));i<strlen(seq1);i++,temp_dirC+=strlen(seq2))
        dirC[i] = temp_dirC;

    A[0][0] = B[0][0] = C[0][0]= score[0][0] = 0;
    dir[0][0] = 'o';  

/*-----------------do first "-" (row)-----------------*/
    for(i =1;i<strlen(seq2);i++)
    {
        A[0][i] = -1000000;
        B[0][i] = -1000000; 
        C[0][i] = gap_opp+gap_exp*i;
        score[0][i] = gap_opp+gap_exp*i;
        dir[0][i] = 'C';
        dirA[0][i] = 'C';
        dirB[0][i] = 'C';
        dirC[0][i] = 'C';
    }
    
/*-----------------do first "|" (column)-----------------*/
    for(i =1;i<strlen(seq1);i++)
    {
        A[i][0] = -1000000;
        B[i][0] = gap_opp+gap_exp*i; 
        C[i][0] = -1000000; 
        score[i][0] = gap_opp+gap_exp*i;
        dir[i][0] = 'B';
        dirA[i][0] = 'B';
        dirB[i][0] = 'B';
        dirC[i][0] = 'B';
    }
#ifdef DEBUG
    //printf("%d %d %d %d %c\n",A[0][0],B[0][0],C[0][0],score[0][0],dir[0][0]);
    FILE *debug_file = fopen("test_result.txt","w");
    int index_row,index_column;
    fprintf(debug_file,"A array is \n");
    for(index_column = 0;index_column<strlen(seq1);index_column++)
    {
        for(index_row = 0;index_row<strlen(seq2);index_row++)
        {
            fprintf(debug_file,"%d ",A[index_column][index_row]);
        }
        fprintf(debug_file,"\n");
    }
    fprintf(debug_file,"B array is \n");
    for(index_column = 0;index_column<strlen(seq1);index_column++)
    {
        for(index_row = 0;index_row<strlen(seq2);index_row++)
        {
            fprintf(debug_file,"%d ",B[index_column][index_row]);
        }
        fprintf(debug_file,"\n");
    }
        fprintf(debug_file,"C array is \n");
    for(index_column = 0;index_column<strlen(seq1);index_column++)
    {
        for(index_row = 0;index_row<strlen(seq2);index_row++)
        {
            fprintf(debug_file,"%d ",C[index_column][index_row]);
        }
        fprintf(debug_file,"\n");
    }
        fprintf(debug_file,"score array is \n");
    for(index_column = 0;index_column<strlen(seq1);index_column++)
    {
        for(index_row = 0;index_row<strlen(seq2);index_row++)
        {
            fprintf(debug_file,"%d ",score[index_column][index_row]);
        }
        fprintf(debug_file,"\n");
    }
        fprintf(debug_file,"dir array is \n");
    for(index_column = 0;index_column<strlen(seq1);index_column++)
    {
        for(index_row = 0;index_row<strlen(seq2);index_row++)
        {
            fprintf(debug_file,"%d ",dir[index_column][index_row]);
        }
        fprintf(debug_file,"\n");
    }
        fprintf(debug_file,"dirA array is \n");
    for(index_column = 0;index_column<strlen(seq1);index_column++)
    {
        for(index_row = 0;index_row<strlen(seq2);index_row++)
        {
            fprintf(debug_file,"%c ",dirA[index_column][index_row]);
        }
        fprintf(debug_file,"\n");
    }
        fprintf(debug_file,"dirB array is \n");
    for(index_column = 0;index_column<strlen(seq1);index_column++)
    {
        for(index_row = 0;index_row<strlen(seq2);index_row++)
        {
            fprintf(debug_file,"%c ",dirB[index_column][index_row]);
        }
        fprintf(debug_file,"\n");
    }
        fprintf(debug_file,"dirC array is \n");
    for(index_column = 0;index_column<strlen(seq1);index_column++)
    {
        for(index_row = 0;index_row<strlen(seq2);index_row++)
        {
            fprintf(debug_file,"%c ",dirC[index_column][index_row]);
        }
        fprintf(debug_file,"\n");
    }
    fclose(debug_file);
#endif

}
/*-------------------Determine whether file is existed or not?--------------------*/
int file_exists(const char* path)
{
    FILE *fptr = fopen(path,"r");
    if(fptr!=NULL)
    {
        fclose(fptr);
        return 1;
    }
    else
        return 0;
        
}
