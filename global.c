/*************************************************************************
	> File Name: global.c
	> Author: LeePoHan
	> Mail: windvergil@gmail.com
	> Created Time: Wed Jul  1 13:34:18 2015
 ************************************************************************/

#include<stdio.h>
#include <stdlib.h>
#include <string.h>



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
enum{TEST_PRINT,TEST_ERROR,TEST_FUNCTION};
#define DEBUG TEST_PRINT


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
/************************************************************************/

/*-------------------------Setting function-------------------------------*/
void input_function(char*);
void print_result(char*,int);
int max_compare_function(int,int);
int** readmat(char*);
void setting_scoring_function();
/************************************************************************/
int main(int argc,char* argv[])
{
    char *file = argv[1];
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
    char *temp;

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
 *      current state, I don't use it! 
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
            /*
            if(loop ==0)
            {
                seq1=(char*)malloc(strlen(temp)*sizeof(char));
                strcpy(seq1,temp);
            }
            else if(loop==1)
            {

                seq2=(char*)malloc(strlen(temp)*sizeof(char));
                strcpy(seq2,temp);
            }
            else if(loop==2)
            {
                     
                matfile=(char*)malloc(strlen(temp)*sizeof(char));
                strcpy(matfile,temp);
            }
            else if(loop==3)
            {
                gap_opp = atoi(temp);
            }
            else if(loop==4)
            {

                gap_exp = atoi(temp);
            }
            else
            {
                
                gap_suboptimal= atoi(temp);
            }*/
            loop++;
        }
    }

    scoring_matrix= readmat(matfile);
#if DEBUG==TEST_PRINT
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


    free(scoring_matrix);

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


/*---------------Testing------------------------*/
    char *temp = strtok(line," ");

    int test_index=1;
    do
    {
        
        temp = strtok(NULL," ");
        if(temp!=NULL)
            test_index++;
    }while(temp!=NULL);

    printf("test_index = %d\n",test_index);
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
    
    while(!feof(fptr))
    {
        array_row = 0; 
        char *temp = strtok(line," ");
        temp = strtok(NULL," ");
        while(temp!=NULL)
        {
//            data_array[array_col][array_row]= atoi(temp);
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
#if DEBUG==TEST_PRINT

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

