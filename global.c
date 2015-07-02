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

/************************************************************************/

/*-------------------------Setting function-------------------------------*/
void input_function(char*);
void print_result(char*,int);
int max_compare_function(int,int);
void readmat(char*);
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
            }
            loop++;
        }
    }
#if DEBUG==TEST_PRINT
    printf("seq1 = %s\n",seq1);
    printf("seq2 = %s\n",seq2);
    printf("matfile = %s\n",matfile);
    printf("opp = %d\n",gap_opp);
    printf("exp = %d\n",gap_exp);
    printf("suboptimal = %d\n",gap_suboptimal);
#endif
   readmat(matfile); 


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

void readmat(char *file_name)
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
    

    printf("%s\n",line);
    fclose(fptr);
}

