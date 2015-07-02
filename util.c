/*************************************************************************
	> File Name: util.c
	> Author:watchlee 
	> Mail: windvergil@gmail.com
	> Created Time: Thu Jul  2 13:18:35 2015
 ************************************************************************/

#include<stdio.h>
#include <stdlib.h>

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
        puts(file_name);
    }
    return ;
}


