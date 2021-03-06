#include <stdio.h>
#include <stdlib.h>

#include "ReadOutH.h"

void ReadOutH(FILE *fp, int *n, int *m, int *row_w, int *col_w, int **variable, int **check) //pointer pointing to another pointer
{
    char content[300];
    char temp[6];
    int num_line = 0;
    int i = 0, j = 0, num = 0, l = 0;

    while (!feof(fp))
    {
        num_line++ ;
        if ((fgets(content, 300, fp)) == NULL)
        {
            printf("there is an error when reading this line from this file\n");
            return;
        }
        i = 0;
        j = 0;
        num = 0;
        if (num_line ==1)
        {
            while (content[i] != '\0')
            {
                if (content[i] == ' ')
                {
                    *m = atoi(temp);
                    for (l = 0; l < 6; l++)
                    	temp[l] = '\0';
                    j = 0;
                }
                else
                {
                    temp[j] = content[i];
                    j++;
                }
                i++;
            }
            *n = atoi(temp);
            for (l = 0; l < 6; l++)
            	temp[l] = '\0';
            continue;
        }

        if (num_line ==2)
        {
            while (content[i] != '\0')
            {
                if (content[i] == ' ')
                {
                    *col_w = atoi(temp);
                    for (l = 0; l < 6; l++)
                    	temp[l] = '\0';
                    j = 0;
                }
                else
                {
                    temp[j] = content[i];
                    j++;
                }
                i++;
            }
            *row_w = atoi(temp);
            *variable = (int *)calloc((*m)*(*col_w), sizeof(int));
            *check = (int *)calloc((*n)*(*row_w), sizeof(int));
            for (l = 0; l < 6; l++)
            	temp[l] = '\0';
            continue;
        }

        if (num_line > 2 && num_line <= 2 + (*m))
        {
            while (content[i] != '\0')
            {
                if (isspace(content[i]))
                {
                    *((*variable) + ((num_line - 3)* (*col_w) + num)) = atoi(temp);
                    for (l = 0; l < 6; l++)
                    	temp[l] = '\0';
                    j = 0;
                    num++;
                }
                else
                {
                    temp[j] = content[i];
                    j++;
                }

                i++;
            }
   //         *(variable + ((num_line - 3)*(*col_w) + num)) = atoi(temp);
  //          for (l = 0; l < 6; l++)
    //        	temp[l] = '\0';
        }

        if (num_line > 2 + *m)
        {
            while (content[i] != '\0')
            {
                if (isspace(content[i]))
                {
                    *((*check) + ((num_line - 3 - *m)* (*row_w) + num)) = atoi(temp);
                    for (l = 0; l < 6; l++)
                    	temp[l] = '\0';
                    j = 0;
                    num++;
                }
                else
                {
                    temp[j] = content[i];
                    j++;
                }

                i++;
            }
//            *(check + ((num_line - 3 - *m)*(*row_w) + num)) = atoi(temp);
//            for (l = 0; l < 6; l++)
  //          	temp[l] = '\0';
        }

    }

    *((*check) + ((num_line - 3 - *m)* (*row_w) + num)) = atoi(temp); // there is no space sign in the last line 
/*
	for (i = 0; i < *m; i++)
    {
        for (j = 0; j < *col_w; j++)
            printf("%d ", *((*variable) + (i*(*col_w) + j)));
        printf("\n");
    }

printf("\n");

    
    for (i = 0; i < *n; i++)
    {
        for (j = 0; j < *row_w; j++)
            printf("%d ", *((*check) + (i*(*row_w) + j)));
        printf("\n");
    }
*/
    return;
}
