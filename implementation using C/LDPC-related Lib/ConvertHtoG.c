#include "ConvertHtoG.h"
#include "Htrsf.h"
#include <stdlib.h>
#include <stdio.h>

void ConvertHtoG(char *H, int n, int m, char **G, int *dim)
{
    char *full_H = (char *)calloc(n*m, sizeof(char));
    int *indpdt_cols = NULL;
    int *other_cols = NULL;
    int *inv_permutation = (int *)calloc(m, sizeof(int));
    int current_row = 0, sum = 0, k = 0;
    int flag = 0, num = 0, i = 0;
    int row = 0, col = 0, Rank_H = 0;
    char *Gt = NULL;


    // firstly, removing all the all-zero rows in the parity-check matrix H to get full-rank full_H
    for (row = 0; row < n; row++)
    {
        sum = 0;
        for (col = 0; col < m; col++)
            sum += *(H+(row*m+col));
        if (sum != 0)
        {
            for (col = 0; col < m; col++)
            {
                *(full_H + (current_row*m + col)) = *(H + (row*m + col));
            }
            current_row ++;
        }
    }

    Rank_H = current_row;  //the Rank of H is actually equivalent to the number of rows in full_H;
    
    indpdt_cols = (int *)calloc(Rank_H, sizeof(int));

    for (num=0; num < Rank_H; num++)
        if (*(full_H+(num*m+num)) == 1)
        {
  //          printf("num=========%d\n", num);
  //              getchar();
            *(indpdt_cols + num) = num;
            *(inv_permutation + num) = num + m - Rank_H;

        }
        else
        {
            for (col=0; col < m; col++)
            {
                flag = 1;
                for (i = 0; i < num; i++)
                    if (col == *(indpdt_cols + i))
                    {
                        flag = 0;
                        break;
                    }
                if (flag == 0)
                    continue;

                if  (*(full_H + (num*m+col)) == 1)
                {
                    *(indpdt_cols + num) = col;
                    *(inv_permutation + col) = num + m - Rank_H;
                    break;
                }
            }
        }

    num = 0;
    other_cols = (int *)calloc(m - Rank_H, sizeof(int));

    for (col = 0; col <m; col++)
    {
        flag = 1;
        for (i = 0; i < Rank_H; i++)
            if ( col == *(indpdt_cols + i))
            {
                flag = 0;
                break;
            }
        if (flag)
        {
            *(inv_permutation + col) = num;
            *(other_cols + num) = col;
            num ++;
        }

    }


//    for (i = 0; i < Rank_H; i++)
//        printf("i====%d, indpdt_col = %d\n", i, *(indpdt_cols + i));


   
    k = m - Rank_H;

    Gt =  (char *)calloc(k*m, sizeof(char));


    // Construct the permutated generator matrix Gt
    for (row = 0; row < k; row++)
    {
        *(Gt + (row*m + row)) = 1;
        for (num = 0; num < Rank_H; num++)
            *(Gt + (row*m + k + num)) = *(full_H + (num*m + *(other_cols + row)));
    }

    *G =  (char *)calloc(k*m, sizeof(char));
      
    for (col = 0; col < m; col++)
        for (row = 0; row < k; row++)
        {
            *((*G) + (row * m + col)) = *(Gt + (row * m + (*(inv_permutation + col))));
     //       printf("Gt----------%d\n", row * m + *(inv_permutation + col));
        }
   
   
    *dim = k;


    free(indpdt_cols);
    free(other_cols);
    free(inv_permutation);
    free(full_H);
    free(Gt);
    return;
}





 