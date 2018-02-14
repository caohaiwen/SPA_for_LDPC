#include "Htrsf.h"
#include <stdlib.h>
#include <stdio.h>

void Htrsf(char **H, int n, int m)  //Assume that n <= m;
{
    int i = 0, j = 0, k = 0;
    int *zero_rows = (int *)calloc(n, sizeof(int));
    int num_zeros = 0, current_row = 0, num = 0;
    int col = 0, row = 0, coll = 0;
    char temp = 0, flag = 1;

    if (zero_rows == NULL)
    {
        printf("calloc error\n");
        return;
    }

    while ((i < n) && (j < m))
    {
        // firstly, find the index of largest element(i.e. one) in the remainder of column j.
        k = i;
        flag = 1;
        for (row = i; row < n ; row++)
            if (*((*H)+(row*m+j)) == 1)
            {
                flag = 0;
                k = row;
                break;
            }

        if (flag==1) // the remainders of column j are all zero, then continue to the j+1 column;
        {
            *(zero_rows + num_zeros) = i;
            num_zeros ++;
            i ++;
            j ++;
            continue;
        }

        // otherwise, swap the i-th row and k-th row
        for (col = j; col < m; col++)
        {
            temp = (*((*H) + (i*m + col)));
            *((*H) + (i*m + col)) = *((*H) + (k*m + col));
            *((*H) + (k*m + col)) = temp;
          //  printf("temppp ==== %d\n", temp);
        }

        //subtract the pivot row form all the other rows;
        for (row=0; row < n; row++)
            if ( (row != i) && (*((*H) + (row*m + j)) == 1))
            {
                for (col = j; col < m; col++)
                {
                    (*((*H) + (row*m + col))) += (*((*H) + (i*m + col)));
                    (*((*H) + (row*m + col))) %= 2;
                }
            }
        i++;
        j++;
    }

    for (num=0; num < num_zeros; num++)
    {
        current_row = (*(zero_rows + num));
        k = current_row;
        for (col = n; col < m; col++)
            if (*((*H)+(current_row*m+col)) == 1)
            {
                for (row = 0; row < n; row++)
                    if ((row != current_row) && ((*((*H)+(row*m + col))) == 1))
                    {
                        for (coll = col; coll < m; coll++)
                        {
                            *((*H)+(row*m + coll)) += (*((*H)+(current_row*m + coll)));
                            *((*H)+(row*m + coll)) %= 2;
                        }
                    }
                break;
            }

    }
  //  printf("H--------3\n");
    free(zero_rows);
    return;
}
