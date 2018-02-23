#include "SPA.h"
#include <stdlib.h>
#include <math.h>

/* SPA implementation
input : cm_int: the intrinsic information from the channel model;
        n,m: the size of parity-check matrix, i.e., #check nodes and #variable nodes;
        col_w: degree of variable nodes in the tanner graph
        row_w: degree of check nodes in the tanner graph
        max_iterations : the maximum iterations during the iterative decoding in the tanner graph.
        decoded_x: the output of the decoder, i.e. , estimate of the received signal.
output : whether it can convergent to a codeword within the max_iterations.
*/


/*
    int abs(int i) ------> input: integer; output: integer
    double fabs(double x) ------> input float ; output: float
*/

char SPA(double *cm_int, int n, int m, int row_w, int col_w, int *variable, int *check, int max_iterations, char **decoded_x)
{
    double *llr_rl = (double *)calloc(n*m, sizeof(double)) ;
    double *llr_lr = (double *)calloc(m*n, sizeof(double)) ;
    int degree_variable = col_w, degree_check = row_w;
    double *beliefs = (double *)calloc(m, sizeof(double)) ;
    char flag = 1;
    int its = 0, col = 0, row = 0, i = 0,checkN = 0, variableN = 0, t = 0, next_checkN = 0, next_variableN = 0, sum = 0;
    double tanhVaule = 1, value = 0;

    for (its = 0; its < max_iterations; its++)
    {
        flag  = 1;

    /* messages pass from left(variable nodes) to right(check nodes) */
        for (col = 0; col < m; col++)            //enumerating the variables one by one
            for (i = 0; i < degree_variable; i++)
            {
                checkN = *(variable + (col*degree_variable + i)) - 1;    //the i-th check node of the col-th variable node

            /* compute the message from variable node col to its i-th check node. */
                *(llr_lr + (m*checkN + col)) = *(cm_int + col) ;

            //    *(llr_lr + (col*degree_variable + i)) = *(cm_int + col) ;

                for (t = i+1; t < i + degree_variable ; t++ )
                {
                    next_checkN = *(variable + (col*degree_variable + t % degree_variable)) - 1;
                    *(llr_lr + (m*checkN + col)) += *(llr_rl + (m*next_checkN + col));

                }


            }


    /* messages pass from right(check nodes) to left(variable nodes) */
        for (row = 0; row < n; row++)   //enumerating the check nodes one by one
            for (i = 0; i < degree_check; i++)
            {
                variableN = *(check + (row*degree_check + i)) - 1;
                if (variableN == -1)
                    continue;

            /* compute the message from the check node row to its i-th variable node. */
                tanhVaule = 1.0;
                for (t = i+1; t < i + degree_check; t++)
                {
                    next_variableN = *(check + (row*degree_check + t % degree_check)) - 1;
                    if (next_variableN == -1)
                        continue;

                    value = *(llr_lr + (row*m + next_variableN)) / 2.0 ;

                    if (fabs(value) < 17.5)               //otherwise, tanh(tanhVaule) = +1/-1
                        tanhVaule *= tanh(value);
                    else if (value < -17.5)
                        tanhVaule = -tanhVaule;
                }

                if (fabs(tanhVaule - 1.0) <= 1e-15)
                    value = 17.5;
                else if (fabs(tanhVaule - (-1.0)) <= 1e-15)
                    value = -17.5;
                else
                    value = atanh(tanhVaule);
                *(llr_rl + (row*m +variableN)) = 2 * value;
            }
    /* compute the current beliefs in each variable node */
        for (col = 0; col < m; col++)
        {
            *(beliefs + col) = *(cm_int + col);
            for (i = 0; i < degree_variable; i++)
            {
                checkN = *(variable + (col*degree_variable + i)) - 1;
                if (checkN == -1)
                    continue;
                *(beliefs + col) += *(llr_rl + (checkN*m + col));
            }
        }

        for (col = 0; col < m ; col++)
            if ( *(beliefs + col) >= 0)
                *(*decoded_x + col) = 0;
            else
                *(*decoded_x + col) = 1;

    /* check whether this word is a real codeword or not */
        for (row = 0; row < n; row++)
        {
            sum = 0;
            for (i = 0; i < degree_check; i++)
            {
                variableN = *(check + (row*degree_check + i)) - 1;
                sum += *(*decoded_x + variableN);
            }

            if (sum % 2 == 1)
            {
                flag = 0;
                break;
            }
        }

        if (flag==1)
        {
            free(llr_rl);
            free(llr_lr);
            free(beliefs);
            return 1;
        }
    }
    free(llr_rl);
    free(llr_lr);
    free(beliefs);
    return 0;
}
