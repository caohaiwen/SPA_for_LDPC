#include "SPA_BEC.h"
#include <stdlib.h>
#include <stdio.h>

char SPA_BEC(char *y_r, int n, int m, int row_w, int col_w, int *variable, int *check, int max_iterations, char **decoded_x, int *P)
{
    char *llr_rl = (char *)calloc(n*row_w, sizeof(char)) ;
    char *llr_lr = (char *)calloc(n*col_w, sizeof(char)) ;
    int degree_variable = col_w, degree_check = row_w;
    int its = 0, col = 0, row = 0, i = 0,checkN = 0, variableN = 0, t = 0, next_checkN = 0, next_variableN = 0,  cnt = 0, prv_cnt = 0;
    int edge = 0, next_edge = 0;

/* iteration0 : from left to right is just the message from the channel. */
    for (col = 0; col < m; col++)
        for (i = 0; i < degree_variable; i++)
        {
            checkN = *(variable + (col*degree_variable + i)) - 1;
            if (checkN == -1)
                continue;
            edge = *(P + (m*checkN + col));
            *(llr_lr + edge) = *(y_r + col);
            *((*decoded_x) + col) = *(y_r + col);
        }


    for (its = 0; its < max_iterations; its++)
    {

    /* messages pass from right(check nodes) to left(variable nodes) */
        for (row = 0; row < n; row++)   //enumerating the check nodes one by one

            for (i = 0; i < degree_check; i++)
            {
                variableN = *(check + (row*degree_check + i)) - 1;
                if (variableN == -1)
                    continue;
                edge = *(P + (row*m +variableN));
            /*
                the outgoing message from the check node is an erasure if there exists any erasure from the incoming
                message, otherwise, all the incoming messages are not erased and then mod (all the sum of incoming messages)
            */
                *(llr_rl + edge) = 0;
                for (t = i+1; t < i + degree_check; t++)
                {
                    next_variableN = *(check + (row*degree_check + t % degree_check)) - 1;
                    if (next_variableN == -1)
                        continue;
                    next_edge = *(P + (row*m + next_variableN));

                    if (*(llr_lr + next_edge) != 2)
                    {
                        *(llr_rl + edge) += (*(llr_lr + next_edge));
                        *(llr_rl + edge) %= 2;
                    }
                    else
                    {
                        *(llr_rl + edge) = 2;
                        break;
                    }

                }
            }
     /* messages pass from left(variable nodes) to right(check nodes) */
        for (col = 0; col < m; col++)   //enumerating the variables one by one
        {
            for (i = 0; i < degree_variable; i++)
            {
                checkN = *(variable + (col*degree_variable + i)) - 1;    //the i-th check node of the col-th variable node
                if (checkN == -1)
                    continue;
                edge = *(P + (m*checkN + col));

                if (*((*decoded_x) + col) != 2)
                    *(llr_lr + edge) = *((*decoded_x) + col) ;
                else
                {
                    /*
                        the outgoing message from the variable node is an erasure if all incoming messages are erasures;
                    */
                    for (t = i+1; t < i + degree_variable ; t++ )
                    {
                        next_checkN = *(variable + (col*degree_variable + t % degree_variable)) - 1;
                        if (next_checkN == -1)
                            continue;
                        next_edge = *(P + (m*next_checkN + col));

                        if (*(llr_rl + next_edge) != 2)
                        {
                            *(llr_lr + edge) = *(llr_rl + next_edge);
                            *((*decoded_x) + col) = *(llr_rl + next_edge);
                            break;
                        }
                    }
                }
            }
            if (*((*decoded_x) + col) != 2)
                cnt ++;
        }

        if (cnt == m)
        {
            free(llr_rl);
            free(llr_lr);
            return 1;
        }
        else if (cnt == prv_cnt)
            break;

        prv_cnt = cnt;
        cnt = 0;
    }

    free(llr_rl);
    free(llr_lr);
    return 0;
}
