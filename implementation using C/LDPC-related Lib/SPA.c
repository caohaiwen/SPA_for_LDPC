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
        P: the matrix whose entries are the labeled edge numbers from 1 to #total edges or 0 for which position has no edge.
output : whether it can convergent to a codeword within the max_iterations.
*/


/*
    int abs(int i) ------> input: integer; output: integer
    double fabs(double x) ------> input float ; output: float
*/

char SPA(double *cm_int, int n, int m, int row_w, int col_w, int *variable, int *check, int max_iterations, char **decoded_x, int *P)
{
    double *llr_rl = (double *)calloc(n*row_w, sizeof(double)) ;
    double *llr_lr = (double *)calloc(m*col_w, sizeof(double)) ;
    int degree_variable = col_w, degree_check = row_w;
    double *beliefs = (double *)calloc(m, sizeof(double)) ;
    char flag = 1;
    int its = 0, col = 0, row = 0, i = 0,checkN = 0, variableN = 0, sum = 0;
    double tanhVaule = 1, value = 0;
    int edge = 0;
    double *open_variableN = (double *)calloc(2*(degree_variable - 1), sizeof(double));
    double *open_checkN = (double *)calloc(2*(degree_check - 1), sizeof(double));
    int num = 0;

    for (its = 0; its < max_iterations; its++)
    {
        flag  = 1;

    /* messages pass from left(variable nodes) to right(check nodes) */


        for (col = 0; col < m; col++)            //enumerating the variables one by one
        {
            /* opening the box, that is, we split the Degree-degree_variable variable node into degree_variable Degree-2 nodes */

            /* compute all the messages (after opening the box) from top to down */
            for (i = 0; i < degree_variable -1; i++)
            {
                checkN = *(variable + (col*degree_variable + i)) - 1;    //the i-th check node of the col-th variable node
                edge = *(P + (m*checkN + col));
                if (i == 0)
                    *(open_variableN + i) = (*(cm_int + col)) + (*(llr_rl + edge));  //messages from top to down
                else
                    *(open_variableN + i) = *(open_variableN + i -1) + *(llr_rl + edge);

            }

        /* compute all the messages (after opening the box) from down to top, and compute the outgoing messages */
            num = degree_variable - 1;
            for (i = degree_variable - 1; i > 0; i--)
            {
                checkN = *(variable + (col*degree_variable + i)) - 1;    //the i-th check node of the col-th variable node
                edge = *(P + (m*checkN + col));
                /*messages from down to top, and merge the upper and lower messages to compute the outgoing messages. */
                if (i == degree_variable -1)          //messages from down to top;
                {
                    *(open_variableN + num) = *(llr_rl + edge);
                    *(llr_lr + edge) = *(open_variableN + degree_variable - 2);
                }
                else
                {
                    *(open_variableN + num) = *(open_variableN + num - 1) + *(llr_rl + edge);
                    *(llr_lr +edge) = *(open_variableN + num - 1) + *(open_variableN + 2*(degree_variable - 2) - (num -1));
                }

                num ++;
            }

            /* compute the first outgoing message */
            checkN = *(variable + (col*degree_variable)) - 1;
            edge = *(P + (m*checkN + col));
            *(llr_lr + edge) = *(open_variableN + num - 1) + *(cm_int + col); // the intrinsic information from the channel

        }


    /* messages pass from right(check nodes) to left(variable nodes) */

        for (row = 0; row < n; row++)   //enumerating the check nodes one by one
        {

            /* open the box : messages passing from top to down */
            for (i = 0; i < degree_check - 1; i++)
            {
                variableN = *(check + (row*degree_check + i)) - 1;
            //adapted to the almost regular LDPC, that is, there exist some 0s indicating this check node has less than row_w degrees.
                if (variableN == -1)
                    continue;
                edge = *(P + (row*m +variableN));

                value = *(llr_lr +edge) / 2.0;
                tanhVaule = 1.0;
                if (fabs(value) < 17.5)               //otherwise, tanh(tanhVaule) = +1/-1
                    tanhVaule = tanh(value);
                else if (value < -17.5)
                    tanhVaule = -1.0;

                if (i == 0)
                    *(open_checkN + i) = tanhVaule;
                else
                    *(open_checkN + i) = *(open_checkN + i -1) * tanhVaule ;

            }

            /*open the box : messages passing from down to top , and also compute all the outgoing messages*/
            num = degree_check - 1;
            for (i = degree_check - 1; i > 0; i--)
            {
                variableN = *(check + (row*degree_check + i)) - 1;
                if (variableN == -1)
                    continue;
                edge = *(P + (row*m +variableN));

                value = *(llr_lr +edge) / 2.0;
                tanhVaule = 1.0;
                if (fabs(value) < 17.5)               //otherwise, tanh(tanhVaule) = +1/-1
                    tanhVaule = tanh(value);
                else if (value < -17.5)
                    tanhVaule = -1;
                /* messages from down to top, and merge the upper and lower messages to compute the outgoing messages. */
                if (i == degree_check -1)          //messages from down to top;
                {
                    *(open_checkN + num) = tanhVaule;
                    tanhVaule = *(open_checkN + degree_check - 2);
                }
                else
                {
                    *(open_checkN + num) = *(open_checkN + num - 1) * tanhVaule;
                    tanhVaule = *(open_checkN + num - 1) * (*(open_checkN + 2*(degree_check - 2) - (num -1)));
                }

                if (fabs(tanhVaule - 1.0) <= 1e-15)
                    value = 17.5;
                else if (fabs(tanhVaule - (-1.0)) <= 1e-15)
                    value = -17.5;
                else
                    value = atanh(tanhVaule);

                *(llr_rl + edge) = 2 * value;
                num ++;
            }

             /* compute the first outgoing message */
            variableN = *(check + (row * degree_check)) - 1;
            if (variableN == -1)
                continue;
            edge = *(P + (row*m +variableN));

            tanhVaule = *(open_checkN + num - 1);
            if (fabs(tanhVaule - 1.0) <= 1e-15)
                value = 17.5;
            else if (fabs(tanhVaule - (-1.0)) <= 1e-15)
                value = -17.5;
            else
                value = atanh(tanhVaule);

            *(llr_rl + edge) =  2 * value;
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
                edge = *(P + (checkN*m + col));
                *(beliefs + col) += *(llr_rl + edge);
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
            free(open_variableN);
            free(open_checkN);
            return 1;
        }
    }
    free(llr_rl);
    free(llr_lr);
    free(beliefs);
    free(open_variableN);
    free(open_checkN);
    return 0;
}
