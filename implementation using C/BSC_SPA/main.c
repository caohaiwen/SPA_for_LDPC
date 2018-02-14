#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "SPA.h"
#include "ConvertHtoG.h"
#include "ReadOutH.h"
#include "Htrsf.h"

int main()
{
    FILE *fp = NULL, *WBER = NULL;
    int n = 0, m = 0, row_w = 0, col_w = 0;
    int *variable = NULL, *check = NULL;
    int i = 0, j = 0 , l = 0;
    char *H = NULL, *G = NULL;
    int k = 0;
    float cross_prob = 0, p = 0;
    char *u = NULL;
    char *x_s = NULL, *y_r = NULL;
    float *cm_int = NULL;
    int sum = 0, num = 0, failure = 0, trial = 0;
    char *decoded_x = NULL;
    char convergence = 1;
    int total_trial = 0, max_iterations = 100, error_bits = 0, hd = 0;
    time_t t;

    fp = fopen("PCMatrix(4095.738.4.102 (N=4095,K=3357,M=738,R= 0.82)).txt", "r");

    if (fp == NULL)
    {
        printf("cannot open this file-_-\n");
        return 0;
    }

    ReadOutH(fp, &n, &m, &row_w, &col_w, &variable, &check);
    fclose(fp);

    H = (char *)calloc(n*m, sizeof(char));
    for (i = 0; i < n; i++)
        for (j = 0; j < row_w; j++)
        {
            if ((*(check + i*row_w + j)) == 0)
                continue;
            *(H + i*m + (*(check + i*row_w + j)) - 1) = 1;
        }

    Htrsf(&H, n, m);

    ConvertHtoG(H, n, m, &G, &k);

    u = (char *)calloc(k, sizeof(char));
    x_s = (char *)calloc(m, sizeof(char));
    y_r = (char *)calloc(m, sizeof(char));
    cm_int = (float *)calloc(m, sizeof(float));

    total_trial = 10000;
    num = 0;

    decoded_x = (char *)calloc(m, sizeof(char));

    WBER = fopen("WBER under BSC(N=4095).txt", "w");

    if (WBER == NULL)
    {
        printf("cannot open this file to record the word error rate and bit error rate-_-\n");
        return 0;
    }

    for (i = 1; i < 20; i++)
    {
        cross_prob = i*0.02;
        failure = 0;
        error_bits = 0;
        num ++;
        srand((unsigned) time(&t));

        for (trial = 0; trial < total_trial; trial++)
        {

            for (j = 0; j < k; j++) //randomly generating information bits
                *(u + j) = rand() % 2;

            for (j = 0; j < m; j++) // encoding----->u*G & distorted by bsc
            {
                sum = 0;
                p = ((float)(rand())) / RAND_MAX;
                for (l = 0; l < k; l++)
                    sum += ((*(u + l)) * (*(G + (l*m + j))));
                *(x_s + j) = sum % 2;

                if (p <= cross_prob)
                    *(y_r + j) = 1 - (*(x_s + j));
                else
                    *(y_r + j) = (*(x_s + j));

                if (*(y_r + j) == 1)
                    *(cm_int + j) = log(cross_prob / (1 - cross_prob));
                else
                    *(cm_int + j) = log((1 - cross_prob) / cross_prob);
            }

            convergence = SPA(cm_int, n, m, row_w, col_w, variable, check, max_iterations, &decoded_x);

//if the algorithm doesn't convergent to a codeword, then we assume that it's all-zero codewords;
            if (convergence == 0)
            {
                failure ++;
                hd = 0;
                for (j = 0; j < m; j++)
                    hd += (*(x_s + j));
                error_bits += hd;

            }
            else
            {
                hd = 0;
                for (j = 0; j < m; j++)
                    if ((*(x_s + j)) != (*(decoded_x + j)))
                        hd++;
                if (hd != 0)
                {
                    failure ++;
                    error_bits += hd;
                }
            }

//            memset(cm_int, 0, m*sizeof(float));
//            memset(decoded_x, 0, m*sizeof(char));
        }

        fprintf(WBER, "cross_prob == %f, word error rate == %f\n",  cross_prob, (float)failure / total_trial);
        fprintf(WBER, "cross_prob == %f, bit error rate == %f\n",  cross_prob, (float)error_bits / (total_trial*m));
        if (abs(failure - total_trial) <= total_trial*0.001)
            break;
    }
    fclose(WBER);

    free(cm_int);
    free(u);
    free(decoded_x);
    free(x_s);
    free(y_r);
    free(G);
    free(H);
    free(variable);
    free(check);

    return 0 ;

}
