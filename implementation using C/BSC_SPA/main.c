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
    double cross_prob = 0, p = 0;
    char *u = NULL;
    char *x_s = NULL, *y_r = NULL;
    double *cm_int = NULL;
    int sum = 0, num = 0, failure = 0, trial = 0;
    char *decoded_x = NULL;
    char convergence = 1;
    int total_trial = 0, max_iterations = 100, error_bits = 0, hd = 0;
    time_t t;

    fp = fopen("PCMatrix(816.3.174 (N=816,K=408,M=408,R=0.5)).txt", "r");

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

    free(H);

    u = (char *)calloc(k, sizeof(char));
    x_s = (char *)calloc(m, sizeof(char));
    y_r = (char *)calloc(m, sizeof(char));
    cm_int = (double *)calloc(m, sizeof(double));

    total_trial = 10000;
    num = 0;

    decoded_x = (char *)calloc(m, sizeof(char));

    WBER = fopen("WBER under BSC(N=816).txt", "w");

    if (WBER == NULL)
    {
        printf("cannot open this file to record the word error rate and bit error rate-_-\n");
        return 0;
    }

    for (i = 1; i <= 25; i++)
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
                p = ((double)(rand())) / RAND_MAX;
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
        }

        fprintf(WBER, "%lf, %lf\n",  cross_prob, (double)failure / total_trial);
        fprintf(WBER, "%lf, %lf\n",  cross_prob, (double)error_bits / (total_trial*m));
    }
    fclose(WBER);

    free(cm_int);
    free(u);
    free(decoded_x);
    free(x_s);
    free(y_r);
    free(G);
    free(variable);
    free(check);

    return 0 ;

}
