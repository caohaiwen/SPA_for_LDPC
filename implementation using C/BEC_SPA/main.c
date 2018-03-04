#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "SPA_BEC.h"
#include "ConvertHtoG.h"
#include "ReadOutH.h"
#include "Htrsf.h"
#include "twister.h"

int main(int argc, char *argv[])  //using windows command to pass parameters about input file name and output file name
{
    FILE *fp = NULL, *WBER = NULL;
    int n = 0, m = 0, row_w = 0, col_w = 0;
    int *variable = NULL, *check = NULL;
    int i = 0, j = 0 , l = 0;
    char *H = NULL, *G = NULL;
    int k = 0;
    double erasure_prob = 0, p = 0;
    char *u = NULL;
    char *x_s = NULL, *y_r = NULL;
    int sum = 0, num = 0, failure = 0, trial = 0;
    char *decoded_x = NULL;
    char convergence = 1;
    int total_trial = 0, max_iterations = 100, error_bits = 0, hd = 0;
    time_t t;
    unsigned long seed;
    int edge = 0;
    int *P = NULL;

    fp = fopen(argv[1], "r");

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
            *(P + i*m + (*(check + i*row_w + j)) - 1) = edge;
            edge ++;
        }

    Htrsf(&H, n, m);
    ConvertHtoG(H, n, m, &G, &k);

    free(H);

    u = (char *)calloc(k, sizeof(char));
    x_s = (char *)calloc(m, sizeof(char));
    y_r = (char *)calloc(m, sizeof(char));

    total_trial = 100000;
    num = 0;

    decoded_x = (char *)calloc(m, sizeof(char));

    WBER = fopen(argv[2], "w");

    if (WBER == NULL)
    {
        printf("cannot open this file to record the word error rate and bit error rate-_-\n");
            return 0;
    }

    seed = (unsigned) time(&t);
    seed += ((seed+1) % 2 ); // make sure that seed is odd

    for (num = 3; num <= 7; num++)
    {
        max_iterations = (int)pow(2, num);
        for (i = 1; i <= 26; i++)
        {
            erasure_prob = i * 0.02;
            failure = 0;
            error_bits = 0;
            seedMT(seed);

            for (trial = 0; trial < total_trial; trial++)
            {

                for (j = 0; j < k; j++) //randomly generating information bits
                    *(u + j) =  randomMT()  % 2;

                for (j = 0; j < m; j++) // encoding----->u*G & distorted by bsc
                {
                    sum = 0;
                    p = ((double)( randomMT() )) / Max_RandMT;
                    for (l = 0; l < k; l++)
                        sum += ((*(u + l)) * (*(G + (l*m + j))));
                    *(x_s + j) = sum % 2;

                    if ( p <= erasure_prob )
                        *(y_r + j) = 2;
                    else
                        *(y_r + j) = (*(x_s + j));
                }

                convergence = SPA_BEC(y_r, n, m, row_w, col_w, variable, check, max_iterations, &decoded_x, P);
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
            fprintf(WBER, "max_iterations = %d, erasure_prob == %lf, word error rate == %lf\n",  max_iterations, erasure_prob, (double)failure / total_trial);
            fprintf(WBER, "max_iterations = %d, erasure_prob == %lf, bit error rate == %lf\n",  max_iterations, erasure_prob, (double)error_bits / (total_trial*m));
        }
    }
    fclose(WBER);
    fclose(WBER);

    free(P);
    free(u);
    free(decoded_x);
    free(x_s);
    free(y_r);
    free(G);
    free(variable);
    free(check);

    return 0 ;

}
