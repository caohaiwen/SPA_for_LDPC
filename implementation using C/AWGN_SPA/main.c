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
#include "Gaussian.h"
#include "twister.h"



// tanh(17.5) = 0.999999999999999 (%.15lf) under "double" whose precision is about 15~16 decimal digits.
// RAND_MAX =32767

int main()
{
    FILE *fp = NULL, *WBER = NULL;
    int n = 0, m = 0, row_w = 0, col_w = 0;
    int *variable = NULL, *check = NULL;
    int i = 0, j = 0 , l = 0;
    char *H = NULL, *G = NULL;
    int k = 0;
    double SNR_b = 0, noise_var = 0, noise = 0, R = 0;
    char *u = NULL;
    char *x_s = NULL;
    double *cm_int = NULL, *y_r = NULL, *x_ss = NULL;
    int sum = 0, num = 0, failure = 0, trial = 0;
    char *decoded_x = NULL;
    char convergence = 1;
    int *P = NULL;
    int total_trial = 0, max_iterations = 100, error_bits = 0, hd = 0;
    time_t t;
    int E_b = 1;    //the input energy per information symbol.---> E_c = R*E_b;
    unsigned long seed = 0;
    int edge = 0;

    fp = fopen("PCMatrix(96.3.963 (N=96,K=48,M=48,R=0.5)).txt", "r");

    if (fp == NULL)
    {
        printf("cannot open this file-_-\n");
        return 0;
    }
    ReadOutH(fp, &n, &m, &row_w, &col_w, &variable, &check);
    fclose(fp);

    H = (char *)calloc(n*m, sizeof(char));
    P = (int *)calloc(n*m, sizeof(int));

    for (i = 0; i < n; i++)
        for (j = 0; j < row_w; j++)
        {
            if ((*(check + i*row_w + j)) == 0)
                continue;
            *(H + i*m + (*(check + i*row_w + j)) - 1) = 1;
            *(P + i*m + (*(check + i*row_w + j)) - 1) = edge;
            edge++;
        }

    Htrsf(&H, n, m);

    ConvertHtoG(H, n, m, &G, &k);

    free(H);

    R = (double)k / m;

    u = (char *)calloc(k, sizeof(char));
    x_s = (char *)calloc(m, sizeof(char));
    y_r = (double *)calloc(m, sizeof(double));
    x_ss = (double *)calloc(m, sizeof(double));
    cm_int = (double *)calloc(m, sizeof(double));

    total_trial = 100000;
    num = 0;
    decoded_x = (char *)calloc(m, sizeof(char));

    WBER = fopen("WBER under AWGN Channel(N = 96).txt", "w");

    if (WBER == NULL)
    {
        printf("cannot open this file to record the word error rate and bit error rate-_-\n");
        return 0;
    }

    seed =(unsigned) time(&t);
    if (seed % 2==0)
        seed ++;


    for (i = 0; i < 20; i++)
    {
        SNR_b = i*0.4;
        failure = 0;
        error_bits = 0;
        num ++;
        seedMT(seed);
        for (trial = 0; trial < total_trial; trial++)
        {

            for (j = 0; j < k; j++) //randomly generating information bits
                *(u + j) = randomMT() % 2;

            for (j = 0; j < m; j++) // encoding----->u*G & distorted by AWGN Channel
            {
                sum = 0;
                for (l = 0; l < k; l++)
                    sum += ((*(u + l)) * (*(G + (l*m + j))));
                *(x_s + j) = sum % 2;

            /* N0 = 2*sigma^2, SNR_b = E_b / N0. */
                *(x_ss + j) = sqrt(R * E_b) * (-2 * (*(x_s + j)) + 1); //BPSK Modulation

                noise_var = E_b / (2.0 * pow(10, (SNR_b / 10.0)));
                noise = sqrt(noise_var) * Gaussian();

                *(y_r + j) = (*(x_ss + j)) + noise;
                *(cm_int + j) = 4.0 * ((sqrt(R*E_b)) / (2.0*noise_var)) * (*(y_r + j)) ;

            }

            convergence = SPA(cm_int, n, m, row_w, col_w, variable, check, max_iterations, &decoded_x, P);

//if the algorithm doesn't convergent to a codeword, then we assume that it's all-zero codewords;
            if (convergence == 0)
            {
                failure ++;
                hd = 0;
                for (j = 0; j < m; j++)
                    hd += *(x_s + j);
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
        fprintf(WBER, "%lf %lf\n",  SNR_b, (double)failure / total_trial);
        fprintf(WBER, "%lf %lf\n",  SNR_b, (double)error_bits / (total_trial*m));
    }

    fclose(WBER);

    free(P);
    free(cm_int);
    free(u);
    free(decoded_x);
    free(x_ss);
    free(x_s);
    free(y_r);
    free(G);
    free(variable);
    free(check);

    return 0 ;

}
