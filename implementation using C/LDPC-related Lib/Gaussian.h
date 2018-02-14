#ifndef GAUSSIAN_H_INCLUDED
#define GAUSSIAN_H_INCLUDED
/*
    Using Box-Muller transform to obtain the standard Gaussian Distribution from the (0,1) uniform Distribution
*/

#define PI 3.14159265359

#include <stdlib.h>
#include <math.h>

float Gaussian();


inline float Gaussian()
{
    float U = 0, V = 0, z = 0;
    U = rand() / (RAND_MAX + 1.0);
    V = rand() / (RAND_MAX + 1.0);
    z = sqrt(-2.0 * log(U)) * sin(2.0 * PI * V);

    return z;
}

#endif // GAUSSIAN_H_INCLUDED
