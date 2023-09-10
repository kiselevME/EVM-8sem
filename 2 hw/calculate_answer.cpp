#include <iostream>
#include "header.h"

using namespace std;

void calculate_answer(int n, double a, double* A, double* B, double* C, double* F, double* Y, int mode_flg)
{
    if (mode_flg == 0) // C_k*Y_k + A_k*Y_{k-1} = F_k -- выражаем Y_k
    {
        Y[0] = 1;
        for (int k = 1; k <= n; ++k)
        {
            Y[k] = (F[k] - A[k]*Y[k-1])/C[k];
        }
    }
    else // B_k*Y_k + C_k*Y_{k-1} + A_k*Y_{k-2} = F_k -- выражаем Y_k
    {
        Y[0] = 1;
        Y[1] = 1 - a/((double) n);
        for (int k = 2; k <= n; ++k)
        {
            Y[k] = - (C[k]*Y[k-1] + A[k]*Y[k-2] - F[k]) / B[k];
        }
    }
}
