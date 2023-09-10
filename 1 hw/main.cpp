#include <iostream>
// #include "header.h"
#include <cmath>

using namespace std;

double prod(int N, double* coef, double* Y_m, double* Y_n) // скалярное произведение 2-х векторов разм. N+1
{
    double prod = 0;
    for (int i = 0; i < N+1; ++i)
    {
        prod += coef[i]*Y_m[i]*Y_n[i]; 
    }

    return prod;
}

double max_prod(int N)
{
    double* coef;
    coef = new double[N+1];
    double* Y_m;
    Y_m = new double[N+1];
    double* Y_n;
    Y_n = new double[N+1];

    for (int i = 0; i < N+1; ++i) // 0 1 1 ... 1 0
    {
        if ((i == 0) || (i == N))
            coef[i] = 0;
        else 
            coef[i] = 4/(2*N-1);
    }

    double max_prod = 0;
    for (int n = 0; n < N-1; ++n)
        {
            for (int m = 0; m < N-1; ++m)
            {
                if (m != n)
                {
                    for (int k = 0; k < N+1; ++k)
                    {
                        Y_m[k] = sin((M_PI*(2*m+1)*k)/(2*N-1));
                        Y_n[k] = sin((M_PI*(2*n+1)*k)/(2*N-1));
                    }
                    max_prod = max(max_prod, prod(N, coef, Y_m, Y_n));
                }
            }
        }

    return max_prod;
}

double max_norm(int N)
{
    double h = (double) 1/N;
    double* coef;
    coef = new double[N+1];
    double* Y_m;
    Y_m = new double[N+1];
    double* Y_n;
    Y_n = new double[N+1];

    for (int i = 0; i < N+1; ++i) // 0 1 1 ... 1 0
    {
        if ((i == 0) || (i == N))
            coef[i] = 0;
        else 
            coef[i] = 4/(2*N-1);
    }

    double max_prod = 0;
    for (int m = 0; m < N-1; ++m)
        {
            for (int k = 0; k < N+1; ++k)
            {
                Y_m[k] = sin((M_PI*(2*m+1)*k)/(2*N-1));
            }

            double* vec;
            vec = new double[N+1];
            for (int k = 0; k <= N; ++k)
            {
                if ((k == 0) || (k==N))
                    vec[k] = 1*Y_m[k];
                else
                    vec[k] = 1/(h*h)*Y_m[k-1] - 2/(h*h)*Y_m[k] + 1/(h*h)*Y_m[k+1];

                vec[k] += 4/(h*h)*sin(M_PI*(2*m+1)/(2*(2*N-1))) * Y_m[k];
            }

            max_prod = max(max_prod, sqrt(prod(N, coef, vec, vec)) / (4/(h*h)*pow(sin(M_PI*(2*m+1)/(2*(2*N-1))), 2)) );
            // cout << (4/(h*h)*pow(sin(M_PI*(2*m+1)/(2*(2*N-1))), 2)) << endl;
        }

    return max_prod;
}

int main()
{
    int N;
    for (int k = 2; k < 13; ++k)
    {
        N = pow(2, k);

        cout << "Максимальное произведение для k=" << k << ": "<< max_prod(N) << endl;   
        cout << "Максимальное норма для k=" << k << ": "<< max_norm(N) << endl << endl;  
    }

    return 0;
}