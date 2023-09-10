#include <iostream>
#include "header.h"
#include <math.h>

using namespace std;


void enter_matrix(int n, int form, double a, double* A, double* B, double* C, double* F) //Функция, вызываемая из main. Возваращает 3 массива со значениями a, b, c
{
    double h = ((double) 1)/((double) n);

    switch (form) 
    {
        case 1:
        {
            for (int i = 1; i < n+1; ++i)
            {
                A[i] = ((double) -1)/h + a;
                B[i] = (double) 0;
                C[i] = ((double) 1)/h;
                F[i] = (double) 0;
            }
        break;
        }

        case 2:
        {
            for (int i = 1; i < n+1; ++i)
            {
                A[i] = ((double) -1)/h;
                B[i] = (double) 0;
                C[i] = ((double) 1)/h + a;
                F[i] = (double) 0;
            }
        break;
        }

        case 3:
        {
            for (int i = 1; i < n+1; ++i)
            {
                A[i] = ((double) -1)/h + a * 0.5;
                B[i] = (double) 0;
                C[i] = ((double) 1)/h + a * 0.5;
                F[i] = (double) 0;
            }
        break;
        }

        case 4:
        {
            for (int i = 2; i < n+1; ++i)
            {
                A[i] = ((double) -1)/(2*h);
                B[i] = ((double) 1)/(2*h);
                C[i] = (double) a;
                F[i] = (double) 0;
            }
        break;
        }

        case 5:
        {
            for (int i = 2; i < n+1; ++i)
            {
                A[i] = ((double) 0.5)/h;
                B[i] = ((double) 1.5)/h + a;
                C[i] = ((double) -2)/h;
                F[i] = (double) 0;
            }
        break;
        }

        case 6:
        {
            for (int i = 2; i < n+1; ++i)
            {
                A[i] = ((double) -1.5)/h + a;
                B[i] = ((double) -0.5)/h;
                C[i] = ((double) 2)/h;
                F[i] = (double) 0;
            }
        break;
        }

        default:
            cout << "Error: incorrect formule";
        break;
    }
}