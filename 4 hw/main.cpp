#include <iostream>
// #include "header.h"
#include <cmath>
#include <fstream>

using namespace std;


double f(double t, double x, int type)
{
    switch (type)
    {
        case 0:
            return 0;
            break;
        // case 1:
        //     return x;
        //     break;
    }

    return -1;
}

double* progonka(int n, double* A, double* B, double* C, double* F)
{
    // Считаем полином
    double* alpha;
    alpha = new double [n+1];
    double* beta;
    beta = new double [n+1];

    alpha[1] = B[0]/C[0];
    beta[1] = F[0]/C[0];
    for (int k = 1; k < n; ++k)
    {
        alpha[k+1] = B[k] / (C[k] - A[k]*alpha[k]);
        beta[k+1] = (F[k] + A[k]*beta[k]) / (C[k] - A[k]*alpha[k]);
    }
        
    double* Y;
    Y = new double[n+1];

    Y[n] = (F[n] + A[n]*beta[n]) / (C[n] - A[n]*alpha[n]);
    for (int k = n-1; k >= 0; --k)
    {
        Y[k] = alpha[k+1]*Y[k+1] + beta[k+1];
    }

return Y;
}



void solver(double** u_nm, double* x_m, double* t_n, int M, int N, int sheme, int type_u0, int type_f)
{
    double h = ((double) 1) / ((double) M);
    double tau = ((double) 1) / ((double) N);
    double ro = tau/pow(h,2);

    // слой n = 0
    switch (type_u0)
    {
        case 0:
        {
            for (int j = 0; j <= M; ++j)
                u_nm[0][j] = x_m[j];
            break;
        }
        case 1:
        {
            for (int j = 0; j <= M; ++j)
                u_nm[0][j] = exp(x_m[j]) - 1;
            break;
        }
        case 2:
        {
            for (int j = 0; j <= M; ++j)
                u_nm[0][j] = 0;
            break;
        }
        case 3:
        {
            for (int j = 0; j <= M; ++j)
                u_nm[0][j] = sin(3.0/2.0*M_PI*x_m[j]);
            break;
        }
    } 

        // вычисляем остальные точки
    if (sheme == 0)
    {
        for (int i = 1; i <= N; ++i)
        {
            u_nm[i][0] = 0;

            for (int j = 1; j <= M-1; ++j)
                u_nm[i][j] = ro*u_nm[i-1][j+1] + (1-2*ro)*u_nm[i-1][j] + ro*u_nm[i-1][j-1] + tau*f(t_n[i], x_m[j], type_f);

            u_nm[i][M] = 4.0/3.0*u_nm[i][M-1] - 1.0/3.0*u_nm[i][M-2];
        }
    }
    else if (sheme == 1)
    {
        // double* tmp;
        // tmp = new double[M-1];

        double* A;
        double* B;
        double* C;
        double* F;
        A = new double[M+1];
        B = new double[M+1];
        C = new double[M+1];
        F = new double[M+1];

        for (int i = 1; i <= N; ++i)
        {
            u_nm[i][0] = 0;

            B[0] = 0;
            C[0] = 1;
            F[0] = 0;
            for (int j = 1; j <= M-1; ++j)
            {
                A[j] = 1;
                B[j] = 1;
                C[j] = 2 + pow(h,2)/tau;
                F[j] = pow(h,2)/tau*u_nm[i-1][j] + tau*f(t_n[i+1], x_m[j], type_f);
            }
            A[M] = 1;
            C[M] = 1;
            F[M] = 0;

            u_nm[i] = progonka(M, A, B, C, F);

            // u_nm[i][M] = (1 - 2*tau/(h*h))*u_nm[i-1][M] + 2*tau/(h*h)*u_nm[i-1][M-1];
            u_nm[i][M] = 4.0/3.0*u_nm[i][M-1] - 1.0/3.0*u_nm[i][M-2];
        }
    }
}


int main()
{
    int sheme;
    cout << "Введите тип схемы (0 – явная, 1 – неявная): ";
    cin >> sheme;

    int type_u0;
    cout << "Введите тип u^0(x) (0 – x, 1 – exp(x)-1, 2 – 0, 3 – sin(3*pi*x/2)): ";
    cin >> type_u0;

    int type_f;
    cout << "Введите тип f(t,x) (0 – 0): ";
    cin >> type_f;

    int M;
    cout << "Введите M (число узлов по x): ";
    cin >> M;

    int N;
    cout << "Введите N (число узлов по t): ";
    cin >> N;

    double h = ((double) 1) / ((double) M);
    double tau = ((double) 1) / ((double) N);

    double** u_nm;
    u_nm = new double*[N+1];
    double* x_m;
    x_m = new double[M+1];
    double* t_n;
    t_n = new double[N+1];
    for (int i = 0; i <= N; ++i)
    {
        u_nm[i] = new double[M+1]; // определяю u
        t_n[i] = tau*i;
    }

    for (int j = 0; j <= M; ++j)
        x_m[j] = h*j;

    solver(u_nm, x_m, t_n, M, N, sheme, type_u0, type_f);

    // for (int i = 0; i <= N; ++i)
    // {
    //     for (int j = 0; j <= M; ++j)
    //         cout << u_nm[i][j] << ' ';
    //     cout << endl;
    // }
    
    // for (int i = 0; i <= N; ++i)
    // {
    //     for (int j = 0; j <= M; ++j)
    //         cout << exp(-pow(3.0/2.0*M_PI,2)*t_n[i])*sin(3.0/2.0*M_PI*x_m[j]) << ' ';
    //     cout << endl;
    // }

    double max_abs_error = 0;
    double abs_error;
    for (int i = 0; i <= N; ++i)
    {
        for (int j = 0; j <= M; ++j)
            switch (type_f)
            {
                case 0:
                    abs_error = abs(u_nm[i][j] - exp(-pow(3.0/2.0*M_PI,2)*t_n[i])*sin(3.0/2.0*M_PI*x_m[j]));
                    break;
                // case 1:
                //     abs_error = abs(u_nm[i][j] - ((1.0/2.0 - pow(x_m[j],3)/6) + exp(-pow(3.0/2.0*M_PI,2)*t_n[i])*sin(3.0/2.0*M_PI*x_m[j])));
                //     break;
            }
            if (abs_error > max_abs_error)
                max_abs_error = abs_error;
    }
    cout << "Максимальная абсолютная ошибка: " << max_abs_error;

    // fstream fout;
    // fout.open("data.dat", fstream::out);
    // if (!fout)
    // {
    //     cout << "Can't open file";
    // }
    // else
    // {
    //     for (int i = 0; i <= N; ++i)
    //     {
    //         fout << x_n[i][0] << ' ' << y_n[i][0] << endl;
    //     }
    // }
    // fout.close();

    delete [] u_nm;
    delete [] x_m;
    delete [] t_n;

    return 0;
}