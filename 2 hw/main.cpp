#include <iostream>
#include "header.h"
#include <math.h>

using namespace std;

double launch(int n, double a, int form)
{
    // cout << "n=" << n << "    A=" << a << "    form=" << form << endl;

    int N = (int) pow(10, n);

    double* A;
    A = new double[N+1];
    double* B;
    B = new double[N+1];
    double* C;
    C = new double[N+1];
    double* F;
    F = new double[N+1];
    enter_matrix(N, form, a, A, B, C, F);

    // cout << "Введенная матрица: " << endl;
    // for (int i = 0; i <= N; ++i)
    // {
    //     for (int j = 0; j <= N; ++j)
    //     {
    //         if (i == j+1)
    //             cout << -A[i] << ' ';
    //         else if (i == j)
    //             cout << C[i] << ' ';
    //         else if (i == j-1)
    //             cout << -B[i] << ' ';
    //         else
    //             cout << 0 << ' ';
    //     }   
    //     cout << endl;
    // }

    double* Y;
    Y = new double[N+1];
    if (form < 4)
        calculate_answer(N, a, A, B, C, F, Y, 0);
    else
        calculate_answer(N, a, A, B, C, F, Y, 1);

    // cout << "Полученный ответ: " << endl;
    // for (int i = 0; i <= N; ++i)
    // {
    //     cout << Y[i] << ' ';
    // }

    double h = ((double) 1)/((double) N);
    double max = 0;
    for (int i = 0; i <= N; ++i)
    {
        if (abs(exp(-a*((double) i)*h) - Y[i]) > max)
            max = abs(exp(-a*((double) i)*h) - Y[i]);
    }
    // cout << "Норма: " << max << endl << endl;

    delete []A;
    delete []B;
    delete []C;
    delete []F;

    return max;
}




int main()
{
    int mode;
    cout << "Введите тип исполнения(0 – с ручным выбором n и A; 1 - для всех n, всех A): ";
    cin >> mode;

    if (mode == 1)
    {
        double e1, e2, e3, e6;
        cout << "№    " << "E_1    " << "E_2    " << "E_3    " << "E_6    " << "m    " << "A    " << endl;
        for (int form = 1; form <= 6; ++form)
        {
            e1 = launch(1, 1, form);
            e2 = launch(2, 1, form);
            e3 = launch(3, 1, form);
            e6 = launch(6, 1, form);
            cout << form << "    " << e1 << "    " << e2 << "    " << e3 << "    " << e6 << "    " << log(e2 / e3) / log(10) << "    " << 1 << endl;
            e1 = launch(1, 10, form);
            e2 = launch(2, 10, form);
            e3 = launch(3, 10, form);
            e6 = launch(6, 10, form);
            cout << form << "    " << e1 << "    " << e2 << "    " << e3 << "    " << e6 << "    " << log(e2 / e3) / log(10) << "    " << 10 << endl;
            e1 = launch(1, 1000, form);
            e2 = launch(2, 1000, form);
            e3 = launch(3, 1000, form);
            e6 = launch(6, 1000, form);
            cout << form << "    " << e1 << "    " << e2 << "    " << e3 << "    " << e6 << "    " << log(e2 / e3) / log(10) << "    " << 1000 << endl;
        }
        return 0;
    }
    else
    {
        int n;
        cout << "Введите n: ";
        cin >> n;

        double a;
        cout << "Введите A: ";
        cin >> a;

        int form;
        cout << "Введите номер формулы (от 1 до 6): ";
        cin >> form;

        cout << "Норма: " << launch(n, a, form) << endl;
        return 0;
    }
}