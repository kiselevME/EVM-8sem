#include <iostream>
// #include "header.h"
#include <cmath>
#include <fstream>

using namespace std;


double calculate_f(double x, double* y, int type, int n_eq)
{
    if (type == 1)
    {
        if (n_eq == 0)
            return x;
        else if (n_eq == 1)
            return x*x;
        else if (n_eq == 2)
            return x*x*x;
        else if (n_eq == 3)
            return x*x*x*x;
    }
    else if (type == 2)
        return exp(x);
    else if (type == 3)
        return 1 + 2*y[n_eq] - x*x;
    else if (type == 4)
    {
        if (n_eq == 0)
            return y[1];
        else if (n_eq == 1)
            return -y[0];
    }

    return -1;
}


double calculate_y0(int type, int n_eq)
{
    if (type == 1)
    {
        if (n_eq == 0)
            return 0;
        else if (n_eq == 1)
            return 0;
        else if (n_eq == 2)
            return 0;
        else if (n_eq == 3)
            return 0;
    }
    else if (type == 2)
        return 1;
    else if (type == 3)
        return 1;
    else if (type == 4)
    {
        if (n_eq == 0)
            return 0;
        else if (n_eq == 1)
            return 1;
    }

    return -1;
}


void solver(double** y_n, double** x_n, int dim, int type, int N)
{
    double* k;
    k = new double[dim*4]; // вектор 4 * dim
    double* y_;
    y_ = new double[dim];

    double h = ((double) 1) / ((double) N);
    for (int i = 0; i < N; ++i)  // по x_i
    {
        for (int j = 0; j < dim; ++j)
            k[dim*0+j] = h*calculate_f(x_n[i][j], y_n[i], type, j);    // k[0][j]
        for (int j = 0; j < dim; ++j)
            y_[j] = y_n[i][j] + k[dim*0+j]/2.0;
        for (int j = 0; j < dim; ++j)
            k[dim*1+j] = h*calculate_f(x_n[i][j]+h/2.0, y_, type, j);  // k[1][j]
        for (int j = 0; j < dim; ++j)
            y_[j] = y_n[i][j] + k[dim*1+j]/2.0;
        for (int j = 0; j < dim; ++j)
            k[dim*2+j] = h*calculate_f(x_n[i][j]+h/2.0, y_, type, j);  // k[2][j]
        for (int j = 0; j < dim; ++j)
            y_[j] = y_n[i][j] + k[dim*2+j];                            // k[3][j]
        for (int j = 0; j < dim; ++j)
        {
            k[dim*3+j] = h*calculate_f(x_n[i][j]+h, y_, type, j);                
            y_n[i+1][j] = y_n[i][j] + 1.0/6.0 * (k[dim*0+j] + 2*k[dim*1+j] + 2*k[dim*2+j] + k[dim*3+j]);
        }
    }
}


int main()
{
    int type;
    cout << "Введите номер системы ODE(1 – полиномы f(x,y) = [x, x^2, x^3, x^4]; 2 – f(x,y) = exp(x)); 3 – f(x,y) = 1 + 2y - x^2; 4 – f(x, y) = [y_1, -y_0]): ";
    cin >> type;

    int N;
    cout << "Введите N: ";
    cin >> N;

    double h = ((double) 1) / ((double) N);
    
    int dim;
    if (type == 1)
        dim = 4;
    else if (type == 2)
        dim = 1;
    else if (type == 3)
        dim = 1;
    else if (type == 4)
        dim = 2;

    double** y_n;
    y_n = new double*[N+1];
    double** x_n;
    x_n = new double*[N+1];
    for (int i = 0; i <= N; ++i)
    {
        y_n[i] = new double[dim];
        x_n[i] = new double[dim];
        for (int j = 0; j < dim; ++j)
            x_n[i][j] = ((double) i) * h; // беру сетку на [0;1]^dim
    }

    for (int j = 0; j < dim; ++j)
        y_n[0][j] = calculate_y0(type, j);

    solver(y_n, x_n, dim, type, N);

    // for (int i = 0; i <= N; ++i)
    // {
    //     for (int j = 0; j < dim; ++j)
    //     {
    //         cout << y_n[i][j] << ' ';
    //     }
    //     cout << endl;
    // }

    cout << "Порядок сходимости: " << endl;
    double max_diff = 0;
    for (int i = 0; i <= N; ++i)
    {
        double real_answer;
        for (int j = 0; j < dim; ++j)
        {
            if (type == 1)
            {
                if (j == 0)
                    real_answer = (1.0/2.0)*pow(x_n[i][j], 2);
                else if (j == 1)
                    real_answer = (1.0/3.0)*pow(x_n[i][j], 3);
                else if (j == 2)
                    real_answer = (1.0/4.0)*pow(x_n[i][j], 4);
                else if (j == 3)
                    real_answer = (1.0/5.0)*pow(x_n[i][j], 5);
            }
            else if (type == 2)
                real_answer = exp(x_n[i][j]);
            else if (type == 3)
                real_answer = ((2*x_n[i][j]*x_n[i][j] + 2*x_n[i][j] - 1)*exp(-2*x_n[i][j]) + 5) / (4*exp(-2*x_n[i][j]));
            else if (type == 4)
            {
                if (j == 0)
                    real_answer = sin(x_n[i][j]);
                else if (j == 1)
                    real_answer = cos(x_n[i][j]);
            }

            if (abs(real_answer - y_n[i][j]) > max_diff)
                max_diff = abs(real_answer - y_n[i][j]);
        }
    }

    cout << "max_diff: " << max_diff << endl;
    cout << "max_diff / h: " << max_diff/h << endl;
    cout << "max_diff / h^2: " << max_diff/(h*h) << endl;
    cout << "max_diff / h^3: " << max_diff/(h*h*h) << endl;
    cout << "max_diff / h^4: " << max_diff/(h*h*h*h) << endl;
    // cout << "max_diff / h^5: " << max_diff/(h*h*h*h*h) << endl;

    fstream fout;
    fout.open("data.dat", fstream::out);
    if (!fout)
    {
        cout << "Can't open file";
    }
    else
    {
        for (int i = 0; i <= N; ++i)
        {
            fout << x_n[i][0] << ' ' << y_n[i][0] << endl;
        }
    }
    fout.close();

    delete [] y_n;
    delete [] x_n;

    return 0;
}