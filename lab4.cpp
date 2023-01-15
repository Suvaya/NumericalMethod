#include <iostream>
#include <cmath>
using namespace std;

const int N = 3; // Number of equations

/*
//gauss elimination
void GaussianElimination(double a[][N], double b[], double x[]) {
    // Forward elimination
    for (int k = 0; k < N - 1; k++) {
        for (int i = k + 1; i < N; i++) {
            double factor = a[i][k] / a[k][k];
            b[i] -= factor * b[k];
            for (int j = k; j < N; j++) {
                a[i][j] -= factor * a[k][j];
            }
        }
    }

    // Back substitution
    for (int i = N - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < N; j++) {
            sum += a[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / a[i][i];
    }
}

int main() {
    double a[N][N] = { { 3, 2, -1 }, { 2, -2, 4 }, { -1, 0.5, -1 } };
    double b[N] = { 1, -2, 0 };
    double x[N];

    GaussianElimination(a, b, x);

    cout << "The solution is:" << endl;
    for (int i = 0; i < N; i++) {
        cout << "x[" << i << "] = " << x[i] << endl;
    }

    return 0;
}

*/
/*
//gauss jordan

void GaussianJordan(double a[][N], double b[], double x[]) {
    // Forward elimination and backward elimination
    for (int k = 0; k < N; k++) {
        // pivot
        double pivot = a[k][k];
        for (int j = 0; j < N; j++) {
            a[k][j] /= pivot;
        }
        b[k] /= pivot;

        // elimination
        for (int i = 0; i < N; i++) {
            if (i != k) {
                double factor = a[i][k];
                for (int j = 0; j < N; j++) {
                    a[i][j] -= factor * a[k][j];
                }
                b[i] -= factor * b[k];
            }
        }
    }

    // Copy solution to x[]
    for (int i = 0; i < N; i++) {
        x[i] = b[i];
    }
}

int main() {
    double a[N][N] = { { 3, 2, -1 }, { 2, -2, 4 }, { -1, 0.5, -1 } };
    double b[N] = { 1, -2, 0 };
    double x[N];

    GaussianJordan(a, b, x);

    cout << "The solution is:" << endl;
    for (int i = 0; i < N; i++) {
        cout << "x[" << i << "] = " << x[i] << endl;
    }

    return 0;
}

*/

/*
//Jacobi iteration

const double TOL = 1e-6; // Tolerance for convergence

void Jacobi(double a[][N], double b[], double x[]) {
    double x_prev[N];
    double error;

    do {
        for (int i = 0; i < N; i++) {
            x_prev[i] = x[i];
        }

        for (int i = 0; i < N; i++) {
            double sum = 0;
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    sum += a[i][j] * x_prev[j];
                }
            }
            x[i] = (b[i] - sum) / a[i][i];
        }

        error = 0;
        for (int i = 0; i < N; i++) {
            error += abs(x[i] - x_prev[i]);
        }
    } while (error > TOL);
}

int main() {
    double a[N][N] = { { 3, 2, -1 }, { 2, -2, 4 }, { -1, 0.5, -1 } };
    double b[N] = { 1, -2, 0 };
    double x[N] = { 0, 0, 0 }; // initial guess

    Jacobi(a, b, x);

    cout << "The solution is:" << endl;
    for (int i = 0; i < N; i++) {
        cout << "x[" << i << "] = " << x[i] << endl;
    }

    return 0;
}

*/

/*
//Gauss Seidel
const double TOL = 1e-6; // Tolerance for convergence

void GaussSeidel(double a[][N], double b[], double x[]) {
    double error;

    do {
        error = 0;
        for (int i = 0; i < N; i++) {
            double sum = 0;
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    sum += a[i][j] * x[j];
                }
            }
            double prev_x = x[i];
            x[i] = (b[i] - sum) / a[i][i];
            error += abs(prev_x - x[i]);
        }
    } while (error > TOL);
}

int main() {
    double a[N][N] = { { 3, 2, -1 }, { 2, -2, 4 }, { -1, 0.5, -1 } };
    double b[N] = { 1, -2, 0 };
    double x[N] = { 0, 0, 0 }; // initial guess

    GaussSeidel(a, b, x);

    cout << "The solution is:" << endl;
    for (int i = 0; i < N; i++) {
        cout << "x[" << i << "] = " << x[i] << endl;
    }

    return 0;
}

*/

/*
//Eigen values using power method

const double TOL = 1e-6; // Tolerance for convergence

void powerMethod(double a[][N], double &eigenvalue, double x[]) {
    double lambda_prev;
    double error;

    do {
        // matrix-vector multiplication
        double y[N];
        for (int i = 0; i < N; i++) {
            y[i] = 0;
            for (int j = 0; j < N; j++) {
                y[i] += a[i][j] * x[j];
            }
        }

        // normalize y and update x
        double norm = 0;
        for (int i = 0; i < N; i++) {
            norm += y[i] * y[i];
        }
        norm = sqrt(norm);
        for (int i = 0; i < N; i++) {
            x[i] = y[i] / norm;
        }

        // calculate eigenvalue
        lambda_prev = eigenvalue;
        eigenvalue = 0;
        for (int i = 0; i < N; i++) {
            eigenvalue += x[i] * y[i];
        }

        // check for convergence
        error = abs(eigenvalue - lambda_prev);
    } while (error > TOL);
}

int main() {
    double a[N][N] = { { 2, 1, 1 }, { 1, 2, 1 }, { 1, 1, 2 } };
    double eigenvalue;
    double x[N] = { 1, 1, 1 }; // initial guess

    powerMethod(a, eigenvalue, x);

    cout << "The eigenvalue is: " << eigenvalue << endl;

    return 0;
}

*/