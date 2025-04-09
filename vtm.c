#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double cnorm(double* y, double* exact, int size) {
    double max = 0.0;
    for (int i = 0; i < size; i++) {
        double err = fabs(y[i] - exact[i]);
        if (err > max) {
            max = err;
        }
    }
    return max;
}

double l2_norm(double* y, double* exact, int size, double h) {
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        double err = y[i] - exact[i];
        sum += err * err;
    }
    return sqrt(h * sum);
}

int main(void) {
    int Ns[] = {10, 20, 40, 80, 160};
    int num_cases = sizeof(Ns) / sizeof(Ns[0]);

    FILE* fout = fopen("errors.csv", "w");

    for (int k = 0; k < num_cases; k++) {
        int N = Ns[k];
        double a = 0.0;
        double b = sin(1.0);
        double h = 1.0 / N;

        double* x = malloc((N + 1) * sizeof(double));
        double* f = malloc((N - 1) * sizeof(double));
        double* rhs = malloc((N - 1) * sizeof(double));
        double* alpha = malloc((N - 1) * sizeof(double));
        double* beta = malloc((N - 1) * sizeof(double));
        double* y = malloc((N + 1) * sizeof(double));
        double* u = malloc((N + 1) * sizeof(double));

        clock_t start_time = clock();

        for (int i = 0; i <= N; i++) {
            x[i] = i * h;
            u[i] = sin(x[i]);
        }

        for (int i = 1; i < N; i++) {
            f[i - 1] = sin(x[i]);
        }

        for (int i = 0; i < N - 1; i++) {
            rhs[i] = f[i];
        }
        rhs[0] += a / (h * h);
        rhs[N - 2] += b / (h * h);

        double A = -1.0 / (h * h);
        double B = 2.0 / (h * h);
        double C = -1.0 / (h * h);

        alpha[0] = -C / B;
        beta[0] = rhs[0] / B;

        for (int i = 1; i < N - 1; i++) {
            double denom = B + A * alpha[i - 1];
            alpha[i] = -C / denom;
            beta[i] = (rhs[i] - A * beta[i - 1]) / denom;
        }

        y[0] = a;
        y[N] = b;
        y[N - 1] = beta[N - 2];

        for (int i = N - 3; i >= 0; i--) {
            y[i + 1] = alpha[i] * y[i + 2] + beta[i];
        }

        double C_norm = cnorm(y, u, N + 1);
        double L2_norm = l2_norm(y, u, N + 1, h);

        clock_t end_time = clock();
        double time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

        fprintf(fout, "%d,%f,%e,%e,%f\n", N, h, C_norm, L2_norm, time);

        free(x); free(f); free(rhs); free(alpha); free(beta); free(y); free(u);
    }

    fclose(fout);
    return 0;
}
