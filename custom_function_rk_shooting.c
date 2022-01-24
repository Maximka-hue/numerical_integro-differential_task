//2019 year
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>


#define nn 5
#define kk 2
#define ALPHA 10

double f(double x, double* y, double shift_y, int j);
double increment(double x_old, double* y_old, double h, int j);
double* RK(FILE *file, bool print, double x_l, double x_r, double* y_0, double* y_new, double tol);
double* newton(double x_l, double x_r, double* y_0, double* y_T, double* y_new, double* alpha, double tol);
double* X(double x_l, double x_r, double* y_0, double* y_T, double* y_new,  double* X_tmp, double tol);
double* solution(double a[10][10], int var, double *x);
double norma(double* V, int dim);
int main(void) {
    FILE *file = fopen("/Users/kirill/CLionProjects/task_2/out_43/rk_out.txt", "w");
    double x_r = M_PI/2, x_l = 0., R,R0, R1;
    double *y_new = (double* ) malloc((nn) * sizeof(double));
    double *y_0 = (double* ) malloc((nn) * sizeof(double));
    double *y_T = (double* ) malloc((kk) * sizeof(double));
    double *y = (double* ) malloc((nn) * sizeof(double));
    double *alpha = (double* ) malloc((kk) * sizeof(double));
    double tol = pow(10, -10);
    //R = 0.00000315013342683939
    //R = 0.00000456708272662080

    y_0[1] = 0; // x2(0)
    y_0[2] = 0; // p1(0)
    y_0[4] = 0; // B0(0)

    y_T[0] = 0; // x1(pi/2)
    y_T[1] = 0; // p2(pi/2)

    alpha[0] = 1; // x1(0)
    alpha[1] = 1 ; // p2(0)

//    alpha = newton(x_l, x_r, y_0, y_T, y_new, alpha, tol);
    y_0[0] = alpha[0];
    y_0[3] = alpha[1];
    printf("\nshooting parameters: alpha_1 = %f, alpha_2 = %f\n", alpha[0], alpha[1]);

    y = RK(file, true, x_l, x_r, y_0, y_new, tol);
    printf("x1 =  %.10f, x2 =  %.10f , p1 =  %.10f , p2 =  %.10f, B0 =  %.10f\n", y[0], y[1], y[2], y[3], y[4]);
    R0 = fabs(y[0] - y_T[0]); // ошибка по x1(pi/2)
    R1 = fabs(y[3] - y_T[1]); // ошибка по p2(pi/2)
    R = sqrt(R0*R0 + R1*R1);
    printf("R0 = %.10f R1 = %.10f\n", R0, R1);
    printf("R = %.20f \n", R);
    free(y_new); free(y_0);
    fclose(file);;
    return 0;
}

double f(double x, double* y, double shift_y, int j) {
    double x1 =  y[0] + shift_y, x2 =  y[1] + shift_y, p1 =  y[2] + shift_y, p2 =  y[3] + shift_y;

    x = x;
    x2 = x2;
    p2 = p2;
    p1 = p1;
    x1 = x1;

    if (j == 0) {
        return x2;
    }
    if (j == 1) {
        return p2 - x1 * exp(- ALPHA * x);
    }
    if (j == 2) {
        return  p2  * exp(- ALPHA * x);
    }
    if (j == 3) {
        return - p1;
    }
    if (j == 4) {
        return p2 * p2; // B0
    }

    return -1;
}
double increment(double x_old, double* y_old, double h, int j) {
    double k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;
    double y_new = 0;
    double C2=1./18, C3=1./12, C4=1./8, C5=5./16, C6=3./8, C7=59./400, C8=93./200,
            C9=5490023248./9719169821., C10=13./20, C11=1201146811./1299019798,
            C12=1, C13=1;
    double A21 = C2, A31=1./48, A32=1./16, A41=1./32, A42=0, A43=3./32, A51=5./16, A52=0, A53 = - 75./64, A54 =- A53,
            A61=3./80, A62 = 0, A63 = 0, A64=3./16, A65=3./20, A71=29443841./614563906, A72 = 0, A73 = 0, A74=77736538./692538347, A75=-28693883./1125000000,A76=23124283./1800000000,
            A81=16016141./946692911, A82 = 0, A83 = 0, A84=61564180./158732637, A85=22789713./633445777, A86=545815736./2771057229.,A87=-180193667./1043307555,
            A91=39632708./573591083, A92 = 0, A93 = 0, A94=-433636366./683701615, A95=-421739975./2616292301., A96=100302831./723423059, A97=790204164./839813087,A98=800635310./3783071287,
            A101=246121993./1340847787, A102 = 0, A103 = 0, A104=-37695042795./15268766246., A105=-309121744./1061227803, A106=-12992083./490766935, A107=6005943493./2108947869, A108=393006217./1396673457, A109=123872331./1001029789,
            A111=-1028468189./846180014, A112 = 0, A113 = 0, A114=8478235783./508512852, A115=1311729495./1432422823, A116=-10304129995./1701304382, A117=-48777925059./3047939560, A118=15336726248./1032824649, A119=-45442868181./3398467696, A1110 = 3065993473./597172653,
            A121=185892177./718116043, A122 = 0, A123 = 0, A124=-3185094517./667107341, A125=-477755414./1098053517, A126=-703635378./230739211, A127=5731566787./1027545527, A128=5232866602./850066563, A129=-4093664535./808688257, A1210=3962137247./1805957418, A1211=65686358./487910083,
            A131=403863854./491063109, A132 = 0, A133 = 0, A134=-5068492393./434740067, A135=-411421997./543043805, A136=652783627./914296604, A137=11173962825./925320556, A138=-13158990841./6184727034, A139=3936647629./1978049680, A1310=-160528059./685178525, A1311 =248638103./1413531060, A1312 = 0;
    double B1=14005451./335480064, B2 =0, B3 = 0, B4 = 0, B5 = 0, B6=-59238493./1068277825, B7=181606767./758867731, B8=561292985./797845732, B9=-1041891430./1371343529, B10=760417239./1151165299, B11=118820643./751138087, B12=-528747749./2220607170., B13=1./4;


    k1 = f(x_old, y_old, 0., j);
    k2 = f(x_old + C2 * h, y_old, h * (A21 * k1), j);
    k3 = f(x_old + C3 * h, y_old, h * (A31 * k1 + A32 * k2), j);
    k4 = f(x_old + C4 * h, y_old, h * (A41 * k1 + A42 * k2 + A43 * k3), j);
    k5 = f(x_old + C5 * h, y_old, h * (A51 * k1 + A52 * k2 + A53 * k3 + A54 * k4), j);
    k6 = f(x_old + C6 * h, y_old, h * (A61 * k1 + A62 * k2 + A63 * k3 + A64 * k4 + A65 * k5), j);
    k7 = f(x_old + C7 * h, y_old, h * (A71 * k1 + A72 * k2 + A73 * k3 + A74 * k4 + A75 * k5 + A76 * k6), j);
    k8 = f(x_old + C8 * h, y_old, h * (A81 * k1 + A82 * k2 + A83 * k3 + A84 * k4 + A85 * k5 + A86 * k6 + A87 * k7), j);
    k9 = f(x_old + C9 * h, y_old, h * (A91 * k1 + A92 * k2 + A93 * k3 + A94 * k4 + A95 * k5 + A96 * k6 + A97 * k7 + A98 * k8), j);
    k10 = f(x_old + C10 * h, y_old, h * (A101 * k1 + A102 * k2 + A103 * k3 + A104 * k4 + A105 * k5 + A106 * k6 + A107 * k7 + A108 * k8 + A109 * k9), j);
    k11 = f(x_old + C11 * h, y_old, h * (A111 * k1 + A112 * k2 + A113 * k3 + A114 * k4 + A115 * k5 + A116 * k6 + A117 * k7 + A118 * k8 + A119 * k9 + A1110 * k10), j);
    k12 = f(x_old + C12 * h, y_old, h * (A121 * k1 + A122 * k2 + A123 * k3 + A124 * k4 + A125 * k5 + A126 * k6 + A127 * k7 + A128 * k8 + A129 * k9 + A1210 * k10 + A1211 * k11), j);
    k13 = f(x_old + C13 * h, y_old, h * (A131 * k1 + A132 * k2 + A133 * k3 + A134 * k4 + A135 * k5 + A136 * k6 + A137 * k7 + A138 * k8 + A139 * k9 + A1310 * k10 + A1311 * k11 + A1312 * k12), j);
    y_new = h * (B1*k1 + B2*k2 + B3*k3 + B4*k4 + B5*k5 + B6*k6 + B7*k7 + B8*k8 + B9*k9 + B10*k10 + B11*k11 + B12*k12 + B13*k13);
    return y_new;
}
double* RK(FILE *file, bool print, double x_l, double x_r, double* y_0, double* y_new, double tol) {
    clock_t begin = clock();
    clock_t end;
    double time_spent;
    double *y_old = (double* ) malloc((nn) * sizeof(double));
    double *w_old = (double* ) malloc((nn) * sizeof(double));
    double *w_new = (double* ) malloc((nn) * sizeof(double));
    double h_new, h_last, h = (x_r - x_l)/1000;
    double x_old, x_new, x_old_2h, x_new_2h, h_last_2h;
    double fac = 0.91, facmin = 1./2.5, facmax = 1.5;
    double err, global_old = 0, global_new = 0, l = 0;
    int i = 1;
    l=l;
    global_new = global_new;
    global_old = global_old;

    for (int s = 0; s < nn; s++) {
        y_new[s] = 0;
        y_old[s] = 0;
        w_new[s] = 0;
        w_old[s] = 0;
    }

    x_old = x_l;
    x_old_2h = x_l;

    for(int j = 0 ; j < nn; j++) {
        y_old[j] = y_0[j];
        w_old[j] = y_0[j];
    }
    if(print) {
        fprintf(file, "%f %f %f %f %f %f\n",x_old, w_old[0], w_old[1], w_old[2], w_old[3], w_old[4]);
    }
    while (x_old < x_r- 1.5*h) {
        x_new = x_old + h;
        for(int j = 0 ; j < nn; j++) y_new[j] = y_old[j] + increment(x_old, y_old, h, j);

        x_old = x_new;
        for(int j = 0 ; j < nn; j++) y_old[j] = y_new[j];

        x_new = x_old + h;
        for(int j = 0 ; j < nn; j++) y_new[j] = y_old[j] + increment(x_old, y_old, h, j);

        x_new_2h = x_old_2h + 2 * h;
        for(int j = 0 ; j < nn; j++) w_new[j] = w_old[j] + increment(x_old_2h, w_old, 2 * h, j);

        err = fmax(fmax( fabs(y_new[0] - w_new[0]), fabs(y_new[1] - w_new[1])),fmax(fabs(y_new[2] - w_new[2]), fabs(y_new[3] - w_new[3])))/31;
        h_new = h * fmin(facmax, fmax(facmin, fac * pow(tol/err, 1./6)));

        if (err <= tol) {
            h = h_new;
            l = 0;
            global_new = err + global_old * exp(h * l);
            global_old = global_new;

            if (print && i%100 == 0) {
                fprintf(file, "%f %f %f %f %f %f\n", x_new, y_new[0], y_new[1], y_new[2], y_new[3], y_new[4]);
            }
        }
        if (err > tol) {
            x_old = x_old_2h;
            for(int j = 0 ; j < nn; j++) y_old[j] = w_old[j];
            h = h_new;
            continue;
        }

        x_old = x_new;
        for(int j = 0 ; j < nn; j++) y_old[j] = y_new[j];

        x_old_2h = x_new_2h;
        for(int j = 0 ; j < nn; j++) w_old[j] = w_new[j];

        i++;

    }

    h_last = x_r - x_old;
    x_new = x_r;
    for(int j = 0 ; j < nn; j++) y_new[j] = y_old[j] + increment(x_old, y_old, h_last, j);

    h_last_2h = x_r - x_old_2h;
    x_new_2h = x_r;
    for(int j = 0 ; j < nn; j++) w_new[j] = w_old[j] + increment(x_old_2h, w_old, h_last_2h, j);

    //printf("global = %.20f\n",global_new);
    if (print && i%100 == 0) {
        fprintf(file, "%f %f %f %f %f %f\n", x_new, y_new[0], y_new[1], y_new[2], y_new[3], y_new[4]);
    }
    free(y_old);free(w_old);free(w_new);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if (print) printf("RK ok and took %.4f sec\n",time_spent);
    return y_new;
}

double* newton(double x_l, double x_r, double* y_0, double* y_T, double* y_new, double* alpha, double tol) {
    double *alpha_old = (double* ) malloc((kk) * sizeof(double));
    double *alpha_new = (double* ) malloc((kk) * sizeof(double));
    double *X_tmp = (double* ) malloc((kk) * sizeof(double));
    double *X_cur = (double* ) malloc((kk) * sizeof(double));
    double *vect_solution = (double* ) malloc((kk) * sizeof(double));
    double norma_old, norma_new, delta = pow(10, -8);
    double A[10][10];
    int k = 0;
    alpha_old[0] = alpha[0];
    alpha_old[1] = alpha[1];
    y_0[0] = alpha_old[0];
    y_0[3] = alpha_old[1];
    X_cur[0] = X(x_l,  x_r, y_0, y_T, y_new, X_tmp, tol)[0];
    X_cur[1] = X(x_l,  x_r, y_0, y_T, y_new, X_tmp, tol)[1];
    norma_new = norma(X_cur, kk);
    while (norma_new >= pow(10, -6)){//for(int k = 0; k < 4; k++)
        printf("iteration number: %d\n",k + 1);
        y_0[0] = alpha_old[0];
        y_0[3] = alpha_old[1];
        X_cur[0] = X(x_l,  x_r, y_0, y_T, y_new, X_tmp, tol)[0];
        X_cur[1] = X(x_l,  x_r, y_0, y_T, y_new, X_tmp, tol)[1];
        norma_old = norma(X_cur, kk);
        printf("norma_old = %.20f \n", norma_old);
        //первый столбец
        alpha_new[0] = alpha_old[0] + delta;
        alpha_new[1] = alpha_old[1];

        y_0[0] = alpha_new[0];
        y_0[3] = alpha_new[1];

        X_tmp[0] = X(x_l,  x_r, y_0, y_T, y_new, X_tmp, tol)[0];
        X_tmp[1] = X(x_l,  x_r, y_0, y_T, y_new, X_tmp, tol)[1];

        A[0][0] = (X_tmp[0] - X_cur[0])/delta;
        A[1][0] = (X_tmp[1] - X_cur[1])/delta;
        //второй столбец
        alpha_new[0] = alpha_old[0];
        alpha_new[1] = alpha_old[1] + delta;

        y_0[0] = alpha_new[0];
        y_0[3] = alpha_new[1];

        X_tmp[0] = X(x_l,  x_r, y_0, y_T, y_new, X_tmp, tol)[0];
        X_tmp[1] = X(x_l,  x_r, y_0, y_T, y_new, X_tmp, tol)[1];

        A[0][1] = (X_tmp[0] - X_cur[0])/delta;
        A[1][1] = (X_tmp[1] - X_cur[1])/delta;

        //правая часть
        A[0][2] = - X_cur[0];
        A[1][2] = - X_cur[1];

        //решение системы
        vect_solution[0] = solution(A,kk,vect_solution)[0];
        vect_solution[1] = solution(A,kk,vect_solution)[1];

        //printf("Vect_solution = %.20f \n", vect_solution[0]);
        alpha_new[0] = alpha_old[0] + vect_solution[0];
        alpha_new[1] = alpha_old[1] + vect_solution[1];
        printf("alpha_new = (%.20f, %.20f) \n", alpha_new[0], alpha_new[1]);

        y_0[0] = alpha_new[0];
        y_0[3] = alpha_new[1];

        X_cur[0] = X(x_l,  x_r, y_0, y_T, y_new, X_tmp, tol)[0];
        X_cur[1] = X(x_l,  x_r, y_0, y_T, y_new, X_tmp, tol)[1];
        norma_new = norma(X_cur, kk);
        printf("norma_new = %.20f \n", norma_new);
        while (norma_new > norma_old) {
            printf("issaev-sonin...\n");
            vect_solution[0] = vect_solution[0]/2;
            vect_solution[1] = vect_solution[1]/2;
            alpha_new[0] = alpha_old[0] + vect_solution[0];
            alpha_new[0] = alpha_old[0] + vect_solution[0];

            y_0[0] = alpha_new[0];
            y_0[3] = alpha_new[1];

            X_cur[0] = X(x_l,  x_r, y_0, y_T, y_new, X_tmp, tol)[0];
            X_cur[1] = X(x_l,  x_r, y_0, y_T, y_new, X_tmp, tol)[1];
            norma_new = norma(X_cur, kk);
            //printf("alpha_new = %.20f \n", alpha_new[0]);
        }
        alpha_old[0] = alpha_new[0];
        alpha_old[1] = alpha_new[1];
        k++;
    }
    alpha[0] = alpha_new[0];
    alpha[1] = alpha_new[1];

    //printf("Vcur = %.20f \nVcur = %.20f \n", norma_new, norma_old);
    free(alpha_new);free(alpha_old);free(X_cur);free(X_tmp);
    return alpha;
}
double* X(double x_l, double x_r, double* y_0, double* y_T, double* y_new,  double* X_tmp, double tol) {
    FILE* file_hz = fopen("/Users/kirill/CLionProjects/task_2/out_43/hz.txt", "w");
    X_tmp[0] = RK(file_hz, false, x_l, x_r, y_0, y_new, tol)[0] - y_T[0];
    X_tmp[1] = RK(file_hz, false, x_l, x_r, y_0, y_new, tol)[3] - y_T[1];
    return X_tmp;
}

double norma(double* V, int dim) {
    double sum = 0;
    for (int k = 0; k < dim; k++) {
        sum += V[k] * V[k];
    }
    return sqrt(sum);
}
double* solution(double a[10][10], int var, double *x) {
    int k, i,  j;
    double l;

    for ( k = 0;k < var;k++ )
    {
        for ( i = 0;i <= var;i++ )
        {
            l = a[ i ][ k ];

            for ( j = 0;j <= var;j++ )
            {
                if ( i != k )
                    a[i][j] = (a[k][k]*a[i][j])-(l*a[k][j]);
            }
        }
    }

    //

    for ( i = 0;i < var;i++ )
    {
        x[i] = ( float ) a[ i ][ var ] / ( float ) a[ i ][ i ];
    }

    return x;

}
