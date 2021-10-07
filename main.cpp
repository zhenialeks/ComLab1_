#include <iostream>
#include <cmath>
#include <vector>
#include <bits/stdc++.h>

using namespace std;

double * phi_;

double phi(double x){
    return exp(1 / x) + cos(2 * x);
}

double f(double x, double y, u_int8_t k){
    return -exp(1 / x) / pow(x, 2) - 2*sin(2*x) + k*(y-phi(x));
}

double ** compute_phi(double a, double b, uint64_t N){
    auto * x = new double[N+1];
    auto * y = new double[N+1];


    double h = (b - a) / N;
    for (uint64_t i = 0; i < N+1; i++){
        x[i] = a + h*i;
        y[i] = phi(a + h*i);
    }
    auto ** answer =  new double * [2];
    answer[0] = x;
    answer[1] = y;
    return answer;
}

double ** compute_Euler(double a, double b, double y_a, uint64_t N,
                                      uint8_t k, double * phi_){

    auto * x = new double[N + 1];
    auto * y = new double[N + 1];
    auto * loc_err = new double[N + 1];

    x[0] = a;
    y[0] = y_a;
    loc_err[0] = 0;
    double h = (b - a) / N;
    double x_prev, y_prev;

    for (uint64_t i = 1; i < N+1; i++){
        x_prev = x[i - 1];
        y_prev = y[i-1];

        x[i] = a + h * i;
        y[i] = y_prev + h * f(x_prev, y_prev, k);
        loc_err[i] = abs(phi_[i] - y[i]);
    }

    auto ** answer = new double*[3];
    answer[0] = x;
    answer[1] = y;
    answer[2] = loc_err;
    return answer;
}

double ** compute_Adams(double a, double b, double y_a, uint64_t N,
                               uint8_t k){
    double h = (b - a) / N;
    double ** tmp = compute_Euler(a, a + h * 3, y_a, 3, k, phi_);

    auto * x = (double*) realloc(tmp[0], sizeof(double)*(N+1));
    auto * y = (double *) realloc(tmp[1], sizeof(double)*(N+1));
    auto *loc_err = (double*) realloc(tmp[2], sizeof(double)*(N+1));

    double f_[N+1];
    for (uint i = 0; i < 4; i++)
        f_[i] = f(a + h * i, y[i], k);

    for (uint64_t i = 4; i < N+1; i++){
        x[i] = a + h * i;
        double tmp_y = y[i-1] + h/24 * (44*f_[i-1] - 59*f_[i-2] + 37*f_[i-3] -9*f_[i-4]);
        y[i] = tmp_y;
        f_[i] = f(x[i], y[i], k);
        loc_err[i] = abs(phi_[i] - tmp_y);
    }

    tmp[0] = x;
    tmp[1] = y;
    tmp[2] = loc_err;

    return tmp;
}

double ** compute_Adams_Bashfort_Moulten(double a, double b, double y_a, uint64_t N,
                               uint8_t k){
    double h = (b - a) / N;
    double ** tmp = compute_Euler(a, a + h * 3, y_a, 3, k, phi_);

    auto * x = (double*) realloc(tmp[0], sizeof(double)*(N+1));
    auto * y = (double *) realloc(tmp[1], sizeof(double)*(N+1));
    auto *loc_err = (double*) realloc(tmp[2], sizeof(double)*(N+1));

    double f_[N+1];
    for (uint8_t i = 0; i < 4; i++)
        f_[i] = f(a + h*i, y[i], k);

    for (uint64_t i = 4; i < N + 1; i++){
        x[i] = a + h * i;
        double y_tmp = y[i-1] + h/24*(55*f_[i-1] - 59*f_[i-2] + 37*f_[i-3] - 9*f_[i-4]);
        y[i] = y[i-1] + h/24*(9*f(x[i-1], y_tmp, k) + 19*f_[i-1] - 5*f_[i-2] + f_[i-3]);
        loc_err[i] = abs(y_tmp - phi_[i]);
    }

    tmp[0] = x;
    tmp[1] = y;
    tmp[2] = loc_err;

    return tmp;

}

double** compute_Runge_Kutt(double a, double b, double y_a, uint64_t N,
                                                uint8_t k, uint8_t order=1) {
    double h = (b - a) / N;

    auto * x = new double[(int) N+1];
    auto * y = new double[(int) N+1];
    auto * loc_err = new double[ (int) N+1];

    x[0] = a;
    y[0] = y_a;
    loc_err[0] = 0;

    if (order == 1){
        for (uint64_t i = 1; i < N+1; i++){
            double k1 = h*f(x[i-1], y[i-1], k);
            double k2 = h*f(x[i-1] + h/2, y[i-1] + k1/2, k);
            double k3 = h*f(x[i-1] + h, y[i-1] - k1 + 2*k2, k);
            x[i] = x[i-1] + h;
            y[i] = y[i-1] + (k1 + 4*k2 + k3)/6;
            loc_err[i] = abs(phi(x[i-1]) - y[i-1]);
        }
    } else if (order == 2){
        for (uint64_t i = 1; i < N+1; i++){
            double k1 = h * f(x[i-1], y[i-1], k);
            double k2 = h * f(x[i-1] + h/2, y[i-1] + k1/2, k);
            double k3 = h*f(x[i-1]+h/2, y[i-1] + k2/2, k);
            double k4 = h*f(x[i-1]+h, y[i-1]+k3, k);
            x[i] = x[i-1] + h;
            y[i] = y[i-1] + (k1 + 2*k2 + 2*k3 + k4) / 6;
            loc_err[i] = abs(phi_[i] - y[i]);
        }
    }

    auto ** answer = new double * [3];
    answer[0] = x;
    answer[1] = y;
    answer[2] = loc_err;
    return answer;

}

double** compute_Gir(double a, double b, double y_a, uint64_t N, uint8_t k){
    double h = (b - a) / N;
    double ** tmp = compute_Runge_Kutt(a, a + h * 3, y_a, 3, k, 2);

    auto *x = (double *) realloc(tmp[0], sizeof(double)*(N+1));
    auto *y = (double *) realloc(tmp[1], sizeof(double)*(N+1));
    auto *loc_err = (double *) realloc(tmp[2], sizeof(double)*(N+1));

    for (uint64_t i = 4; i < N+1; i++){
        double new_x = a + h*i;
        double new_y = (12*h*(f(x[i-1], 0, 0) - k*phi(new_x)) + 48*y[i-1] - 36*y[i-1]
                + 16*y[i-3] - 3*y[i-4]) / (25 - 12 * h * k);
        y[i] = new_y;
        x[i] = new_x;
        loc_err[i] = abs(new_y - phi_[i]);
    }

    auto ** answer = new double * [3];
    answer[0] = x;
    answer[1] = y;
    answer[2] = loc_err;

    return answer;
}

bool write_down_results(uint64_t N, double *const *result, uint8_t res_dim, const string& prefix) {
    ofstream fout;
    fout.open("./" + prefix + to_string(N) + ".txt", ios_base::trunc | ios_base::out);
    if (fout.is_open()){
        fout << "x\ty\tloc_err\tN=" << N << '\n';
        for (uint64_t i = 0; i < N+1; i++) {
            if (res_dim == 1)
                fout << result[0][i] << endl;
            else {
                fout << result[0][0] << '\t';
                for (uint8_t j = 1; j < res_dim; j++)
                    fout << result[0][j] << "\t";
                fout << endl;
            }
        }
    } else
        return false;
    fout.close();
    return true;
}

int main() {
    uint64_t N = 10;
    double a = 1;
    double b = 12;
    uint8_t k = 0;

    double * y_1, *x_1;

    // Generate reference data (phi)
    double y_a = phi(a);
    double ** reference_phi = compute_phi(a, b, 100);

    //Generate helpful data
    double ** tmp_phi = compute_phi(a, b, N);
    phi_ = tmp_phi[1];

    //Main computations
    double ** euler_res = compute_Euler(a, b, y_a, N, k, phi_);
    double ** adams_res = compute_Adams(a, b, y_a, N, k);
    double ** r_k_res = compute_Runge_Kutt(a, b, y_a, N, k);
    double ** a_b_m_res = compute_Adams_Bashfort_Moulten(a, b, y_a, N, k);
    double ** gir_res = compute_Gir(a, b, y_a, N, k);

    x_1 = euler_res[0];
    y_1 = euler_res[1];


    // Compute error using Runge-Kutt method
    auto * phi_rk = new double[2*N + 1];
    for (uint64_t i = 0; i < 2*N + 1; i++)
        phi_rk[i] = phi(a + ((b - a) / (2. * N) * i));

    double p = 1;
    double rk_err[N+1];
    double ** tmp_euler_res = compute_Euler(a, b, y_a, 2 * N, k, phi_rk);
    tmp_euler_res[0] = (double *) realloc(tmp_euler_res[0], sizeof(double)*(2*N+1));
    tmp_euler_res[1] = (double *) realloc(tmp_euler_res[1], sizeof(double)*(2*N+1));
    tmp_euler_res[2] = (double *) realloc(tmp_euler_res[2], sizeof(double)*(2*N+1));

    for (uint64_t i = 0; i < N+1; i++)
        rk_err[i] = pow(2, p)/(pow(2, p) - 1) * (tmp_euler_res[1][i * 2] - y_1[i]);
    double global_rk_err = *max_element(rk_err, rk_err + N + 1);


    // Write the data down!
    struct data_entry{
        string method_name;
        double **res;
        uint8_t res_dim;
        uint64_t N;
    };

    data_entry data_set[6] = {{"euler",                  euler_res, 3, N},
                              {"adams",                  adams_res, 3, N},
                              {"runge_kutt",             r_k_res, 3, N},
                              {"adams_bashfort_moulten", a_b_m_res, 3, N},
                              {"gir",                    gir_res, 3, N},
                              {"function",               reference_phi, 2, N}};

    for (const auto& data: data_set)
        write_down_results(N, data.res, data.res_dim, data.method_name);


    //Plotting
//    system("plotting_script.sh");

    return 0;
}
