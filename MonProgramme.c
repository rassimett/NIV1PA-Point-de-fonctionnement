#include <stdio.h>
#include <math.h>

const float n= 0.68, E= 1.0, R= 100.0, Is= 1e-15, V0= 0.025; 

// Coding the function: f(U) = E - U - R * Is * (exp(U*(n/V0)) - 1) 
float f(float U) {
    float expo = expf(U * (n / V0));
    return E - U - R * Is * (expo - 1);
}

// Coding the derivative of the function: f'(U) = -1 - (R * Is * (n / V0)) * expo
float df(float U) {
    float expo = expf(U * (n / V0));
    return -1 - (R * Is * (n / V0)) * expo;
}

// Newton method. Returns root estimate. iterations is output. If it fails, returns NAN.
float newton(float U0, int max_iter, float eps, int *iterations) {
    float U = U0;
    for (int i = 0; i < max_iter; i++) {
        float fU = f(U);
        float dfU = df(U);
        if (fabs(dfU) <  1e-12) {
            *iterations = i + 1;
            return NAN; 
        }
        float U_next = U - fU / dfU;
        if (fabs(U_next - U) < eps) {
            *iterations = i + 1;
            return U_next; // Converged
        }
        U = U_next;
    }
    *iterations = max_iter;
    return U;
}

// Bisection method. Returns root estimate or NAN.
float bisection(float a, float b, float eps, int max_iter, int *iterations) {
    float fa = f(a);
    float fb = f(b);
    if (fa * fb > 0) {
        *iterations = 0;
        return NAN; // No root in [a, b]
    }
    for (int i = 0; i < max_iter; i++) {
        float c = (a + b) / 2;
        float fc = f(c);
        if ((b - a) / 2 < eps) {
            *iterations = i + 1;
            return c; // Converged
        }
        if (fa * fc <= 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    *iterations = max_iter;
    return (a + b) / 2; // Return midpoint as best estimate
}

// Writing my results
void write_iv_file(const char *filename){
    FILE *fp = fopen(filename, "w");
    if(!fp) {
        perror("fopen");
        return;
    }
    fprintf(fp, "# U[V]    I_diode[A]     I_generator[A]\n");
    for(int i=0; i<=100; i++ ){
        float U = i * 0.1;
        float Id = Is * (expf(U * (n / V0)) - 1);
        float Ig = (E - U) / R;
        fprintf(fp, "%.2f %.12e %.12e\n", U, Id, Ig);
    }
    fclose(fp);
}

int main(void){
    

}

