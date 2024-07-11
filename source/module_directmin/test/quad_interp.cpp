#include <iostream>
#include <cmath>
#include <algorithm>

// Function to calculate quadratic interpolation value
double quadraticInterpol(double x_0, double x_1, double f_0, double f_1, double df_0, double df_1, double x) {
    if (x_0 > x_1) {
        std::swap(x_0, x_1);
        std::swap(f_0, f_1);
        std::swap(df_0, df_1);
    }

    double r = x_1 - x_0;
    double a = (df_1 - df_0) / r;
    double b = df_0 - a * x_0;
    double c = f_0 - df_0 * x_0 + 0.5 * a * x_0 * x_0;

    double interp_value = a * x * x + b * x + c;
    return interp_value;
}

// Test function
int main() {
    double x_0 = 0.0, x_1 = 1.0;
    double f_0 = 0.0, f_1 = 1.0;
    double df_0 = 1.0, df_1 = 2.0;
    double x = 0.5;

    double interp_value = quadraticInterpol(x_0, x_1, f_0, f_1, df_0, df_1, x);
    std::cout << "Quadratic interpolation value at x = " << x << ": " << interp_value << std::endl;

    return 0;
}
