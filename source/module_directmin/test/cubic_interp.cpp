#include <iostream>
#include <cmath>
#include <algorithm>

// Function to calculate cubic interpolation
double cubicFunction(double a, double b, double c, double d, double x) {
    return a * x * x * x + b * x * x + c * x + d;
}

// Function to find the minimum point of a cubic interpolation
double minimumCubicInterpol(double x_0, double x_1, double f_0, double f_1, double df_0, double df_1) {
    if (x_0 > x_1) {
        std::swap(x_0, x_1);
        std::swap(f_0, f_1);
        std::swap(df_0, df_1);
    }

    double r = x_1 - x_0;
    double a = -2.0 * (f_1 - f_0) / (r * r * r) + (df_1 + df_0) / (r * r);
    double b = 3.0 * (f_1 - f_0) / (r * r) - (df_1 + 2.0 * df_0) / r;
    double c = df_0;
    double d = f_0;
    double D = b * b - 3.0 * a * c;

    double x_min;
    if (D < 0.0) {
        if (f_0 < f_1) {
            x_min = x_0;
        } else {
            x_min = x_1;
        }
    } else {
        double sqrtD = std::sqrt(D);
        double r0 = (-b + sqrtD) / (3.0 * a) + x_0;
        if (x_0 < r0 && r0 < x_1) {
            double f_r0 = cubicFunction(a, b, c, d, r0 - x_0);
            if (f_0 > f_r0 && f_1 > f_r0) {
                x_min = r0;
            } else {
                if (f_0 < f_1) {
                    x_min = x_0;
                } else {
                    x_min = x_1;
                }
            }
        } else {
            if (f_0 < f_1) {
                x_min = x_0;
            } else {
                x_min = x_1;
            }
        }
    }

    return x_min;
}

// Test the function
int main() {
    double x_0 = 0.0, x_1 = 1.0;
    double f_0 = 0.0, f_1 = 1.0;
    double df_0 = -1.0, df_1 = 2.0;

    double x_min = minimumCubicInterpol(x_0, x_1, f_0, f_1, df_0, df_1);
    std::cout << "Minimum x: " << x_min << std::endl;

    return 0;
}
