#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double rosenbrock(double x1, double x2)
{
    return (100 * (x2 - x1 * x1) * (x2 - x1 * x1) + (1 - x1) * (1 - x1));
}

int main()
{
    // Initial position
    double a_n[2] = {-1.5,-1};
    // sigma
    double sigma = 0.05;
    // rho
    double rho = 0.1;

    double ak;
    double a_new[2];
    double grad_x1 = (-400 * a_n[0] * a_n[1] + 400 * a_n[0] * a_n[0] * a_n[0] + 2 * a_n[0] - 2);
    double grad_x2 = (200 * a_n[1] - 200 * a_n[0] * a_n[0]);
    double fk;
    unsigned short int stop = 0;
    unsigned long int iter = 0;
    while (stop == 0) 
    {
        iter += 1;
        ak = 1;
        a_new[0] = a_n[0] - ak * grad_x1;
        a_new[1] = a_n[1] - ak * grad_x2;
        fk = rosenbrock(a_n[0], a_n[1]);
        while (rosenbrock(a_new[0], a_new[1]) > (fk + sigma * ak * (grad_x1 * (-grad_x1) + grad_x2 * (-grad_x2))))                                  
        {
            ak = ak * rho;
            a_new[0] = a_n[0] - ak * grad_x1;
            a_new[1] = a_n[1] - ak * grad_x2;
        }
        a_n[0] = a_new[0];
        a_n[1] = a_new[1];
        grad_x1 = (-400 * a_n[0] * a_n[1] + 400 * a_n[0] * a_n[0] * a_n[0] + 2 * a_n[0] - 2);
        grad_x2 = (200 * a_n[1] - 200 * a_n[0] * a_n[0]);
        if ( (fabs(grad_x1) < 1e-7) && (fabs(grad_x2) < 1e-7))
        {
            break;
        }    
    }
    printf("Minimum in the function found in %lu iterations in: %lf, %lf \n", iter, a_n[0], a_n[1]);
    return 0;
}