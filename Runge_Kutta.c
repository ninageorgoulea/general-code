// C program to implement Runge Kutta 4th order method
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

// Driver method
int main()
{
	clock_t start, end;
    double cpu_time_used;
	int i;
	i = 0;
	int n = 20000;
	double t[n], x[n], u[n];
	double q1 = 1.60217*pow(10,-19);
    double q2 = 1.60217*pow(10,-19);
    double k = 8.987551*pow(10,9);
    double m = 6.63*pow(10,-26);
	t[0] = 0.0;
	x[0] = 2.0*pow(10,-6);
	u[0] = -2.1983*pow(10,4);
	double h = 1.0*pow(10,-14);
	double k11, k12, k21, k22, k31, k32, k41, k42;
	FILE* fp;
	fp=fopen("rk4_1.txt","w");
	start = clock();
		while(t[i]<0.0000000002)
		{
		//Print values of x, u and t
        printf("x = %0.12lf, u = %0.12lf, t = %0.14lf\n", x[i], u[i], t[i]);
        fprintf(fp,"x = %0.12lf, u = %0.12lf, t = %0.14lf\n", x[i], u[i], t[i]);

        // Apply Runge Kutta Formulas to find next value of x and u
        k11 = h*u[i];
        k12 = h*(q1*q2*k)/(m*pow(x[i],2));

		k21 = h*(u[i]+0.5*k12);
		k22 = h*(q1*q2*k)/(m*pow(x[i]+0.5*k11,2));

		k31 = h*(u[i]+0.5*k22);
		k32 = h*(q1*q2*k)/(m*pow(x[i]+0.5*k21,2));

        k41 = h*(u[i]+k32);
        k42 = h*(q1*q2*k)/(m*pow(x[i]+k31,2));

        // Update next value of x, u and t
        x[i+1] = x[i] + (1.0/6.0)*(k11 + 2.0*k21 + 2.0*k31 + k41);
        u[i+1] = u[i] + (1.0/6.0)*(k12 + 2.0*k22 + 2.0*k32 + k42);
        t[i+1] = t[i] + h;
        i=i+1;
		}
	end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("The time is %lf", cpu_time_used);
	fclose(fp);
    return 0;
}

