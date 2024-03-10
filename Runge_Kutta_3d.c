// C program to implement Runge Kutta 4th order method
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

// A sample differential equation "du/dt=(q/m)*E"
float dudt(float tf, float xf, float uxf, float yf, float uyf, float m, float k, float q1, float q2)
{
    return (q1*q2*k)/(m*(pow(xf,2)+pow(yf,2)));
}
int n=100000;

// Driver method
int main()
{
	clock_t start, end;
    double cpu_time_used;
	int i;
	i=0;
	float t[n], x[n], ux[n], y[n], uy[n];
	t[0]=0.0;
	x[0]=1.0;
	ux[0]=-80.0;
	y[0]=1.0;
	uy[0]=-80.0;
	float h=0.0001;
	float k11x, k12x, k21x, k22x, k31x, k32x, k41x, k42x, k11y, k12y, k21y, k22y, k31y, k32y, k41y, k42y;
	FILE* fp;
	fp=fopen("rk4_2d.txt","w");
	fprintf(fp, "   x            ux           y           uy             t\n");
	start = clock();
	while(x[i]>2*h && y[i]>2*h)
	{
		//Print values of x, ux, y, uy and t
        printf(" x = %0.6f, ux = %0.6f, y = %0.6f, uy = %0.6f t = %0.3f\n", x[i], ux[i], y[i], uy[i], t[i]);
        fprintf(fp,"%0.6f         %0.6f     %0.6f       %0.6f         %0.3f\n", x[i], ux[i], y[i], uy[i], t[i]);

        // Apply Runge Kutta Formulas to find next value of x, ux, y and uy
        k11x = h*ux[i];
        k12x = h*dudt(t[i], x[i], ux[i], y[i], uy[i], 1.0, 1.0, 1.0, 1.0);

        k11y = h*uy[i];
        k12y = h*dudt(t[i], x[i], ux[i], y[i], uy[i], 1.0, 1.0, 1.0, 1.0);

		k21x = h*(ux[i]+0.5*k12x);
		k22x = h*dudt(t[i] + 0.5*h, x[i] + 0.5*k11x, ux[i] + 0.5*k12x,  y[i] + 0.5*k11y, uy[i] + 0.5*k12y, 1.0, 1.0, 1.0, 1.0);

		k21y = h*(uy[i]+0.5*k12y);
		k22y = h*dudt(t[i] + 0.5*h, x[i] + 0.5*k11x, ux[i] + 0.5*k12x, y[i] + 0.5*k11y, uy[i] + 0.5*k12y, 1.0, 1.0, 1.0, 1.0);

		k31x = h*(ux[i]+0.5*k22x);
		k32x = h*dudt(t[i] + 0.5*h, x[i] + 0.5*k21x, ux[i] + 0.5*k22x, y[i] + 0.5*k21y, uy[i] + 0.5*k22y, 1.0, 1.0, 1.0, 1.0);

		k31y = h*(uy[i]+0.5*k22y);
		k32y = h*dudt(t[i] + 0.5*h, x[i] + 0.5*k21x, ux[i] + 0.5*k22x, y[i] + 0.5*k21y, uy[i] + 0.5*k22y, 1.0, 1.0, 1.0, 1.0);

        k41x = h*(ux[i]+k32x);
        k42x = h*dudt(t[i] + h, x[i] + k31x, ux[i] + k32x, y[i] + k31y, uy[i] + k32y, 1.0, 1.0, 1.0, 1.0);

        k41y = h*(uy[i]+k32y);
        k42y = h*dudt(t[i] + h, x[i] + k31x, ux[i] + k32x, y[i] + k31y, uy[i] + k32y, 1.0, 1.0, 1.0, 1.0);

        // Update next value of x, ux, y, uy and t
        x[i+1] = x[i] + (1.0/6.0)*(k11x + 2.0*k21x + 2.0*k31x + k41x);
        ux[i+1] = ux[i] + (1.0/6.0)*(k12x + 2.0*k22x + 2.0*k32x + k42x);
        y[i+1] = y[i] + (1.0/6.0)*(k11y + 2.0*k21y + 2.0*k31y + k41y);
        uy[i+1] = uy[i] + (1.0/6.0)*(k12y + 2.0*k22y + 2.0*k32y + k42y);
        t[i+1]=t[i]+h;

        i=i+1;
	}
	end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("The time is %lf", cpu_time_used);
	fclose(fp);
    return 0;
}
