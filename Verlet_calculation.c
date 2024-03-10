// C program to implement Verlet Algorithm
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define n 40000000
// n = the elements of the relevant vector
#define k 8.987551E9
// k = Coulomb constant in N*m^2*C^(-2)
double q1 = 1.60217E-19;
// q1 = the charge of the moving ion Ar or electron in C
double me = 9.10938356E-31;
// me = the electron mass in kg
double mar = 6.63E-26;
// mar = the ion Ar mass in kg
double xc[500] = {0.0};
// xc[] = the relevant vector that corresponds to the equidistant spaces of the surface PMMA in the x-direction in m
double yc[500] = {0.0};
// yc[] = the relevant vector that corresponds to the equidistant spaces of the surface PMMA in the y-direction in m
double q2[500] = {0.0};
// q2 = the charge vector which corresponds to the relevant equidistant spaces of the PMMA surface in C
double t[n];
// t[] = the array of time in s
double x[n];
// x[] = the array of the position of the moving ion Ar or electron in the x-direction in m
double ux[n];
// ux[] = the array of the velocity of the moving ion Ar or electron in the x-direction in m/s
double y[n];
// y[] = the array of the position of the moving ion Ar or electron in the y-direction in m
double uy[n];
// x[] = the array of the velocity of the moving ion Ar or electron in the y-direction in m/s
double h = 1E-12;
// h = the time step in s
int i;
// i = the relevant component of each vector
double L1 = 2.75E-6;
// L1 = the length of the PMMA surface

// Fx = the electrostatic force of the moving ion Ar or electron with the PMMA surface along the x-direction in N
static double Fx(double x1, double x2, double x3, double x4, double x5, double x6, double q1_, double q2_)
{
	double r;
	r = (q1_*q2_*k*(x1))/pow((pow(x1,2)+pow(x2,2)),1.5)+(q1_*q2_*k*(x3))/pow((pow(x3,2)+pow(x4,2)),1.5)+(q1_*q2_*k*(x5))/pow((pow(x5,2)+pow(x6,2)),1.5);
	return (r);
}

// Fy = the electrostatic force of the moving ion Ar or electron with the PMMA surface along the y-direction in N
static double Fy(double x1, double x2, double x3, double x4, double x5, double x6, double q1_, double q2_)
{
	double r;
	r = (q1_*q2_*k*(x2))/pow((pow(x1,2)+pow(x2,2)),1.5)+(q1_*q2_*k*(x4))/pow((pow(x3,2)+pow(x4,2)),1.5)+(q1_*q2_*k*(x6))/pow((pow(x5,2)+pow(x6,2)),1.5);
	return (r);
}

// sumax = the acceleration of the moving ion Ar or electron along the x-direction in m/s^2
static double sumax(double x1, double x2, double m)
{
	int j;
	double sum = 0.0;
	for(j=0;j<500;j++)
	{
		sum+=Fx(x1-xc[j], x2-yc[j], x1-xc[j]-L1, x2-yc[j], x1-xc[j]+L1, x2-yc[j], q1, q2[j])/m;
	}
	return (sum);
}

// sumay = the acceleration of the moving ion Ar or electron along the y-direction in m/s^2
static double sumay(double x1, double x2, double m)
{
	int j;
	double sum = 0.0;
	for(j=0;j<500;j++)
	{
		sum+=Fy(x1-xc[j], x2-yc[j], x1-xc[j]-L1, x2-yc[j], x1-xc[j]+L1, x2-yc[j], q1, q2[j])/m;
	}
	return (sum);
}

// E = the function that calculates total energy of the moving ion Ar or electron in eV
static double E(double u1, double u2, double x1, double x2, double m)
{
	double r, q;
	int j;
	double sum = 0.0;
	r = (m*(pow(u1,2)+pow(u2,2)))/2.0;
	for(j=0;j<500;j++)
	{
		sum+=(k*q1*q2[j])*(1/sqrt(pow(x1-xc[j],2)+pow(x2-yc[j],2))+1/sqrt(pow(x1-xc[j]-L1,2)+pow(x2-yc[j],2))+1/sqrt(pow(x1-xc[j]+L1,2)+pow(x2-yc[j],2)));
	}
	q = (r + sum)*(6.24150913*10E18);
	return (q);
}

// randfrom() = fuction that selects random values from the interval (min, max)
static double randfrom(double min, double max)
{
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand()/div);
}

// Driver method
int main()
{
	int j, j_, s, loop, stop;
	// j = the equidistant number of spaces that the PMMA surface is consisted of
	// j_ = the number of the ions Ar or electrons that interact with the PMMA surface
	// stop = a parameter that determines whether to stop the calculations or not
	// loop corresponds to the relevant numbers of the components of the q2[] vector
	stop = 2;
	xc[0] = 6E-9;
	for(j=0;j<500;j++)
	{
		xc[j+1]=xc[j]+12E-9;
	}
	clock_t start, end;
    double cpu_time_used;
    // cpu_time_used = the time that is being spent in order to be computed the relevant calculations in s

	i = 0;
	FILE* fp;
	fp = fopen("verlet_program.txt","w");
	start = clock();
	for(j_=1;j_<10;j_++)
	{
		printf("The number of throw is j = %d\n",j_);
		fprintf(fp,"The number of throw is j = %d\n",j_);
		i = 0;
		stop = 2;
		t[0] = 0.0;
		x[0] = randfrom(0, 2.75E-6);
		ux[0] = randfrom(-170.0, 170.0);
		y[0] = 2.2E-6;
		uy[0] = randfrom(-170, 0.0);
		xc[0] = 6E-9;
		while(stop!=1 && y[i]<=2.2E-6)
		{

		// Print values of x, y, ux, uy, t and E
		if(i%100==0.0)
		{
        printf(" %1.12lf        %1.12lf       %1.12lf        %1.12lf        %1.10lf       %1.10lf       %1.12lf     %1.21lf\n", x[i], y[i], ux[i], uy[i], sumax(x[i], y[i], mar), sumay(x[i], y[i], mar), t[i], E(ux[i], uy[i], x[i], y[i], mar));
        fprintf(fp," %1.12lf         %1.12lf        %1.12lf         %1.12lf       %1.10lf      %1.10lf      %1.12lf       %1.21lf\n", x[i], y[i], ux[i], uy[i], sumax(x[i], y[i], mar), sumay(x[i], y[i], mar), t[i], E(ux[i], uy[i], x[i], y[i], mar));
    	}

        // Implement x, y, ux and uy
        x[i+1] = x[i] + h*ux[i] + (pow(h,2)*sumax(x[i], y[i], mar))/(2.0);
        y[i+1] = y[i] + h*uy[i] + (pow(h,2)*sumay(x[i], y[i], mar))/(2.0);

        ux[i+1] = ux[i] + (h*(sumax(x[i], y[i], mar)+sumax(x[i+1], y[i+1], mar)))/(2.0);
        uy[i+1] = uy[i] + (h*(sumay(x[i], y[i], mar)+sumay(x[i+1], y[i+1], mar)))/(2.0);

        t[i+1] = t[i] + h;
        i = i + 1;

        if(x[i]<0.0)
		{
		x[i]=x[i]+2.75E-6;
		printf("Left Boundary\n");
		fprintf(fp, "Left Boundary\n");
		}

		else if(x[i]>2.75E-6)
		{
		x[i]=x[i]-2.75E-6;
		printf("Right Boundary\n");
		fprintf(fp, "Right Boundary\n");
		}

		else if(y[i]<-1.6E-6 && x[i]>=1.125E-6 && x[i]<=1.625E-6)
		{
		s = (x[i])/(12E-9);
		q2[s]+=1.60217E-19;
		printf("Trench\n");
		fprintf(fp, "Trench\n");
		stop = 1;
		}

		else if(y[i]<0.0 && (x[i]<1.125E-6 || x[i]>1.625E-6))
		{
		s = (x[i])/(12E-9);
		q2[s]+=1.60217E-19;
		printf("Out of trench\n");
		fprintf(fp, "Out of trench\n");
		stop = 1;
		}
		}
	}
	end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("The time is %lf for %d throws\n", cpu_time_used, j_-1);
    fprintf(fp, "The time is %lf for %d throws\n", cpu_time_used, j_-1);
	for(loop=0; loop<500; loop++)
	{
    	printf("For the loop = %d the charge is q2 = %0.25lf\n", loop, q2[loop]);
    	fprintf(fp,"For the loop = %d the charge is q2[loop] = %0.25lf\n", loop, q2[loop]);
	}
	fclose(fp);
    return 0;
}

