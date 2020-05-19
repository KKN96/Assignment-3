#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#define num 128
int PI=3.1415926535;

double f(double x)
{ 
	if(x==0.0)
		return 1.0;
	else
		return sin(x)/x;
}
double analy(double x)
{ 
	if (x<=1&&x>=-1)
		return sqrt(PI/2.0);
	else
		return 0.0;
}





int main()
{ 
	int i;
	double x_min=-30,x_max=30,dx;
	double ft[num],k[num],tru[num];
	fftw_complex in[num], out[num];
	fftw_plan p;
	p=fftw_plan_dft_1d(num,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	
	dx=(x_max-x_min)/(double)(num-1);

	for (i=0;i<num;i++)
	{	if(i<num/2){
			k[i]=PI*2.0*(double)i/(double)num/dx;}
		else
			{k[i]=PI*2.0*(double)(i-num)/(double)(num)/dx;}
		in[i][0]=f(x_min+dx*(double)i);
		in[i][1]=0.0;}

	FILE*file;
	file = fopen("output.txt","w");
	fftw_execute(p); 

	for(i=0;i<num;i++) 
	{ 
		ft[i]=dx*sqrt(1.0/(2.0*PI))*(cos(-x_min*k[i])*out[i][0]-sin(-x_min*k[i])*out[i][1]);
		tru[i]=analy(k[i]);
	}
	for (i=num/2;i<num;i++){
		fprintf(file,"%f, %f, %f\n", k[i],ft[i],tru[i]);}
		
	for (i=0;i<num/2;i++)
		{
		fprintf(file,"%f, %f, %f\n", k[i],ft[i],tru[i]);}
		
	fclose(file);
	fftw_destroy_plan(p);
	return 0;
}


