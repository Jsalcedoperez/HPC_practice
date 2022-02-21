#include <stdio.h>
#include <stdlib.h>

/* func.c */

float f(float);
float Trap(float, float, int, float);

int main(){

float a = 1.0;
float b = 6.0;
int n = 1024;
float h = (b-a) / n;

float res;

res = Trap(a,b,n,h);

printf("the integral is, according to Simpson's rule, approx: %0.6f\n",res);

return 0;

}


float f(float x)
{

return x*x;

}


float Trap(float a, float b, int n, float h){

float integral;
float x;
int i;

integral = (f(a) + f(b))/2.0;
x = a;
for (i=1; i<= n-1; i++)
{
x = x + h;
integral = integral + f(x);
}
return integral*h;
}
