#include <stdio.h>
#include <stdlib.h>
#include <math.h>


typedef struct exx{
  int mu;
} est;
  
double func(int x, void *params){
  est estructure = *(est*)params;
  estructure.mu = 20;
  int b = estructure.mu;
  return x*b;
}

int main(void){
  int a = 1;
  est parameter;
  parameter.mu =50;
  printf("%f \n",func(2.0,&parameter));
  return 0;
}



