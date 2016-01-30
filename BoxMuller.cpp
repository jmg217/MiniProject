#include <iostream>
#include <cmath>
#include <ctime>
#include <math.h>
#include <stdlib.h>
#include <fstream>

int main(void){

srand((unsigned)time(NULL));

std::ofstream outFile("/Users/tomgeary/Desktop/miniproj/Code/BMnormal.txt");

double U1, U2, R, Y, Z1, Z2;


for(int i=0; i<200000; i++){

U1=((double)rand()/(double)RAND_MAX);
U2=((double)rand()/(double)RAND_MAX);

U1=2*U1-1;
U2=2*U2-1;

R=U1*U1+U2*U2;

if(R<=1){
Y=sqrt((-2*log(R))/R);
Z1=U1*Y;
Z2=U2*Y;

//std::cout<<U1<<"\t"<<U2<<"\t"<<Z1<<"\t"<<Z2<<std::endl;
outFile <<Z1<< std::endl;
//outFile <<Z2<< std::endl;

}
}
//change
outFile.close();
return(0);

}
