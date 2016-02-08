#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>

double boxmuller()
{
double U1, U2, R, Y, Z1, Z2;

for(int i=0; i<100000; i++){

U1=((double)rand()/(double)RAND_MAX);
U2=((double)rand()/(double)RAND_MAX);

U1=2*U1-1;
U2=2*U2-1;

R=U1*U1+U2*U2;

if(R<=1){
Y=sqrt((-2*log(R))/R);
Z1=U1*Y;
Z2=U2*Y;

}
}
//std::cout <<"Z1="<<Z1<< std::endl;
return Z1;
}



int main(){
srand((unsigned)time(NULL));

double X0=100;
double T = 1;
double m =1000;
double delta_t=T/m;
double v_0, V_0, Z, r=0.03, delta=0.05, sigma=0.40, Xi, Xj, w, wdenominator, v_sum, Path_estimator_iterations=500;
double strike=100;

//std::vector< std::vector<double> > X;
std::vector <double> X;

for(int i=0; i<m; i++){

        if(i==0){

                X.push_back(X0);
        
        }

        if(i>0){

                Z=boxmuller();
        //      std::cout<<Z<< std::endl;
        //        Rn=UniRandom(b);
        //      std::cout<<"Rn="<<Rn<<std::endl;
                Xi=X[i-1];
        //      std::cout<<"Xi="<<Xi<<"at t="<<i-1<<  <<std::endl;
                Xj=Xi * (exp ((r-delta-0.5*sigma*sigma)*delta_t + sigma*sqrt(delta_t)*Z));
        //      std::cout<<"Xi="<<Xi<<"at t="<<i-1<<"X(i+1)="<<Xj <<std::endl;
                X.push_back(Xj);
        
        }

}

for(int j=0; j<m; j++){
std::cout<<X[j]<<std::endl;
}

return 0;
}
