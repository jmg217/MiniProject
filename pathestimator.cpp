#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>


double boxmuller();

double Payoff(double x, double k){
double h;
h=k-x;
        if(h<0){
        h=0;
        }

return h;

}


void PathEstimator(double strike, double r, double delta_t, int b, double m, double sigma, double delta, double X0, std::vector< std::vector<double> >& X, std::vector< std::vector< std::vector<double> > >& W){

//srand((unsigned)time(NULL));
double v_0, S_i, Z, C, H;

std::vector<double > S;

for(int i=0; i<m; i++){
	if(i==0){
	Z=boxmuller();
	S_i=X0 * (exp ((r-delta-0.5*sigma*sigma)*delta_t + sigma*sqrt(delta_t)*Z));
	}

	else{
	Z=boxmuller();
	S_i=S[i-1]*(exp ((r-delta-0.5*sigma*sigma)*delta_t + sigma*sqrt(delta_t)*Z));
	}

S.push_back(S_i);
}

//sub optimal path loop
for(int i=0; i<m; i++){

	for(int k=0; k<b; k++){
		
		if(i==0){

		}		

		weight=interpolate_weights();
		C=(1/(double)b)*weight*
	

		if(i==m){
		C=0;
		}
	}

}


/*
for(int i=0; i<m; i++){
std::cout<<S[i]<<std::endl;
}
*/

//return v_0;

}
