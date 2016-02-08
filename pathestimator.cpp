#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>



double density(double Xold, double  Xnew, double sigma, double r, double delta, double delta_t);

double boxmuller();

double Payoff(double x, double k){
double h;
h=x-k;
        if(h<0){
        h=0;
        }

return h;

}

double Product(double S, int j, int i, std::vector< std::vector<double> >& X, int b){

double Product, prod=1;

for (int t=0; t<b; t++){
	if(t<j){
	prod*=(S-X[i][t])/(X[i][j]-X[i][t]);

	}

	if(t>j){
	prod*=(S-X[i][t])/(X[i][j]-X[i][t]);
	}
}

return prod;
}


double interpolate_weights(int i, int k, double S, std::vector< std::vector<double> >& X, std::vector< std::vector< std::vector<double> > >& W, int b, double m){
double w, weightmatrix, Prod;
double sum=0;

for(int j=0; j<b; j++){
weightmatrix=W[i+1][k][j];
Prod=Product(S, j, i, X, b);
//std::cout<<weightmatrix<<"\t"<<Prod<<std::endl;
sum+=Prod*weightmatrix;
//sum+=W[i][k][j]*Product(S, j, i, X, b);
}

w=sum;

return w;
}


double PathEstimator(double strike, double r, double delta_t, int b, double m, double sigma, double delta, double X0, std::vector< std::vector<double> >& X, std::vector< std::vector< std::vector<double> > >& W, std::vector< std::vector<double> >& V){

//srand((unsigned)time(NULL));
double v_0, S_i, Z, C, H, sum, weight, w_s, sum_Z;

//Simulated path for sub optimal stoppage
std::vector<double > S;
//temp vector of S weights
std::vector<double > tempvec;
//weights matrix for S
std::vector< std::vector<double> > S_Weights; 

for(int i=0; i<m; i++){
	if(i==0){
	sum_Z=0;
	for(int z=0; z<7; z++){
                Z=boxmuller();
//		std::cout<<Z<<std::endl;
                sum_Z+=Z;
        }
        Z=(sum_Z)/7;
	//Z=boxmuller();
	S_i=X0 * (exp((r-delta-0.5*sigma*sigma)*delta_t + sigma*sqrt(delta_t)*Z));
	}

	else{
	sum_Z=0;
	for(int u=0; u<7; u++){
                Z=boxmuller();
                sum_Z+=Z;
        }
                Z=(sum_Z)/7;
	//Z=boxmuller();
	S_i=S[i-1]*(exp((r-delta-0.5*sigma*sigma)*delta_t + sigma*sqrt(delta_t)*Z));
	}

S.push_back(S_i);
}

//S weights loop
for(int t=0; t<(m-1); t++){
tempvec.clear();
	for(int h=0; h<b; h++){   //h=k
	sum=0;
	w_s=density(S[t], X[t+1][h], sigma, r, delta, delta_t);
		for(int g=0; g<b; g++){   //g=l
		sum+=(1/((double)b))*density(X[t][g], X[t+1][h], sigma, r, delta, delta_t);
		}
	//std::cout<<w_s<<"\t"<<sum<<std::endl;
	w_s=w_s/sum;	
	tempvec.push_back(w_s);
	}

S_Weights.push_back(tempvec);
}



////////////////////////////
// NEED TO INCLUDE DISCOUNT FACTOR IN CALCULATIONS. 
///////////////////////////
double con_val;
//sub optimal path loop
for(int i=0; i<m; i++){
sum=0;

	if(i==m-1){
	C=0;
	}
	
	else{
		for(int k=0; k<b; k++){	
			
		//weight=interpolate_weights(i, k, S[i], X, W, b, m);
		weight=S_Weights[i][k];
		con_val=V[(m-1)-i-1][k];
		//std::cout<<weight<<"\t"<<con_val<<std::endl;
		sum+=weight*con_val;			
		}
	C=(1/(double)b)*sum;
	}	
	

H=Payoff(S[i], strike)*exp(-r*delta_t*((i+1)));
//std::cout<<i<<"\t"<<H<<"\t"<<C<<std::endl;
//NOT SURE IF THIS IS > OR >=. CHECK THIS!
	if(H>=C || i==m-1){
	v_0=H; 
	break;
	}

}


/*
for(int i=0; i<m+1; i++){
if(i==0){
std::cout<<X0<<std::endl;
//std::cout<< v_0<<std::endl;
}
else{
std::cout<<S[i-1]<<std::endl;
}
}
*/
return v_0;

}
