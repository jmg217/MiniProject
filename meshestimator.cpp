#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>


double payoff(double x, double k){

double h;
h=k-x;
	if(h<0){
	h=0;
	}

return h;
}

double contin_val(){

}

int main(){


double H; //payoff function
double C; //continuation value
double strike=100;

//mesh estimator(high bias)
std::vector< std::vector<double> > V;
// temp vector in Estimator loop
std::vector< double > tempvec;


//Estimator loop
for(int i=m; i=0; i--){
tempvec.clear();
	for(int j=0; j<b; j++){
		if(i==m){
		H=payoff(X[i][j], stirke);
		tempvec.push_back(H);
		}

		else{
		
//continue to develope continuation vale
			for(int l; l<b; l++){
			

			}

		H=payoff(X[i][j], strike);
		

			if(H>C){

			}

			else{
		
			}


		}
	}

V.push_back(tempvec);
}

return 0;
}
