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


double MeshEstimator(double strike, double r, double delta_t, int b, double m,  std::vector< std::vector<double> >& X, std::vector< std::vector< std::vector<double> > >& W, std::vector< std::vector<double> >& V, int num_assets){

double H; //payoff function
double C; //continuation value
double sum, V_0;

// temp vector in Estimator loop
std::vector< double > tempvec;


//Estimator loop

for(int i=0; i<m; i++){
tempvec.clear();

	for(int j=0; j<b; j++){
		if(i==0){
		H=payoff(X[(m-1)-i][j], strike)*exp(-r*delta_t*(m-i));
		tempvec.push_back(H);
		}
	
		else{
//continue to develope continuation vale			
			sum=0;
			for(int k=0; k<b; k++){
			sum+=(W[(m-i)][k][j])*V[i-1][k]; //m-i when i=1 is 10-1=9.when i=9 m-i=1. we get V_0 separately by using W[0][k][j]	
		
			}

			C=(1/((double)b))*sum;
			H=payoff(X[(m-1)-i][j], strike)*exp(-r*delta_t*(m-i));
			//H=0;
			if(H>=C){
			tempvec.push_back(H);
			//std::cout<<"H>=C and X="<<X[(m-1)-i][j]<<std::endl;
			}

			else{
			tempvec.push_back(C);
			//std::cout<<"C>H and X="<<X[(m-1)-i][j]<<std::endl;
			}	
		}
	}
V.push_back(tempvec);
}



sum=0;
for(int k=0; k<b; k++){
sum+=V[m-1][k];        
}

V_0=(1/((double)b))*sum;

/*
for ( std::vector<std::vector<double> >::size_type l = 0; l < V.size(); l++ )
{
   for ( std::vector<double>::size_type k = 0; k < V[l].size(); k++ )
   {
      std::cout << V[l][k] << ' ';
   }
   std::cout << std::endl;
}
*/

/*
std::cout<<"sum="<<sum<<std::endl;

std::cout<<"exp="<<exp(-r*delta_t)<<std::endl;
std::cout<<"r="<<r<<std::endl;
std::cout<<"delta_t="<<delta_t<<std::endl;
std::cout<<"1/b="<<1/b<<std::endl;
*/


/*
std::cout<<"V[9][0]="<<V[9][0]<<std::endl;
std::cout<<"V[9][1]="<<V[9][1]<<std::endl;
std::cout<<"V[9][2]="<<V[9][2]<<std::endl;
*/
//std::cout<<"V_0="<<V_0<<std::endl;

return V_0;
}
