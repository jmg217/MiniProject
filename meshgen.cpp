#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>

#define PI 3.14159265358979323846

double round( double value )
  {
  return floor( value + 0.5 );
  }

double boxmuller(){


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
std::cout <<"Z1="<<Z1<< std::endl;
return Z1;

}


double density(double Xold, double  Xnew, double sigma, double r, double delta){

double f, x;

x=(1/sigma)*(log(Xnew)-log(Xold)-(r-delta-0.5*sigma*sigma));

f= (1/(sigma*Xnew))*(1/(2*PI))*exp(-0.5*x*x);

return f;
}

int main(){
//1 year in 0.1 increments starting at 0 and ending at 1
srand((unsigned)time(NULL));
double X0=100;
double T = 1;
double m =10;
double t=T/m;
int b=3;
double Z, r=0.03, delta=0.05, sigma=0.4, Xi, Xj, w, wdenominator;

//MESH
std::vector< std::vector<double> > X;
//WEIGHTS for step 1 and beyond
std::vector< std::vector< std::vector<double> > > W;

for(int i=0; i<m; i++){


	std::vector< double > myvector;


	if(i==0){

	for(int l=0; l<b; l++){
		Z=boxmuller();
		//std::cout<<"Z="<<Z<< std::endl;
		Xi=X0 * (exp (r-delta-0.5*sigma*sigma + sigma*Z));
		myvector.push_back(Xi);	
	}
	}

	if(i>0){
	
	for(int j=0; j<b; j++){
		Z=boxmuller();
		Xi=X[i-1][j];
		Xj=Xi * (exp (r-delta-0.5*sigma*sigma + sigma*Z));
		myvector.push_back(Xj);
	}
	}

X.push_back(myvector);

}



for(int i=0; i<m; i++){

std::vector< std::vector<double> > dim2temp;
	
	if(i==0){
		for(int k=0; k<b; k++){
        	std:: vector<double> dim1temp;
                        //W[i][0][k]=w;
			//w=density(X0, X[i][k], sigma, r, delta);
			w=1;
		dim1temp.push_back(w);
		dim2temp.push_back(dim1temp);
		}
	}


	if(i>0){
	
		for(int k=0; k<b; k++){	
		std:: vector<double> dim1temp;
		wdenominator=0;
	
			for(int j=0; j<b; j++){
			w=density(X[i-1][j], X[i][k], sigma, r, delta);	
			dim1temp.push_back(w);
			}
			
			//this loop is for getting the denominator
			for(int l=0; l<b; l++){
			wdenominator+=dim1temp[l];	
			}
			
			//devide each element by the denominator
			for(int t=0; t<b; t++){
			dim1temp[t]=dim1temp[t]/wdenominator;
			}		
		dim2temp.push_back(dim1temp); //dim1 is full therefore we add it onto dim2 vector
		}	
	
	}

W.push_back(dim2temp);
}

//std::cout<< X.size() <<"\t"<< X[1].size() << std::endl;

//make this a head file for printing matrices

for ( std::vector<std::vector<double> >::size_type l = 0; l < X.size(); l++ )
{
   for ( std::vector<double>::size_type k = 0; k < X[l].size(); k++ )
   {
      std::cout << X[l][k] << ' ';
   }
   std::cout << std::endl;
}



}


