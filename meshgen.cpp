#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>

#define PI 3.14159265358979323846

double MeshEstimator(double strike, double r, double delta_t, int b, double m,  std::vector< std::vector<double> >& X, std::vector< std::vector< std::vector<double> > >& W, std::vector< std::vector<double> >& V, int num_assets);

double PathEstimator(double strike, double r, double delta_t, int b, double m, double sigma, double delta, double X0, std::vector< std::vector<double> >& X, std::vector< std::vector< std::vector<double> > >& W, std::vector< std::vector<double> >& V, int num_assets);

void print_high_payoff(int b, double m, std::vector< std::vector<double> >& X, std::vector< std::vector<double> >& V);

double round( double value )
{
  return floor( value + 0.5 );
}



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

int UniRandom(int b){
double rn;
rn=((double)rand()/(double)RAND_MAX)*((double)(b-1));

rn=round(rn);

return rn;
}

double density(double Xold, double  Xnew, double sigma, double r, double delta, double delta_t){

double f, x;

x=(1/(sigma*sqrt(delta_t)))*(log(Xnew)-log(Xold)-(r-delta-0.5*sigma*sigma)*delta_t);

f= (1/(sigma*sqrt(delta_t)*Xnew))*(1/(2*PI))*exp(-0.5*x*x);

return f;
}


int main(){
//1 year in 0.1 increments starting at 0 and ending at 1
/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!NOTE: may need to include delta t's in calculation of the mesh since delta t is 0.1 and not 1.!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
srand((unsigned)time(NULL));
double X0=100;
double T = 1;
double m = 10;
double delta_t=T/m;
int b=500, Rn, N=20;
double v_0, V_0, Z, r=0.03, delta=0, sigma=0.4, Xi, Xj, w, wdenominator, v_sum, Path_estimator_iterations=500, vtotal_sum=0, Vtotal_sum=0;
double strike=100, sum_Z=0;
int num_assets=1;

//MESH
std::vector< std::vector<double> > X;
//WEIGHTS for step 1 and beyond
std::vector< std::vector< std::vector<double> > > W;
//temp vector in MeshGen for loop
std::vector< double > myvector;
//2 d temp vector in WeightsGen for loop
std::vector< std::vector<double> > dim2temp;
//1 d vector in Weightsgen for loop
std:: vector<double> dim1temp;
//mesh estimator high bias 2-d matrix
std::vector< std::vector<double> > V;
//V values
std::vector< double > Vvector;
//v values
std::vector< double > vvector;

for(int iterator=0; iterator<N; iterator++){
X.clear();
W.clear();
V.clear();
//MeshGen for loop
for(int i=0; i<m; i++){


	myvector.clear();

	if(i==0){

	for(int l=0; l<b; l++){
		sum_Z=0;
		for(int z=0; z<num_assets; z++){
		Z=boxmuller();
		sum_Z+=Z;
		}
		Z=(sum_Z)/((double)num_assets);
	//	std::cout<<Z<< std::endl;
		Xi=X0 * (exp ((r-delta-0.5*sigma*sigma)*delta_t + sigma*sqrt(delta_t)*Z));
		myvector.push_back(Xi);	
	}
	}

	if(i>0){
	
	for(int j=0; j<b; j++){
		sum_Z=0;
		for(int u=0; u<num_assets; u++){
                Z=boxmuller();
                sum_Z+=Z;
                }
                Z=(sum_Z)/((double)num_assets);
	//	std::cout<<Z<< std::endl;
		Rn=UniRandom(b);
	//	std::cout<<"Rn="<<Rn<<std::endl;
		Xi=X[i-1][Rn];
	//	std::cout<<"Xi="<<Xi<<"at t="<<i-1<<  <<std::endl;
		Xj=Xi * (exp ((r-delta-0.5*sigma*sigma)*delta_t + sigma*sqrt(delta_t)*Z));
	//	std::cout<<"Xi="<<Xi<<"at t="<<i-1<<"X(i+1)="<<Xj <<std::endl;
		myvector.push_back(Xj);
	}
	}

X.push_back(myvector);

}


//WeightsGen for loop 
//NOTE: W^i_(j,k) IS REPRESENTED AT W[i][k][j] where k is at step i+1 and j is at step i.
for(int i=0; i<m; i++){

dim2temp.clear();
	
	if(i==0){
		for(int k=0; k<b; k++){
        	dim1temp.clear();
                        //W[i][0][k]=w;
			//w=density(X0, X[i][k], sigma, r, delta);
			w=1;
		dim1temp.push_back(w);
		dim2temp.push_back(dim1temp);
		}
	}


	if(i>0){
	
		for(int k=0; k<b; k++){	
		dim1temp.clear();
		wdenominator=0;
	
			for(int j=0; j<b; j++){
			w=density(X[i-1][j], X[i][k], sigma, r, delta, delta_t);//step 1 in X is X[0] step 1 in W is W[1]	
			dim1temp.push_back(w);
			wdenominator+=w;
			}
			
			//this loop is for getting the denominator
			/*for(int l=0; l<b; l++){
			wdenominator+=dim1temp[l];	
			//std::cout<<dim1temp[l]<<std::endl;
			}
			*/

			//wdenominator=(1/((double)b))*wdenominator;		
			//std::cout<<wdenominator<<std::endl;

			//devide each element by the denominator
			for(int t=0; t<b; t++){
			dim1temp[t]=(((double)b)*dim1temp[t])/wdenominator;
			}		
		dim2temp.push_back(dim1temp); //dim1 is full therefore we add it onto dim2 vector
		}	
	
	}

W.push_back(dim2temp);
}

double check=0;
//check all the weights from X0 are 1
for(int e=0; e<b; e++){
if(W[0][e][0]!=1){
std::cout<<"there is an error with the weights. check that W[0][k][0]'s =1"<<std::endl;
}
}
//check that the weights going into a node sum to 1
for(int q=1; q<m; q++){ 
for(int a=0; a<b; a++){
check=0;
for(int E=0; E<b; E++){
check+=W[q][a][E];
//std::cout<<W[1][0][E]<<std::endl;
}/*
if(check != b){
std::cout<<"the sum of the weights does not add up. Ignore if check="<<check<<"=b="<<b<<std::endl;
}*/
}
}
V_0=MeshEstimator(strike, r, delta_t, b, m, X, W, V, num_assets);
Vvector.push_back(V_0);
Vtotal_sum+=V_0;

std::cout<<"V_0="<<V_0<<std::endl;

v_sum=0;
for(int f=0; f<Path_estimator_iterations; f++){
v_sum+=PathEstimator(strike, r, delta_t, b,  m, sigma, delta, X0, X, W, V, num_assets);
}

v_0=(1/double(Path_estimator_iterations))*v_sum;
vvector.push_back(v_0);
vtotal_sum+=v_0;

std::cout<<"v_0="<<v_0<<std::endl;
}//this is the end of the loop over the whole process.


V_0=(1/double(N))*Vtotal_sum;
v_0=(1/double(N))*vtotal_sum;

double std_div_V=0, std_div_v=0, squaresumV=0, squaresumv=0, Verror=0, verror=0;

for(int h=0; h<N; h++){
squaresumV+=(Vvector[h]-V_0)*(Vvector[h]-V_0);
squaresumv+=(vvector[h]-v_0)*(vvector[h]-v_0);
}
std_div_V=sqrt((1/double(N))*squaresumV);
std_div_v=sqrt((1/double(N))*squaresumv);

Verror=1.96*std_div_V*(1/sqrt(double(N)));
verror=1.96*std_div_v*(1/sqrt(double(N)));
std::cout<<"V_0="<<V_0<<"\t"<<"V error="<<Verror<<std::endl;
std::cout<<"v_0="<<v_0<<"\t"<<"v error="<<verror<<std::endl;


std::ofstream outFile("/Users/tomgeary/Desktop/miniproj/Code/results.txt");

outFile << N <<"\t"<< b <<"\t"<< Path_estimator_iterations<<"\t"<< V_0 <<"\t"<< v_0 <<"\t"<< Verror+V_0 <<"\t"<< verror+v_0 << std::endl;

outFile.close();

print_high_payoff( b, m, X, V);

//I used the following to check if the pathestimator code was producing random stock prices. I did this by printing out all the values from the pathestimator code.
//PathEstimator( strike, r, delta_t, b, m, sigma, delta, X0);

//std::cout<< X.size() <<"\t"<< X[1].size() << std::endl;

//make this a head file for printing matrices

/*
for ( std::vector<std::vector<double> >::size_type l = 0; l < X.size(); l++ )
{
   for ( std::vector<double>::size_type k = 0; k < X[l].size(); k++ )
   {
      std::cout << X[l][k] << ' ';
   }
   std::cout << std::endl;
}
*/
/*
int i, j, k;

std::cout<< "type in i, k and j"<<std::endl;

std::cin>> i;

std::cin>> k ;

//std::cin>> j;
double sum=0;
for(int d=0; d<b; d++ ){
std::cout<< W[i][k][d]<<std::endl;
sum+=d;
}

std::cout<<"sum="<<sum<<std::endl;
*/

return 0;


}


