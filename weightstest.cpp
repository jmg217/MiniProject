#include <vector>
#include<iostream>

int main(){

double b=3, m=10, w, wdenominator;

std::vector< std::vector< std::vector<double> > > W;
std::vector< std::vector<double> > dim2temp;
std:: vector<double> dim1temp;
for(int i=0; i<m; i++){

dim2temp.clear();

        if(i==0){
                for(int k=0; k<b; k++){
                dim1temp.clear();
                        //W[i][0][k]=w;
                        //w=density(X0, X[i][k], sigma, r, delta);
                        w=i;
                dim1temp.push_back(w);
                dim2temp.push_back(dim1temp);
                }
        }


        if(i>0){

                for(int k=0; k<b; k++){
                dim1temp.clear();
                wdenominator=0;

                        for(int j=0; j<b; j++){
                        w=i;
                        dim1temp.push_back(w);
                        }

                        //this loop is for getting the denominator
                       /* for(int l=0; l<b; l++){
                        wdenominator+=dim1temp[l];
                        }
			

                        //devide each element by the denominator
                        for(int t=0; t<b; t++){
                        dim1temp[t]=dim1temp[t]/wdenominator;
                        }
		*/	
                dim2temp.push_back(dim1temp); //dim1 is full therefore we add it onto dim2 vector
                }

        }
W.push_back(dim2temp);
}

/*
std:: vector<double> dim1temp;
std::vector< std::vector<double> > dim2temp;

for(int i=0; i<3; i++){
dim2temp.clear();
dim1temp.clear();
w=1*(i+1);
dim1temp.push_back(w);
w=2*(i+1);
dim1temp.push_back(w);
dim2temp.push_back(dim1temp);
dim1temp.clear();
w=3*(i+1);
dim1temp.push_back(w);
w=4*(i+1);
dim1temp.push_back(w);
dim2temp.push_back(dim1temp);
W.push_back(dim2temp);
}

int i, j, k;
std::cout<< "type in i, j, k"<<std::endl;

std::cin>> i;

std::cin>> j ;

std::cin>> k;

std::cout<< W[i][j][k]<<std::endl;
*/

int i, j, k;

std::cout<< "type in i, j and k"<<std::endl;

std::cin>> i;

std::cin>> j ;
std::cin>> k ;
std::cout<< W[i][j][k]<<std::endl;



return 0;
}
