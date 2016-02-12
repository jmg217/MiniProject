#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>

void print_high_payoff(int b, double m, std::vector< std::vector<double> >& X, std::vector< std::vector<double> >& V){

std::ofstream outFile("/Users/tomgeary/Desktop/miniproj/Code/highpayoff.txt");

for (int t=0; t<m; t++){

	for (int i=0; i<b; i++){
	outFile << m-t <<"\t"<< X[m-1-t][i] <<"\t"<< V[t][i]<< std::endl;
	}
}

outFile.close();

}


