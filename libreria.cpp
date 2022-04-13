#include "libreria.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TApplication.h>
#include <TGraphErrors.h>
#include <fstream>

using namespace std;



void readfile(string nomefile, int& n, vector<string>& at,string& a1,string& a2, double& b1, double& b2, vector<double>& X, vector<double>& Y, vector<double>& Z){

  ifstream file(nomefile);
  if(!file.good()){
    cout << "Impossibile leggere il file " << endl;
  }
  file >>n>> a1 >> a2 >> b1 >> b2;
  double x,y,z;
  string Ar;
  for(int i=0; i<n; i++){
    file >> Ar >> x >> y >> z;
    at.push_back(Ar);
    X.push_back(x*pow(10,-10));    //converto i dati in sistema internazionale
    Y.push_back(y*pow(10,-10));
    Z.push_back(z*pow(10,-10));
  }
  file.close();
};


void writefile (ofstream& file, int& n, vector<string>& at, string& a1, string& a2, double& b1, double& b2, vector<double>& X, vector<double>& Y, vector<double>& Z){

  file << n << endl;
  file.precision(15);
  file << a1 << '\t'  <<    a2 << '\t' <<  b1 << '\t'  << b2 << endl;
  file.precision(10);
								   
  for(int i=0; i<n; i++){
    file << at[i] << '\t'   << X[i]*pow(10,10) << '\t' <<  Y[i]*pow(10,10) << '\t'  << Z[i]*pow(10,10) << endl;     //ho riconvertito i dati nel file di output in armstrong
  }
};


double energia_pot(int& n,double& epsilon, double& sigma,double& kb, vector<double>& X, vector<double>& Y, vector<double>& Z){
  double U=0;
  double r_ij;

  for(int i=0; i<n-1; i++){
    for(int j=i+1; j<n; j++){
      r_ij = sqrt( pow(X[i] - X[j],2) +  pow(Y[i] - Y[j],2) +  pow(Z[i] - Z[j],2)  );
      U = U + 4*epsilon * ( pow(sigma / r_ij ,12) - pow(sigma / r_ij,6)  );
    }
  }
  return U;
};


double energia_cin(int& n, double& m, vector<double>& VX, vector<double>& VY, vector<double>& VZ){
  double K=0;
  for(int i=0; i<n; i++){
    K += 0.5*m*pow(VX[i],2) + 0.5*m*pow(VY[i],2) + 0.5*m*pow(VZ[i],2);
  }
  return K;
};



void forza(int& n, double& epsilon, double& sigma, vector<double>& X, vector<double>& Y, vector<double>& Z, vector<double> &fx, vector<double> &fy, vector<double> &fz){
  double r_ij;
  fx.clear();
  fy.clear();
  fz.clear();
  vector<double> fx_temp(n,0.0);
  vector<double> fy_temp(n,0.0);
  vector<double> fz_temp(n,0.0);
  double a,b,c;
  for(int i=0;i<n-1;i++){
    a=0; b=0; c=0;
    for(int j=i+1; j<n; j++){
      r_ij = sqrt( pow(X[i] - X[j],2) +  pow(Y[i] - Y[j],2) +  pow(Z[i] - Z[j],2)  );
      a = (4*epsilon*(- (6*pow(sigma,6))/(pow(r_ij,6)) + (12*pow(sigma,12))/(pow(r_ij,12)) )*(X[i]-X[j])/(r_ij*r_ij) );
      b = (4*epsilon*(- (6*pow(sigma,6))/(pow(r_ij,6)) + (12*pow(sigma,12))/(pow(r_ij,12)) )*(Y[i]-Y[j])/(r_ij*r_ij) );
      c = (4*epsilon*( -(6*pow(sigma,6))/(pow(r_ij,6)) + (12*pow(sigma,12))/(pow(r_ij,12)) )*(Z[i]-Z[j])/(r_ij*r_ij) );

      fx_temp[i] += a;
      fx_temp[j] += -a;
      fy_temp[i] += b;
      fy_temp[j] += -b;
      fz_temp[i] += c;
      fz_temp[j] += -c;
    }
  }
  fx = fx_temp;
  fy = fy_temp;
  fz = fz_temp;
};
