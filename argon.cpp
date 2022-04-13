#include <iostream>
#include <cmath>
#include <vector>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TApplication.h>
#include <TGraphErrors.h>
#include <fstream>
#include <TH1D.h>
#include <bits/stdc++.h>
#include "libreria.h"


using namespace std;


namespace dati{
  double kb = 1.38065e-23;   //in Joule/K
  double sigma = 3.4e-10;    //nel SI
  double epsilon = 119* kb;  //epsilon in Joule
  string a1,a2;
  double b1,b2;
  double T=20;    //in kelvin
  double tau = 5e-15;     //5 femptosecondi
  double m = 6.68e-26;    //è la massa in chili

}






//MAIN

int main(){

  
  TApplication app("app",0,NULL);

  vector<double> rx;
  vector<double> ry;
  vector<double> rz;
  vector<string> at;

  int N;

  readfile("initial_position.xyz", N, at,dati::a1, dati::a2,dati::b1,dati::b2, rx, ry,rz);
  const char *path="/home/delle/Computazionale/Quencing/Dati/output.xyz";
  ofstream file(path); 
  writefile(file,N,at,dati::a1, dati::a2, dati::b1, dati::b2,rx,ry,rz);

 
  
  double U;
  U = energia_pot(N,dati::epsilon, dati::sigma, dati::kb, rx, ry, rz);  
  cout.precision(15);
 


  //Generatore
  
  TRandom3 rnd;
  rnd.SetSeed(time(0));
  vector<double> vx_star;
  vector<double> vy_star;
  vector<double> vz_star;
  double v_cm_x, v_cm_y, v_cm_z;
  for(int i=0; i<N; i++){
    
    vx_star.push_back(  sqrt(dati::kb*dati::T/dati::m)*rnd.Gaus(0,1));
    vy_star.push_back(  sqrt(dati::kb*dati::T/dati::m)*rnd.Gaus(0,1));
    vz_star.push_back(  sqrt(dati::kb*dati::T/dati::m)*rnd.Gaus(0,1));
    v_cm_x += vx_star[i];
    v_cm_y += vy_star[i];
    v_cm_z += vz_star[i];
  }
  v_cm_x = v_cm_x/N;
  v_cm_y = v_cm_y/N;
  v_cm_z = v_cm_z/N;

  vector<double> vx;
  vector<double> vy;
  vector<double> vz;
  for(int i=0; i<N; i++){
    vx.push_back(vx_star[i] - v_cm_x  );    //queste sono le velocità traslate
    vy.push_back(vy_star[i] - v_cm_y  );    //quelle che lui chiama con v barra
    vz.push_back(vz_star[i] - v_cm_z  );
  }

  double K_bar;
    K_bar  = energia_cin(N, dati::m, vx,vy,vz);

  double K_0;
  K_0 = 3./2*N*dati::kb*dati::T ;
 

  for (int i=0; i<N; i++) {
    vx[i] = vx[i]*sqrt(K_0/K_bar);
    vy[i] = vy[i]*sqrt(K_0/K_bar);
    vz[i] = vz[i]*sqrt(K_0/K_bar);
  }


  //Calcolo della forza

  vector<double> fx;
  vector<double> fy;
  vector<double> fz;


  

  //VELOCITY VERLET

  double t=0;
  double tmax = 3e-11;
  double dt = dati::tau;
  double ax, ay, az;
  vector<double> ax_temp,  ay_temp,  az_temp;
  TGraph g1;
  TGraph g2, g3, g4;
  int n_passi=0;
  int n_campioni=0;
  double E,K;
  double T_temp;
  double E_med = 0;
  double std_dev = 0;


 
  forza(N,dati::epsilon, dati:: sigma, rx,ry,rz,fx,fy,fz);      //prima del loop altrimenti le calcolo troppe volte


  while(t<tmax){

    E=0;
    K=0;
    ax_temp.clear();
    ay_temp.clear();
    az_temp.clear();
    v_cm_x = 0;
    for(int i=0; i<N; i++){
      ax = fx[i]/dati::m;     // questa è a(t) calcolata con f(t) per la i esima particella
      ax_temp.push_back(ax);
      ay = fy[i]/dati::m;     // questa è a(t) calcolata con f(t) per la i esima particella
      ay_temp.push_back(ay);
      az = fz[i]/dati::m;     // questa è a(t) calcolata con f(t) per la i esima particella
      az_temp.push_back(az);
      rx[i] = rx[i] + vx[i]*dt + 0.5 * ax*dt*dt;    //calcolo nuove posizioni
      ry[i] = ry[i] + vy[i]*dt + 0.5 * ay*dt*dt;
      rz[i] = rz[i] + vz[i]*dt + 0.5 * az*dt*dt;
    }
    
    forza(N,dati::epsilon, dati:: sigma, rx,ry,rz,fx,fy,fz);
    U = energia_pot(N,dati::epsilon, dati::sigma, dati::kb, rx, ry, rz);
	
    for(int i=0; i<N; i++){
      
      if(fx[i]*vx[i] + fy[i]*vy[i] + fz[i]*vz[i] >= 0){
	ax = fx[i]/dati::m ;   //sono le nuove accelerazioni
	ay = fy[i]/dati::m ;
	az = fz[i]/dati::m ; 
	vx[i] = vx[i] + 0.5*( ax_temp[i] + ax )*dt;       //calcolo nuove velocità
	vy[i] = vy[i] + 0.5*( ay_temp[i] + ay )*dt;
	vz[i] = vz[i] + 0.5*( az_temp[i] + az )*dt;
      }

      if(fx[i]*vx[i] + fy[i]*vy[i] + fz[i]*vz[i] < 0){
	vx[i] = 0;
	vy[i] = 0;
	vz[i] = 0;
      }
    }
  
    K = energia_cin(N, dati::m, vx,vy,vz);
    E = K+U;
    T_temp = 2./3 * K /(N*dati::kb);
    E_med += E;
    g1.SetPoint(n_campioni,t,E);
    g2.SetPoint(n_campioni,t,U);
    g3.SetPoint(n_campioni,t,K);
    g4.SetPoint(n_campioni,t,T_temp);
    n_campioni += 1;

    if(n_passi%100==0){
      writefile(file,N,at,dati::a1, dati::a2, dati::b1, dati::b2,rx,ry,rz);
    }

    string name="/home/delle/Computazionale/Quencing/Dati/output_" + std::to_string(n_passi+1) + ".xyz"; 
    ofstream file1(name);
    writefile(file1,N,at,dati::a1, dati::a2, dati::b1, dati::b2,rx,ry,rz);
    file1.close();
      
    n_passi += 1;
    t += dt;
  }
  n_passi -= 1;
  //fine ciclo while

  
  file.close();
  
  E_med = E_med /n_passi;
  cout << "energia media = " << E_med << endl;
  for(int i=0; i<n_campioni; i++){
    std_dev += pow( g1.GetPointY(i) - E_med  ,2);
  }
  std_dev = sqrt(std_dev);
  cout << "std_dev = " << std_dev << endl;
  cout << "numero di passi = " << n_passi << endl;



  //GRAFICI CONTROLLO
    TCanvas c("","",800,600);
   c.Divide(2,2);

    c.cd(1);
  // TCanvas c1;
  g1.SetTitle("Energia");
  g1.Draw("ap");
    c.cd(2);
  //TCanvas c2;
  g2.SetTitle("Energia potenziale");
  g2.Draw("ap");
   c.cd(3);
  // TCanvas c3;
  g3.SetTitle("Energia cinetica");
  g3.Draw("ap");
   c.cd(4);
  //TCanvas c4;
  g4.SetTitle("Temperatura");
  g4.Draw("ap");






  app.Run(true);

  return 0;

}
