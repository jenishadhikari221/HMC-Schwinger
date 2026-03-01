#pragma once
#include<random>
#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>
#include<complex>
#include<functional>
#include <highfive/H5Easy.hpp>

using H5Easy::File ;

using complex = std::complex<double>;

void creat_lattice(std::vector<std::vector<std::vector<complex>>> &lattice,int L, std::mt19937 &gen);
complex plaquette(const std::vector<std::vector<std::vector<complex>>> &lattice,std::vector<int> n,int mu,int nu);
double g_plaquette(const std::vector<std::vector<std::vector<complex>>> &lattice,const double beta);
void generate_Phi(std::vector<std::vector<std::vector<complex>>> & Phi,std::mt19937 &gen);
double delta_f(int n,int m);     
double delta_f(std::vector<int> n,std::vector<int> m);                 
complex fermion(const std::vector<std::vector<std::vector<complex>>> & U_gauge, std::vector<int> n,std::vector<int> n_prime,int alpha,int beta,double m0);
complex fermion_dagger(const std::vector<std::vector<std::vector<complex>>> & U_gauge, std::vector<int> n,std::vector<int> n_prime,int alpha,int beta,double m0);
std::vector<std::vector<std::vector<complex>>> f_M(const std::vector<std::vector<std::vector<complex>>> & Phi,const std::vector<std::vector<std::vector<complex>>> & U_gauge,const double m0);
std::vector<std::vector<std::vector<complex>>> f_Mdag(const std::vector<std::vector<std::vector<complex>>> & Phi,const std::vector<std::vector<std::vector<complex>>> & U_gauge,const double m0);
std::vector<std::vector<std::vector<complex>>> f_Mdag_f_M(const std::vector<std::vector<std::vector<complex>>> & Phi,const std::vector<std::vector<std::vector<complex>>> & U_gauge,const double m0);
double normsquared(const std::vector<std::vector<std::vector<complex>>> & psi);
double scalar_product(std::vector<std::vector<std::vector<complex>>> p, std::vector<std::vector<std::vector<complex>>> t);
void assign_add_mul(std::vector<std::vector<std::vector<complex>>> &x, std::vector<std::vector<std::vector<complex>>> &p,const double alpha);
void assign_mul_add(std::vector<std::vector<std::vector<complex>>> &p, std::vector<std::vector<std::vector<complex>>> &r,const double beta);
std::vector<std::vector<std::vector<complex>>> cg(std::vector<std::vector<std::vector<complex>>> & chi,
    std::vector<std::vector<std::vector<complex>>> & psi, const std::vector<std::vector<std::vector<complex>>> & U_gauge, double m0, size_t Max_Iterations, double epsilon);
double schwinger_action(std::vector<std::vector<std::vector<complex>>> & psi, const std::vector<std::vector<std::vector<complex>>> & U_gauge,double beta,double m0);


std::vector<std::vector<std::vector<complex>>> g_force(const std::vector<std::vector<std::vector<complex>>> & U_gauge,double beta);
double hamiltonian(std::vector<std::vector<std::vector<complex>>> & pi, std::vector<std::vector<std::vector<complex>>> & phi, const std::vector<std::vector<std::vector<complex>>> & U_gauge,double beta,double m0);
std::vector<std::vector<std::vector<complex>>> net_force(std::vector<std::vector<std::vector<complex>>> & psi, const std::vector<std::vector<std::vector<complex>>> & U_gauge,double beta,double m0);
void leapfrog(double beta, double m0,std::vector<std::vector<std::vector<complex>>> phi,std::vector<std::vector<std::vector<complex>>>& init_gauge,std::vector<std::vector<std::vector<complex>>>& init_pi
    ,int MD_steps,int MD_trajectory_length, int MD_direction);

void HMC(std::vector<std::vector<std::vector<complex>>>& U_gauge, int configs, int MD_steps,int MD_trajectory_length,std::mt19937 & gen, double beta,double m0, File & file );




double temporal_bc(int t, int L, int step);
void GetUserParam( int argc, char *argv[], std::string & filename,int & epsilon, int & tau, int & config, int & L);
void check_leapfrog( int argc, char *argv[], std::vector<std::vector<std::vector<complex>>> & U_gauge,std::mt19937 & gen);
complex scalar_product_complex(std::vector<std::vector<std::vector<complex>>> p, std::vector<std::vector<std::vector<complex>>> t);
void check_hermitian(const  std::vector<std::vector<std::vector<complex>>> & lattice,std::mt19937 & gen);
double topological_charge(const std::vector<std::vector<std::vector<complex>>>& U_gauge);