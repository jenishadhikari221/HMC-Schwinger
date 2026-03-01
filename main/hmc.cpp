#include "../header/hmc.hpp"
#include<random>
#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>
#include<functional>
#include <highfive/H5Easy.hpp>

using H5Easy::File ;
using complex = std::complex<double>;




int main(int argc, char *argv[]){
    size_t seed = 1234567890122;
    std::mt19937 gen(seed);
    int L = 5;
    double m0 = 1.0;

    std::string file_path = "../data/";
    std::string filename;
   int configs = 20; int  MD_steps = 25; int MD_trajectory_length = 1; double beta = 1.0;

    
   GetUserParam(argc,argv,filename,MD_trajectory_length,MD_steps,configs,L);
   
    filename = file_path + filename ; 
    File file(filename, File::ReadWrite|File::Truncate);

    std::vector<std::vector<std::vector<complex>>> lattice(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    
    creat_lattice(lattice,L,gen);
    
  // for(int i{0}; i < 10; i++){ 
    HMC(lattice,configs,MD_steps,MD_trajectory_length,gen,beta,m0,file);
   //check_leapfrog(argc,argv,lattice,gen);
   //check_hermitian(lattice,gen);
   
  // }
    
    return 0;
}