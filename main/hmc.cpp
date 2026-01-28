#include "../header/hmc.hpp"
#include<random>
#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>
#include<functional>

using complex = std::complex<double>;

std::string filename;

int main(int argc, char *argv[]){
    size_t seed = 1234567;
    std::mt19937 gen(seed);
    int L = 3;
    double m0 = 1;

    std::vector<std::vector<std::vector<complex>>> lattice(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    
    creat_lattice(lattice,L,gen);
    
    std::vector<std::vector<std::vector<complex>>> test_chi(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    std::vector<std::vector<std::vector<complex>>> test_phi(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));

    generate_Phi(test_chi,gen);
    test_phi = f_M(test_chi,lattice,m0);
    complex test;
    complex test_2;

    //std::vector<std::vector<std::vector<complex>>> inverse = cg(f_M_f_Mdag,test_phi,lattice,m0,100,1e-10);
    
/*/ check if M and M_dagger diagonal element match


    std::vector<std::vector<std::vector<complex>>> v(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2,complex(1.0,0.0))));
    std::vector<std::vector<std::vector<complex>>> ferm ;
    std::vector<std::vector<std::vector<complex>>> fermd;
    for (int t = 0; t < L; t++) {      
        for (int x = 0; x < L; x++) {
            for(int m = 0;m<2;m++){
                ferm = f_M(v,lattice,m0);
                fermd = f_Mdag(v,lattice,m0);
                
                std::cout<<":t"<<t<<",x "<<x<<"  "<<ferm[t][x][m]<< " "<<fermd[t][x][m]<<" \n";
            }
        }
        
    } /**/  
   int configs = 20; int MD_steps = 20; int MD_trajectory_length = 1; double beta = 1;
   HMC(lattice,configs,MD_steps,MD_trajectory_length,gen,beta,m0);
    
    
    return 0;
}