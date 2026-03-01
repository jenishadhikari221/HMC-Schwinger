#include "../header/hmc.hpp"
#include<random>
#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>
#include <cstring>
#include <complex>
#include <functional>
#include <iomanip>
#include <highfive/H5Easy.hpp>

using H5Easy::File;
using complex = std::complex<double>;

// U_mu(x)

void creat_lattice(std::vector<std::vector<std::vector<complex>>> &lattice,int L, std::mt19937 &gen){

    std::normal_distribution<double> dist(0.0,1.0);
    for (int x=0;x<L;x++){
        for(int y=0;y<L;y++){
            for(int z=0;z<2;z++){
                lattice[x][y][z] = std::exp(complex(0.0,1.0*dist(gen)));
            }
        }
    }
}

void generate_Phi(std::vector<std::vector<std::vector<complex>>> & Phi,std::mt19937 &gen){
    int L = Phi.size();
    for(int x = 0;x<L;x++){
        for(int y = 0;y<L;y++){
            for(int z = 0;z<2;z++){
                std::normal_distribution<double> dist(0.0,1.0);
                Phi[x][y][z] = complex(dist(gen),1.0*dist(gen));
            }
        }
    }

}

void generate_Pi(std::vector<std::vector<std::vector<complex>>> & Phi,std::mt19937 &gen){
    int L = Phi.size();
    for(int x = 0;x<L;x++){
        for(int y = 0;y<L;y++){
            for(int z = 0;z<2;z++){
                std::normal_distribution<double> dist(0.0,1.0);
                Phi[x][y][z] = complex(dist(gen),0.0);
            }
        }
    }

}

complex plaquette(const std::vector<std::vector<std::vector<complex>>> &lattice,std::vector<int> n,int mu,int nu){
    int x = n[0];
    int y = n[1];
    int L = lattice.size();
    complex p;
    p = lattice[x][y][mu]*lattice[(x+1)%L][y][nu]*std::conj(lattice[x][(y+1)%L][mu])*std::conj(lattice[x][y][nu]);
    
    return p;
    
}

double g_plaquette(const std::vector<std::vector<std::vector<complex>>> &lattice,const double beta){

    complex sum = 0;
    int N = lattice.size();
    std::vector<int> n;
    int mu = 0; int nu = 1;

    for (int x=0;x<N;x++){
        for(int y=0;y<N;y++){
            n.push_back(x);n.push_back(y);
            sum += plaquette(lattice,n,mu,nu);
            n.clear();
        }

    }
    sum = beta*((double) N * (double) N - sum.real());
    return sum.real();
}
double delta_f(int n,int m){
    if (n==m){
    return 1.0;}
    return 0.0;
}

double delta_f(std::vector<int> n,std::vector<int> m){
    int x = n[0];
    int y = n[1];
    int x_p = m[0];
    int y_p = m[1];
    if (x==x_p && y==y_p){
        return 1.0;}
    return 0.0;
}

double temporal_bc(int t, int L, int step) {
    // step = +1 (forward) or -1 (backward)
    if (step == +1 && t == L-1) return -1.0;
    if (step == -1 && t == 0)   return -1.0;
    return 1.0;
}
complex fermion(const std::vector<std::vector<std::vector<complex>>> & U_gauge, std::vector<int> n,std::vector<int> n_prime,int alpha,int beta,double m0){
    std::vector<std::vector<std::vector<complex>>> sigma(3,std::vector<std::vector<complex>>(2,std::vector<complex>(2)));

    sigma[0] ={{{0.0,0.0},{1.0,0.0}},{{1.0,0.0},{0.0,0.0}}};
    sigma[1] ={{{0.0,0.0},{0.0,-1.0}},{{0.0,1.0},{0.0,0.0}}};
    sigma[2] ={{{1.0,0.0},{0.0,0.0}},{{0.0,0.0},{-1.0,0.0}}};
    int L = U_gauge.size();
    complex out;
    int x = n[0];
    int y = n[1];
    int x_p = n_prime[0];
    int y_p = n_prime[1];
    int mu = 0; int nu =1;

    auto new_n_prime = n_prime;
    int t = n_prime[0];
    double bc = temporal_bc(t, L, +1);

    new_n_prime[0] = (t + 1) % L;

    out = complex((m0+2)*delta_f(alpha,beta)*delta_f(n_prime,n),0.0);
    new_n_prime[mu] = (n_prime[mu] + 1) % L;
    out -=bc*0.5*(complex(1.0,0.0)-sigma[mu][alpha][beta])*U_gauge[x_p][y_p][mu]*complex(delta_f(new_n_prime,n),0.0);
    new_n_prime = n_prime;
    new_n_prime[mu] = (n_prime[mu] - 1 + L) % L;
    bc = temporal_bc(t, L, -1);
    out -=bc*0.5*(complex(1.0,0.0)+sigma[mu][alpha][beta])*std::conj(U_gauge[(x_p-1+L)%L][y_p][mu])*complex(delta_f(new_n_prime,n),0.0);
    
    new_n_prime = n_prime;
    new_n_prime[nu] = (n_prime[nu] + 1) % L;
    out -=0.5*(complex(1.0,0.0)-sigma[nu][alpha][beta])*U_gauge[x_p][y_p][nu]*complex(delta_f(new_n_prime,n),0.0);
    new_n_prime = n_prime;
    new_n_prime[nu] = (n_prime[nu] - 1 + L) % L;
    out -=0.5*(complex(1.0,0.0)+sigma[nu][alpha][beta])*std::conj(U_gauge[(x_p)%L][(y_p-1+L)%L][nu])*complex(delta_f(new_n_prime,n),0.0);

    return out;

}

complex fermion_dagger(const std::vector<std::vector<std::vector<complex>>> & U_gauge, std::vector<int> n,std::vector<int> n_prime,int alpha,int beta,double m0){
    
    std::vector<std::vector<std::vector<complex>>> sigma(3,std::vector<std::vector<complex>>(2,std::vector<complex>(2)));
    sigma[0] ={{{0.0,0.0},{1.0,0.0}},{{1.0,0.0},{0.0,0.0}}};
    sigma[1] ={{{0.0,0.0},{0.0,-1.0}},{{0.0,1.0},{0.0,0.0}}};
    sigma[2] ={{{1.0,0.0},{0.0,0.0}},{{0.0,0.0},{-1.0,0.0}}};
    int L = U_gauge.size();
    complex out;
    int x = n[0];
    int y = n[1];
    int x_p = n_prime[0];
    int y_p = n_prime[1];
    int mu = 0; int nu =1;

    auto new_n_prime = n_prime;

    out = complex((m0+2)*delta_f(alpha,beta)*delta_f(n_prime,n),0.0);
    int t = n_prime[0];
    double bc = temporal_bc(t, L, +1);
    new_n_prime[mu] = (n_prime[mu] + 1) % L;
    out -=bc*0.5*(complex(1.0,0.0)+sigma[mu][alpha][beta])*(U_gauge[x_p][y_p][mu])*complex(delta_f(new_n_prime,n),0.0);
    new_n_prime = n_prime;
    new_n_prime[mu] = (n_prime[mu] - 1 + L) % L;
    bc = temporal_bc(t, L, -1);
    out -=bc*0.5*(complex(1.0,0.0)-sigma[mu][alpha][beta])*std::conj(U_gauge[(x_p-1+L)%L][y_p][mu])*complex(delta_f(new_n_prime,n),0.0);
    
    new_n_prime = n_prime;
    new_n_prime[nu] = (n_prime[nu] +1 ) % L;
    out -=0.5*(complex(1.0,0.0)+sigma[nu][alpha][beta])*(U_gauge[x_p][y_p][nu])*complex(delta_f(new_n_prime,n),0.0);
    new_n_prime = n_prime;
    new_n_prime[nu] = (n_prime[nu] - 1 + L) % L;
    out -=0.5*(complex(1.0,0.0)-sigma[nu][alpha][beta])*std::conj(U_gauge[(x_p)%L][(y_p-1+L)%L][nu])*complex(delta_f(new_n_prime,n),0.0);
    
    return out;

}



std::vector<std::vector<std::vector<complex>>> f_M(const std::vector<std::vector<std::vector<complex>>> & Phi,const std::vector<std::vector<std::vector<complex>>> & U_gauge,const double m0){
    int L = U_gauge.size();
    std::vector<std::vector<std::vector<complex>>> fermion_matrix(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    
    complex temp;
    std::vector<int> n;
    std::vector<int> n_prime;
    for(int x_p = 0; x_p<L;x_p++){
        for(int y_p = 0; y_p<L;y_p++){
            for(int alpha = 0; alpha<2;alpha++){
                temp = complex(0.0,0.0);
                for(int x = 0;x<L;x++){
                    for(int y = 0;y<L;y++){
                        for(int beta = 0;beta<2;beta++){
                            n_prime.push_back(x_p);
                            n_prime.push_back(y_p);
                            n.push_back(x);
                            n.push_back(y);

                            temp += fermion(U_gauge,n,n_prime,alpha,beta,m0)*Phi[x][y][beta];

                            n_prime.clear();n.clear();
                        }
                    }
                }
            fermion_matrix[x_p][y_p][alpha] = temp;
            }
        }
    }

    return fermion_matrix;

}
std::vector<std::vector<std::vector<complex>>> f_Mdag(const std::vector<std::vector<std::vector<complex>>> & Phi,const std::vector<std::vector<std::vector<complex>>> & U_gauge,const double m0){
    int L = U_gauge.size();
    std::vector<std::vector<std::vector<complex>>> fermion_matrix(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    
    complex temp;
    std::vector<int> n;
    std::vector<int> n_prime;
    for(int x_p = 0; x_p<L;x_p++){
        for(int y_p = 0; y_p<L;y_p++){
            for(int alpha = 0; alpha<2;alpha++){
                temp = complex(0.0,0.0);
                for(int x = 0;x<L;x++){
                    for(int y = 0;y<L;y++){
                        for(int beta = 0;beta<2;beta++){
                            n_prime.push_back(x_p);
                            n_prime.push_back(y_p);
                            n.push_back(x);
                            n.push_back(y);

                            temp += std::conj(fermion(U_gauge,n_prime,n,beta,alpha,m0))*Phi[x][y][beta];

                            n_prime.clear();n.clear();
                        }
                    }
                }
            fermion_matrix[x_p][y_p][alpha] = temp;
            }
        }
    }
    return fermion_matrix;
    
}

std::vector<std::vector<std::vector<complex>>> f_Mdag_f_M(const std::vector<std::vector<std::vector<complex>>> & Phi,const std::vector<std::vector<std::vector<complex>>> & U_gauge,const double m0){

    return f_Mdag(f_M(Phi,U_gauge,m0),U_gauge,m0);
}

double normsquared(const std::vector<std::vector<std::vector<complex>>> & psi){
    int L = psi.size();
    complex norm_x = complex(0.0,0.0);
    for(int x = 0;x<L;x++){
        for(int y = 0;y<L;y++){
                norm_x += (psi[x][y][0]*std::conj(psi[x][y][0]) + psi[x][y][1]*std::conj(psi[x][y][1]));
            }
        }
    return norm_x.real();
}

double scalar_product(std::vector<std::vector<std::vector<complex>>> p, std::vector<std::vector<std::vector<complex>>> t){
    int L = p.size();
    complex x_dir = complex(0.0,0.0);
    for(int x = 0;x<L;x++){
        for(int y = 0;y<L;y++){
            x_dir += (std::conj(p[x][y][0])*t[x][y][0]) + (std::conj(p[x][y][1])*t[x][y][1]);
        }
    }
    return x_dir.real();
}

void assign_add_mul(std::vector<std::vector<std::vector<complex>>> &x,std::vector<std::vector<std::vector<complex>>> &p,const double alpha){
    int L = x.size();
    for(int x_p = 0;x_p<L;x_p++){
        for(int y_p = 0;y_p<L;y_p++){
            x[x_p][y_p][0] = x[x_p][y_p][0] + alpha*p[x_p][y_p][0];
            x[x_p][y_p][1] = x[x_p][y_p][1] + alpha*p[x_p][y_p][1];
        }
    }

}

void assign_mul_add(std::vector<std::vector<std::vector<complex>>> &p, std::vector<std::vector<std::vector<complex>>> &r,const double beta){
    int L = p.size();
    for(int x_p = 0;x_p<L;x_p++){
        for(int y_p = 0;y_p<L;y_p++){ // p{i+1} = r{i+1} + beta{i+1} * p{i}
            p[x_p][y_p][0] = r[x_p][y_p][0] + beta*p[x_p][y_p][0];
            p[x_p][y_p][1] = r[x_p][y_p][1] + beta*p[x_p][y_p][1];
        }
    }

}


std::vector<std::vector<std::vector<complex>>> cg( std::vector<std::vector<std::vector<complex>>> & chi, std::vector<std::vector<std::vector<complex>>> & psi, const std::vector<std::vector<std::vector<complex>>> & U_gauge, double m0, size_t Max_Iterations, double epsilon){
	
    double rsqr, err;
    double alpha;
    double beta;
    double max_err;
    double temp_sc_product;
    int L = psi.size();

	// maximum norm squared of r for which the CG stops (remember stopping condition ||r{i}||/||psi|| < epsilon => (r{i},r{i}) < epsilon^2 * (psi,psi))
    max_err = normsquared(psi);
	max_err = std::pow(epsilon,2.0)*max_err;

	// vectors r,p,t,x have the same meaning as defined in lecture
    std::vector<std::vector<std::vector<complex>>> r(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    std::vector<std::vector<std::vector<complex>>> p(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    std::vector<std::vector<std::vector<complex>>> t(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    std::vector<std::vector<std::vector<complex>>> x(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2,complex(0.0,0.0))));  // initial guess for x, zero here
    

    
	//initial residual r{0} = psi - A * x{0} (x{0} = 0 here)
    r = psi;

	// p{0} = r{0}
    p = psi;

    // rsqr = (r{0},r{0})
    rsqr = normsquared(r);

    for(size_t i = 0; i < Max_Iterations; i++) {
    	// t = A * p
    	t = f_Mdag_f_M(p,U_gauge,m0);

    	//alpha{i} = (r{i},r{i})/(p{i},t{i}) (remember we stored (r{i},r{i}) in rsqr)
        temp_sc_product = scalar_product(p,t);
    	alpha = (rsqr/temp_sc_product);

    	// x{i+1} = x{i} + alpha{i} * p{i}
    	assign_add_mul(x,p,alpha);

    	// r{i+1} = r{i} - alpha{i} * t{i}
        alpha = -1.0*alpha;

    	assign_add_mul(r,t,alpha);

    	// err = (r{i+1},r{i+1})
    	err = normsquared(r);

		// std::cout<< " Max_ Error<< "<< max_err << "Error_x[" << i << "] = " << err << std::endl;

    	// check for convergence
    	if(err <= max_err) {
    		//std::cout << "Required Precision Reached: " << std::endl;
    		return x ;
    	}

    	// beta{i+1} = (r{i+1},r{i+1})/(r{i},r{i}) = err/rsqr
    	beta = (err/rsqr);
    	// p{i+1} = r{i+1} + beta{i+1} * p{i}
		assign_mul_add(p,r,beta);

    	// ensure that rsqr is (r{i},r{i}) for the next iteration (i.e. i -> i+1)
    	rsqr = err;
    }

    throw std::runtime_error("Error: CG did not converged: Error = " + std::to_string(rsqr));
}

double schwinger_action(std::vector<std::vector<std::vector<complex>>> & psi, const std::vector<std::vector<std::vector<complex>>> & U_gauge,double beta,double m0){

    double gauge_action = g_plaquette(U_gauge,beta);
    int L = psi.size();
    complex sum = complex(0.0,0.0);
    std::vector<std::vector<std::vector<complex>>> M_Mdag_in = cg(psi,psi,U_gauge,m0,1000,1e-15);
    for(int i = 0; i < L; i++){
        for(int j = 0; j < L; j++){
            for(int alpha = 0; alpha < 2; alpha++){
                    sum += std::conj(psi[i][j][alpha])*M_Mdag_in[i][j][alpha];
            }
        }
    }
    return gauge_action + sum.real();
}

double hamiltonian(std::vector<std::vector<std::vector<complex>>> & pi,
    std::vector<std::vector<std::vector<complex>>> & phi, const std::vector<std::vector<std::vector<complex>>> & U_gauge,double beta,double m0){
    double gauge_action = g_plaquette(U_gauge,beta);
    std::vector<std::vector<std::vector<complex>>> inverse_matrix = cg(phi,phi,U_gauge,m0,1000,1e-8);



    int L = phi.size();
    complex sum = complex(0.0,0.0);
    complex sum_2 = complex(0.0,0.0);
        for(int i = 0; i < L; i++){
            for(int j = 0; j < L; j++){
                sum += std::conj(pi[i][j][0])*pi[i][j][0] + std::conj(pi[i][j][1])*pi[i][j][1]; // Pi^2
                for (int alpha = 0; alpha<2;alpha ++){
                        sum_2 += std::conj(phi[i][j][alpha])*inverse_matrix[i][j][alpha]; // Phi^+ (MM^+)^-1 Phi
                    
                }
                
                //std::cout<<"Sum = "<<sum<<" sum_2 = "<<sum_2<<"\n";
                
            }
        }
        return gauge_action + 0.5*sum.real() + sum_2.real();

}

std::vector<std::vector<std::vector<complex>>> g_force(const std::vector<std::vector<std::vector<complex>>> & U_gauge,double beta){
    int L = U_gauge.size();
    std::vector<std::vector<std::vector<complex>>> K(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));

    
    for(int i = 0; i < L; i++){
        for(int j = 0; j < L; j++){
            
            K[i][j][0]  = -beta*(U_gauge[i][j][0]*(U_gauge[(i+1)%L][j][1]*std::conj(U_gauge[i][(j+1)%L][0])*std::conj(U_gauge[i][j][1])+
                     std::conj(U_gauge[(i+1)%L][(j-1+L)%L][1])*std::conj(U_gauge[i][(j-1+L)%L][0])*U_gauge[i][(j-1+L)%L][1])).imag();
                     
            K[i][j][1]  = -beta*(U_gauge[i][j][1]*(U_gauge[i][(j+1)%L][0]*std::conj(U_gauge[(i+1)%L][(j)%L][1])*std::conj(U_gauge[i][j][0])+
                     std::conj(U_gauge[(i-1+L)%L][(j+1)%L][0])*std::conj(U_gauge[(i-1+L)%L][(j)%L][1])*U_gauge[(i-1+L)%L][(j)%L][0])).imag();
        }
    }
            
    return K;
}


std::vector<std::vector<std::vector<complex>>> net_force(std::vector<std::vector<std::vector<complex>>> & phi, const std::vector<std::vector<std::vector<complex>>> & U_gauge,double BETA,double m0){
    std::vector<std::vector<std::vector<complex>>> sigma(3,std::vector<std::vector<complex>>(2,std::vector<complex>(2)));
    sigma[0] ={{{0.0,0.0},{1.0,0.0}},{{1.0,0.0},{0.0,0.0}}};
    sigma[1] ={{{0.0,0.0},{0.0,-1.0}},{{0.0,1.0},{0.0,0.0}}};
    sigma[2] ={{{1.0,0.0},{0.0,0.0}},{{0.0,0.0},{-1.0,0.0}}};
    

    int L = U_gauge.size();

    std::vector<std::vector<std::vector<complex>>> gauge_force = g_force(U_gauge,BETA);
    std::vector<std::vector<std::vector<complex>>> F(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2,complex(0.0,0.0))));
    std::vector<std::vector<std::vector<complex>>> xi = cg(phi,phi,U_gauge,m0,1000,1e-8);
    std::vector<std::vector<std::vector<complex>>> chi = f_M(xi,U_gauge,m0);

    
    complex temp;
    complex temp_2;
    int t ;
    double bc;

    for(int x =0; x<L;x++){
        for(int y = 0; y < L; y++){
            temp  = complex(0.0,0.0);
            temp_2  = complex(0.0,0.0);
            for(int alpha = 0; alpha < 2; alpha++){
                for(int beta = 0; beta<2 ; beta++){
                    t = x;
                    bc = temporal_bc(t, L, +1);
                    temp += std::conj(chi[x][y][alpha])*(complex(1.0,0.0)-sigma[0][alpha][beta])*U_gauge[x][y][0]*xi[(x+1)%L][y][beta]*bc  //xi[(x+1)%L][y][0] 
                    - std::conj(bc*chi[(x+1)%L][y][alpha])*(complex(1.0,0.0)+sigma[0][alpha][beta])*std::conj(U_gauge[x][y][0])*xi[x][y][beta];


                    temp_2 += std::conj(chi[x][y][alpha])*(complex(1.0,0.0)-sigma[1][alpha][beta])*U_gauge[x][y][1]*xi[x][(y+1)%L][beta]  //xi[(x+1)%L][y][1]  
                    - std::conj(chi[x][(y+1)%L][alpha])*(complex(1.0,0.0)+sigma[1][alpha][beta])*std::conj(U_gauge[x][y][1])*xi[x][y][beta];

                }
            }
            F[x][y][0] =(temp).imag() +    gauge_force[x][y][0];
            F[x][y][1] =(temp_2).imag() +  gauge_force[x][y][1];

        }

    }
    return F;

}

void leapfrog(double beta, double m0,std::vector<std::vector<std::vector<complex>>> phi,std::vector<std::vector<std::vector<complex>>>& init_gauge,std::vector<std::vector<std::vector<complex>>>& init_pi
    ,int MD_steps,int MD_trajectory_length, int MD_direction){

       double MD_step_size = (double) MD_direction * (double) MD_trajectory_length / (double) MD_steps;
       std::vector<std::vector<std::vector<complex>>> pi = init_pi;
       std::vector<std::vector<std::vector<complex>>>  U = init_gauge;
       int L = init_gauge.size();
       std::vector<std::vector<std::vector<complex>>> force;

     
        // U' = e^(i/2*dt*P)*U
        force = net_force(phi,U,beta,m0);
       for(int i = 0;i <L ;i++){
        for(int j = 0; j < L; j++){
            pi[i][j][0] += 0.5*MD_step_size*force[i][j][0].real();
            pi[i][j][1] += 0.5*MD_step_size*force[i][j][1].real();
        }
       }
    
        for(int x = 0; x<MD_steps;x++){
            
            for(int i = 0;i <L ;i++){
                for(int j = 0; j < L; j++){
                    U[i][j][0] *= std::exp(complex(0.0,1.0)*MD_step_size*pi[i][j][0].real());
                    U[i][j][1] *= std::exp(complex(0.0,1.0)*MD_step_size*pi[i][j][1].real());
                }
            }
            force = net_force(phi,U,beta,m0);
            for(int i = 0;i <L ;i++){
                for(int j = 0; j < L; j++){

                    if(x<MD_steps-1){
                        pi[i][j][0] += MD_step_size*force[i][j][0].real();
                        pi[i][j][1] += MD_step_size*force[i][j][1].real();
                    }
                    else{
                        pi[i][j][0] += 0.5*MD_step_size*force[i][j][0].real();
                        pi[i][j][1] += 0.5*MD_step_size*force[i][j][1].real();
                    }

                   
                }
            }
        }

        init_gauge = U;
        init_pi = pi;

}

double topological_charge(const std::vector<std::vector<std::vector<complex>>>& U_gauge){
    int N_t = U_gauge.size();
    int N_x = U_gauge[0].size();
    double sum = 0.0;
    std::vector<int> n;
    int mu = 0; int nu = 1;

    for(int i = 0;i <N_t ;i++){
        for(int j = 0; j < N_x; j++){
            n.push_back(i);
            n.push_back(j);
            sum += std::arg(plaquette(U_gauge,n,mu,nu));
            n.clear();
        }
    }    
    return sum / (2.0 * M_PI);
}

void HMC(std::vector<std::vector<std::vector<complex>>>& U_gauge, int configs, int MD_steps,int MD_trajectory_length,std::mt19937 & gen, double beta,double m0,File & file ){
    int L = U_gauge.size();
    int MD_direction = 1;
    std::vector<std::vector<std::vector<complex>>> init_pi(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    std::vector<std::vector<std::vector<complex>>> init_phi(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    std::vector<std::vector<std::vector<complex>>> New_gauge;
    std::vector<std::vector<std::vector<complex>>> New_pi(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));

    double H_old;
    double H_new;
    double dH;
    double random_num;
    std::uniform_real_distribution<double> dist(0.0,1.0);
    int acceptence =0;
    int rejected =0;
    //std::cout<< g_plaquette(U_gauge,beta)<< "\n";

    

    
    

    for(int i =0;i<configs;i++){
        generate_Phi(init_phi,gen);
        generate_Pi(init_pi,gen);
        init_phi = f_Mdag(init_phi,U_gauge,m0);
        //std::cout<<"Plaquett element [0,0,0] [0,0,1]"<< U_gauge[1][1][0]<<" "<< U_gauge[1][1][1]<<"\n";

        H_old = hamiltonian(init_pi,init_phi,U_gauge,beta,m0);
        
        New_gauge = U_gauge;
        New_pi = init_pi;
        leapfrog(beta,m0,init_phi,New_gauge, New_pi,MD_steps,MD_trajectory_length,MD_direction);

        H_new = hamiltonian(New_pi,init_phi,New_gauge,beta,m0);
        dH = H_new - H_old;
        //std::cout << H_new<<"<-New_H Old_H-> "<< H_old<<" dH = "<< dH << " \n";
        random_num = dist(gen);
        if (random_num < std::exp(-dH)){
            acceptence++;
            U_gauge = New_gauge;
            // std::cout<< std::setprecision(16) << g_plaquette(U_gauge,beta)<< "\n";
            // std::cout<< std::setprecision(16) <<"Topological Charge: "<<topological_charge(U_gauge)<<"\n";
            file.createDataSet("first_sim/field_"+std::to_string(i),U_gauge);
        }
        else{
            rejected++;
                //std::cout<<"Rejected :("<<rejected << "\n";
            }
        
    }
    double accept_rate = (double) acceptence / (double) configs * 100.0;
    std::cout<<"-t "<<MD_steps<<" -e "<<MD_trajectory_length<<" Accept Rate: "<<accept_rate << " %\n";

}




void check_leapfrog( int argc, char *argv[], std::vector<std::vector<std::vector<complex>>> & U_gauge,std::mt19937 & gen){
    int L = U_gauge.size();

    int MD_steps = 50; int MD_trajectory_length = 1; double beta = 1.0;
    int configs;
    std::string none;
    GetUserParam(argc,argv,none,MD_trajectory_length,MD_steps,configs,L);
    
    std::vector<std::vector<std::vector<complex>>> init_phi(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    std::vector<std::vector<std::vector<complex>>> init_pi(L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));

    generate_Phi(init_phi,gen);
    generate_Pi(init_pi,gen);
    std::vector<std::vector<std::vector<complex>>> old_gauge = U_gauge;
    std::vector<std::vector<std::vector<complex>>> old_pi = init_pi;


    init_phi = f_Mdag(init_phi,U_gauge,1.0);
    leapfrog(1.0,1.0,init_phi,U_gauge,init_pi,MD_steps,MD_trajectory_length,1);
    for(int i = 0;i <L ;i++){
        for(int j = 0; j < L; j++){
            init_pi[i][j][0] = -1.0 * init_pi[i][j][0];
            init_pi[i][j][1] = -1.0 * init_pi[i][j][1];
        }
       }
    std::cout<<"-----------------Forward Done Going Backward---------------------\n";
   leapfrog(1.0,1.0,init_phi,U_gauge,init_pi,MD_steps,MD_trajectory_length,1);
    
    for(int i = 0;i <L ;i++){
        for(int j = 0; j < L; j++){
            std::cout<< "U'' - U  = "<<std::setprecision(16)<< U_gauge[i][j][0] - old_gauge[i][j][0] <<" " << U_gauge[i][j][1] - old_gauge[i][j][1]<<"\n";
        }
       }


}

complex scalar_product_complex(std::vector<std::vector<std::vector<complex>>> p, std::vector<std::vector<std::vector<complex>>> t){
    int L = p.size();
    complex x_dir = complex(0.0,0.0);
    for(int x = 0;x<L;x++){
        for(int y = 0;y<L;y++){
            x_dir += (std::conj(p[x][y][0])*t[x][y][0]) + (std::conj(p[x][y][1])*t[x][y][1]);
        }
    }
    return x_dir;
}

void check_hermitian(const  std::vector<std::vector<std::vector<complex>>> & lattice,std::mt19937 & gen){
    int L = lattice.size();
    std::vector<std::vector<std::vector<complex>>> a_field (L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));
    std::vector<std::vector<std::vector<complex>>> b_field (L,std::vector<std::vector<complex>>(L,std::vector<complex>(2)));

    generate_Phi(a_field,gen);
    generate_Phi(b_field,gen);
    std::vector<std::vector<std::vector<complex>>> right = f_M(b_field,lattice,1.0);

    complex first = scalar_product_complex(a_field,right);
    right = f_Mdag(a_field,lattice,1.0);
    complex second = scalar_product_complex(right,b_field);

    complex diff = first - second ;
    
    std::cout<< "diff: "<< diff <<"first "<<first <<"second "<<second <<"\n";

    std::vector<std::vector<std::vector<complex>>> D_d_D = f_Mdag_f_M(a_field,lattice,1.0);
    complex third = scalar_product_complex(D_d_D,a_field);
    std::cout<< "third test for M⁺M, should be real!: "<< third <<"\n";



}

void GetUserParam( int argc, char *argv[], std::string & filename, int & epsilon, int & tau, int & config, int & L){

/* Variablen: */
    int i;
    char* endptr;
    const char usage[] = 
        "Parameters [-e <MD Step> -t <Integration Time> -c <Configuration> -l <Lattice Size>] -f <Filename>";
    const char error_message[] =
        "# FEHLER(GetuserParam): falsche Option: ";

    if (argc>1) { /* falls es ueberhaupt Parameter gibt ... */
        for (i=1; i<argc; i++){
            /* parameter 2 Charakter lang und sollte mit '-' anfaengen ... */
            if ( (std::strlen(argv[i])==2) && (argv[i][0] == '-') ) { 
                switch (argv[i][1]) { 
                    case 'e':
                        epsilon = strtod(argv[++i],&endptr);
                        
                        break;
                    case 't':
                        tau = strtod( argv[++i], &endptr);
                        
                        break;
                    case 'c':
                        config = strtod( argv[++i], &endptr);
                        
                        break;
                    case 'l':
                        L = strtod( argv[++i], &endptr);
                        
                        break;
                    case 'f':
                        filename = argv[++i];                        
                        break;
                    default:
                        break;
                }
            } else {
                std::cout << error_message << std::endl << usage << std::endl;
                exit(1);
            } /* end of: if-then-else */
        } /* end-of: for */ 
    } /* end-of: if */
}