import numpy as np
import matplotlib.pyplot as plt

outfile = open("results.txt","w")

def uncertainty(vector):
    temp = np.sqrt(sum((vector-np.mean(vector))**2)/(len(vector)-1))
    return temp

one_big_exp = np.loadtxt("../data/oneBigexperiment.txt")
radii = one_big_exp[:,0]
radii_square = one_big_exp[:,1]

fig,ax = plt.subplots(1,2,figsize=(10,10))

ax[0].hist(radii)
ax[0].set_xlabel(r"Radius")
ax[0].set_ylabel(r"#",loc='top')
ax[0].set_title("Histogram of Radii")
ax[1].hist(radii_square)
ax[1].set_ylabel(r"#",loc='top')
ax[1].set_xlabel(r"Radius Squared")
ax[1].set_title("Histogram of Radius Squared")

plt.savefig("One_Big_Experiment_Radii_histogram")
plt.close()

indicator = 4*(radii_square <=1)
mean = np.mean(indicator)
std = uncertainty(indicator)
outfile.write(f'Mean of Indicator: {mean:.5} Std: {std:.5}\n')
n, bins, patches = plt.hist(indicator,label='Histgrom of indicator')
plt.vlines(mean,0,n.max(),colors='r',label='Mean of samples')
plt.vlines(np.pi,0,n.max(),colors='b',label=r'$\pi$')
plt.fill_between([mean-std,mean+std],0,8000,alpha=0.2,label='Standard Deviation')
plt.xlabel(r'Indicatior variable')
plt.ylabel(r'#',loc='top')
plt.legend()
plt.savefig('Pi_comparison')
plt.close()

pi_estimator = np.loadtxt("../data/split_100_exp_estimator.txt")
mean_pi = np.mean(pi_estimator)
std_pi = uncertainty(pi_estimator)

outfile.write(f'Mean of split exp: {mean_pi:.5} Std: {std_pi:.5}\n')


n,bins,_ = plt.hist(pi_estimator)
plt.vlines(mean_pi,0,n.max(),colors='r',label='Mean of estimators')
plt.vlines(np.pi,0,n.max(),colors='b',label=r'$\pi$')
plt.fill_between([mean_pi-std_pi,mean_pi+std_pi],0,n.max(),alpha=0.2,label='Standard Deviation')

plt.xlabel(r'$\pi_x$')
plt.ylabel(r'Experiment #',loc='top')
plt.legend()
plt.savefig('split_experiment')
plt.close()

zillionexp = np.loadtxt("../data/zillion_little_exp_estimator.txt")
mean_zillion = np.mean(zillionexp)
std_zillion = uncertainty(zillionexp)
outfile.write(f'Mean of zillion exp: {mean_zillion:.5} Std: {std_zillion:.5}\n')


def open_file(file_path):

    data = np.loadtxt(file_path)
    return data


def Pi_estimator(i,j):
    file_path = f"../data/combination_pairs_{i}_exp_{j}.txt"
    data = open_file(file_path)
    mean = np.mean(data)
    error = uncertainty(data)
    return mean,error

def all_combination():
    combination = np.array([10,100,1000,10000])
    mean_histogram = np.array([])
    error_histogram = np.array([])
    plt.figure(figsize=(10,10))
    for i in combination:
        for j in combination:
            mean,error = Pi_estimator(i,j)
            mean_histogram = np.append(mean_histogram,mean)
            error_histogram = np.append(error_histogram,error)
            plt.errorbar(i*j,mean-np.pi,yerr=error,fmt='+',label=f'P={i} and X={j}')
            outfile.write(f'{i} & {j} & {i*j} & {mean:.5} & {error:.5} \ \ \n')
            
    plt.hlines(0,0,(i*j),label=r'$\pi$')
    plt.xscale('log')
    plt.xlabel(r'XP')
    plt.ylabel(r'$\pi_x$-$\pi$ with uncertainties')
    plt.legend(loc='center left', bbox_to_anchor=(1,0.5))
    plt.savefig('best_XP',bbox_inches='tight')
    plt.close()
    return mean_histogram,error_histogram


mean_estimator,error_estimator = all_combination()
std_pi = uncertainty(mean_estimator)
mean_pi = np.mean(mean_estimator)
outfile.write(f'Mean of longer exp: {mean_pi:.5} Std: {std_pi:.5}\n')

n,counts,_ = plt.hist(mean_estimator,label=r'Estimator of $\pi_x$')
plt.vlines(np.pi,0,n.max(),colors='r',label=r'$\pi$')
plt.vlines(mean_pi,0,n.max(),colors='b',label=r'Final estimate')
plt.fill_between([mean_pi-std_pi,mean_pi+std_pi],0,n.max(),alpha=0.2,label='Uncertainties of \nfinal estimate')
plt.xlabel(r'$\pi_x$')
plt.ylabel(r'Experiment #',loc='top')
plt.legend()
plt.savefig('final_estimate')
plt.close()

def uncertainty_as_function_of(variable):
    variable = variable.upper()
    combination = np.array([10,100,1000,10000])

    for i in combination:
        for j in combination:
            mean,error = Pi_estimator(i,j)
            if variable=="X":
                plt.scatter(j,error,marker='+',label=r'$\Delta$$\pi$ %i'%(i))
            elif variable=="P":
                plt.scatter(i,error,marker='+',label=r'$\Delta$$\pi$ %i'%(j))
            else:
                print("Variable not recognized")
                break 
    if variable=="X":
        plt.xscale('log')
        plt.xlabel(r'X')
        plt.legend(loc='center left', bbox_to_anchor=(1,0.5))
        plt.savefig('error_as_X',bbox_inches='tight')
        plt.close()
    
    if variable=="P":
        plt.xscale('log')
        plt.xlabel(r'P')
        plt.legend(loc='center left', bbox_to_anchor=(1,0.5))
        plt.savefig('error_as_P',bbox_inches='tight')
        plt.close()
    return 

uncertainty_as_function_of("X")
uncertainty_as_function_of("P")  
outfile.close()
