import numpy as np
def linear_function(x,a,b):
    return x *  a + b


def plaquette_sum(U):
    """
    Compute total plaquette sum and average plaquette
    for a 2D U(1) lattice.

    Parameters:
        U : complex ndarray of shape (Nt, Nx, 2)

    Returns:
        plaq_avg : float
    """
    Nt, Nx, Nd = U.shape

    plaq_sum = 0.0

    for t in range(Nt):
        for x in range(Nx):
            # periodic boundary conditions
            tp = (t + 1) % Nt
            xp = (x + 1) % Nx

            U0 = U[t, x, 0]
            U1_t = U[tp, x, 1]
            U0_x = np.conj(U[t, xp, 0])
            U1 = np.conj(U[t, x, 1])

            plaquette = U0 * U1_t * U0_x * U1
            plaq_sum += np.real(plaquette)

    plaq_avg = plaq_sum / (Nt * Nx)

    return plaq_avg

def topological_charge(lattice):
    mu = 0
    nu = 1
    #p = lattice[x][y][mu]*lattice[(x+1)%L][y][nu]*std::conj(lattice[x][(y+1)%L][mu])*std::conj(lattice[x][y][nu]);
    p = lattice[:,:,mu]*np.roll(lattice[:,:,nu],1,0)*np.conj(np.roll(lattice[:,:,mu],1,1))*np.conj(lattice[:,:,nu])
    arg = np.angle(p)
    arg_sum = np.sum(arg) / (2*np.pi)
    return arg_sum

def wilson_loop(U, R, T, mu=0, nu=1):
    Nt, Nx, Nd = U.shape

    vals = np.zeros((Nt, Nx), dtype=np.float64)

    for t0 in range(Nt):
        for x0 in range(Nx):
            W = 1.0 + 0.0j
            t, x = t0, x0

            # +nu (space)
            for _ in range(R):
                W *= U[t, x, nu]
                x = (x + 1) % Nx

            # +mu (time)
            for _ in range(T):
                W *= U[t, x, mu]
                t = (t + 1) % Nt

            # -nu : use conj(U[t, x-1, nu]) then step back
            for _ in range(R):
                xm1 = (x - 1 + Nx) % Nx
                W *= np.conjugate(U[t, xm1, nu])
                x = xm1

            # -mu : use conj(U[t-1, x, mu]) then step back
            for _ in range(T):
                tm1 = (t - 1 + Nt) % Nt
                W *= np.conjugate(U[tm1, x, mu])
                t = tm1

            vals[t0, x0] = np.real(W)

    return np.mean(vals)


def running_mean(data):
    """
    Takes gauge field entire data.\n
    Returns running mean of plaquett sum
    """
    sim,L,L,dim = data.shape
    mean_array = np.zeros(sim)
    plaquett_sum = np.zeros(sim)
    for i in range(sim):
        plaquett_sum[i] = plaquette_sum(data[i,:,:,:])
        if i == 0:
            mean_array[i] = plaquett_sum[i]
            continue
        mean_array[i] = np.sum(plaquett_sum[:i]) /(i+1)
    return mean_array


def data_blocking(data,l):
    N = data.size
    blocked_data = np.array([])
    for i in range(0,N,l):
        a = np.mean(data[i:i+l])
        blocked_data = np.append(blocked_data,a)
    
    return blocked_data

def data_blocking_error(data,length):
    result = np.zeros((3,length))
    for i in range(1,length):
        N =len(data)
        x = []
        for j in range(0,N,i):
            a = data[j:j+i]
            x.append(np.std(a))
        result[0,i] = i
        result[1,i] = np.mean(x)
        result[2,i] = np.std(x)/np.sqrt(len(x))
   
    return result

def bootstrap_data(data):
    N = data.size
    new_data = np.zeros(N)
    for i in range(N):
        new_data[i] = data[np.random.randint(0,N-1)] 
    return new_data

def autocorr(data,lag):
    mean = np.mean(data)
    C_0 = np.mean((data-mean)**2)
    if(C_0 == 0):
        return 0
    if(lag>0):
        C_X = np.mean((data[lag:]-mean)*(data[:-lag]-mean))
    else:
        C_X = np.mean((data-mean)**2)
    return C_X/C_0

def autocorr_array(data,lag):
    x= np.zeros(lag)
    for i in range(lag):
        x[i] = autocorr(data,i)
    return x

def integrated_autocorrelation_time(data, window):
    """
    Compute the integrated autocorrelation time for a given cutoff window.
    """

    data = np.asarray(data)
    N = len(data)

    # subtract mean
    mean = np.mean(data)
    var = np.var(data)

    if var == 0:
        raise ValueError("Variance is zero, autocorrelation undefined.")

    # compute autocorrelation function up to window
    rho = np.zeros(window + 1)

    for t in range(window + 1):
        cov = np.sum((data[:N - t] - mean) * (data[t:] - mean)) / (N - t)
        rho[t] = cov / var

    # integrated autocorrelation time
    tau_int = 0.5 + np.sum(rho[1:])

    return tau_int


def W(data,L,T, tau = 0.5):
    R = np.linspace(1,L,L,dtype=np.int64)
    N = data.shape[0]
    W = np.zeros((N,L))
    for l in range(N):
        for i in R:
           W[l,i-1] = wilson_loop(data[l],i,T)

    pot = np.mean(W,0)
    tau = integrated_autocorrelation_time(W[:,2],10)

    pot_err = np.sqrt( 2 * tau * np.var(W,0) / N)

    return R, pot, pot_err




