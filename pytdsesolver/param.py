import numpy as np


w0 = 0.057
dt = 0.02
pi = np.pi
param = {"n_threads": 72,
         "init_wf": 1,
         "use_potential": 0,
         "geometry": 3,
         "ni": 1000,
         "imin": -100.0,
         "imax": 100.0,
         "nj": 1000,
         "jmin": -100.0,
         "jmax": 100.0,
         "nk": 1000,
         "kmin": -100.0,
         "kmax": 100.0,
         "w0": w0,
         "period": 2*pi/w0,
         "tmax_ev": 4*2*pi/w0,
         "tmax_sim": 5*2*pi/w0,
         "dt": dt,
         "dt_ITP": 0.004,
         "nt": int(5*2*pi/w0/dt),
         "nt_ITP": 200,
         "nt_diag": 100,
         "env": 0,
         "w0Ei": w0,
         "w0Ej": w0, 
         "w0Ek": w0,
         "w0Bi": w0,
         "w0Bj": w0,
         "w0Bk": w0,
         "E0i": 0.000,
         "E0j": 0.000,
         "E0k": 0.067,
         "B0i": 0.000,
         "B0j": 0.000,
         "B0k": 0.000,
         "phiEi": 0.0*pi,
         "phiEj": 0.0*pi,
         "phiEk": 0.0*pi,
         "phiBi": 0.0*pi,
         "phiBj": 0.0*pi,
         "phiBk": 0.0*pi,
         "n_probes": 3,
         "probe_def": ("acc_i,results/acc_i_B000.dat;"
                       "acc_j,results/acc_j_B000.dat;"
                       "acc_k,results/acc_k_B000.dat;"
                       )
        }
