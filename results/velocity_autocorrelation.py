import numpy as np

def q_transverse(k_sq, t):
    q_t = (1./rho)*np.exp(-nu*k_sq*t)
    return q_t

def q_longitudinal(k_sq, t):
    omega = k_sq*nu_wave*np.sqrt(4.*c*c/(k_sq*nu_wave*nu_wave)-0.5)
    q_l = (1./rho)*np.exp(-k_sq*nu_wave*t/2)*(np.cos(omega*t) - np.sqrt(k_sq*nu_wave*nu_wave/(4.*c*c-k_sq*nu_wave*nu_wave))*np.sin(omega*t))
    return q_l




rho = 10
nu = 0.8705
n_v = 0.3/rho
nu_wave = (4./3.)*nu + n_v
T = 1.0
L = 30
V = np.power(L,3)
c = np.sqrt(T)
k_const = 2.*np.pi/L
k_const_sq = np.power(k_const,2)
n_max = 50
N_c = 10
h = 0.1
gamma = 2.*(1.-np.cos(130./180.))*(np.exp(-N_c)+N_c +1)/(3.*N_c)
D = h*T*((1./gamma)-0.5)

c_v = np.zeros(50)
time = np.linspace(0,50)
f = open("cv_T"+str(T)+"_L"+str(L)+".dat", "w")

for t in time:
    for i in range(n_max):
        for j in range(n_max):
            for k in range(n_max):
                if i==0 and j == 0 and k == 0:
                    continue
                k_sq = k_const_sq*(i*i+j*j+k*k)
                c_v[t] += 2.*(2.*q_transverse(k_sq,t)+q_longitudinal(k_sq,t))*np.exp(-k_sq*D*t)
    c_v[t] *= T/V

for t in time:
    f.write("%f %f \n" %(t, c_v[t]))



