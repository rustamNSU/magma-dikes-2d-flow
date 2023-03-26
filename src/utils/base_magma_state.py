import numpy as np

def beta_T(T, TL, TS):
    Tprime = (T - TS) / (TL - TS)
    beta = 1.0 - Tprime
    beta = np.minimum(1.0, beta)
    beta = np.maximum(0.0, beta)
    return beta


def mu_melt(T):
    mu_max = 1e10
    A = -4.55
    B = 5.03e3
    C = 604.0
    Ttmp = np.maximum(C + 50.0, T)
    mu = np.minimum(mu_max, 10.0**(A + B / (Ttmp - C)))
    return mu


def theta(phi):
    phic = 0.6
    phitmp = np.minimum(phic - 1e-4, phi)
    theta = np.power((1.0 - phitmp / phic), -2.5)
    return theta


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    TL = 1250 + 273.15
    TS = 1000 + 273.15
    T = np.linspace(TS - 100, TL + 100, 100)
    fig, ax = plt.subplots(4, 1, sharex=True)
    
    beta = beta_T(T, TL=TL, TS=TS)
    ax[0].plot(T - 273.15, beta, linewidth=2)
    ax[1].semilogy(T - 273.15, theta(beta), linewidth=2)
    ax[2].semilogy(T - 273.15, mu_melt(T), linewidth=2)
    ax[3].semilogy(T - 273.15, theta(beta) * mu_melt(T), linewidth=2)
    ax[0].invert_xaxis()
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    ax[3].grid()
    
    print("beta(1300) = {}".format(beta_T(1300, TL, TS)))
    print("mu_melt(1300) = {}".format(mu_melt(1300)))
    print("theta(1300) = {}".format(theta(beta_T(1300, TL, TS))))    
    plt.show()