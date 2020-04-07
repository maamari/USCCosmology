from __future__ import division

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.integrate import simps, trapz, quad, dblquad, nquad
from scipy import constants
from scipy.special import erfc

import camb
from camb import model, initialpower

import helper_functions as hp

# =======================================================================================
# constants
delta_crit = 1.686
h = 0.697
rho_0 =  160817441770  #42285353995.6    # rho_bar_0; M_sun/Mpc^3
rho_c0 = 1.37e11                         # M_sun/Mpc^3
delta_vir = 200 # 353
D_H = 2997  # Mpc h^-1
#D_H = 2997/0.7  # Mpc
# =======================================================================================

class Cosmos():

    def __init__(self, z, kmin, kmax, N, Y200rhoc = 0.1,
                 H0 = 67, Omegab0 = 0.02256/0.67**2, Omegam0 = (0.1142+0.02256)/0.67**2,
                 ns = 0.971, Tcmb0 = 2.75, fnl = 5,
                 Mlim = None, Ps = None):

        self.z = z
        self.kmin = kmin
        self.kmax = kmax
        self.N = N
        self.Mlim = Mlim
        self.Y200rhoc = Y200rhoc
        self.Ps = Ps
        self.fnl = fnl

        self.Tcmb0 = Tcmb0
        self.H0 = H0
        self.Omegab0 = Omegab0
        self.Omegam0 = Omegam0
        self.ns = ns

    def Md(self):

        DeltaM = 10**0.3

        Md = []
        Md_ = 1e13/DeltaM  # 1e12 prev

        for i in range(100):

            Md_ *= DeltaM
            Md.append(Md_)

        return np.array(Md)


    def Ks(self):

        return np.linspace(self.kmin, self.kmax, self.N)

    #@np.vectorize
    def T_we(self, k, bSingle_k = False):

        omb = self.Omegab0
        om0 = self.Omegam0
        h = self.H0/100.
        theta2p7 = self.Tcmb0 / 2.7

        s = 44.5 * np.log(9.83/(om0*pow(h,2))) / np.sqrt(1.0 + 10.0 * pow(omb * h * h, 0.75))
        alphaGamma = 1.0 - 0.328 * np.log(431.0 * om0 * h * h) * omb / om0 + 0.38 * np.log(22.3 * om0 * h * h) * (omb / om0) * (omb / om0)

        if bSingle_k == True:

            Gamma = om0 * h * (alphaGamma + (1.0 - alphaGamma) / (1.0 + pow(0.43 * k* h * s, 4.0)))
            q = k * theta2p7 * theta2p7 / Gamma
            C0 = 14.2 + 731.0 / (1.0 + 62.5 * q)
            L0 = np.log(2.0 * np.exp(1.0) + 1.8 * q)
            Tk = L0 / (L0 + C0 * q * q)

            return Tk

        elif bSingle_k == False:

            k = self.Ks()

            Gamma = om0 * h * (alphaGamma + (1.0 - alphaGamma) / (1.0 + pow(0.43 * k* h * s, 4.0)))
            q = k * theta2p7 * theta2p7 / Gamma
            C0 = 14.2 + 731.0 / (1.0 + 62.5 * q)
            L0 = np.log(2.0 * np.exp(1.0) + 1.8 * q)
            Tk = L0 / (L0 + C0 * q * q)

            return Tk

    def setPower(self):

        T = self.T_we(k = None, bSingle_k = False)
    	# Equation obtained from: http://www.astro.caltech.edu/~george/ay21/eaa/eaa-powspec.pdf
        Pk = self.Ks()**self.ns * T**2
    	# Normalization obtained from: https://ned.ipac.caltech.edu/level5/March02/White/White3.html
        norm = 2927876.81638

        self.Ps = Pk * norm * self.d()**2

    def P_th(self, k):

        T = self.T_we(k, bSingle_k = True)
    	# Equation obtained from: http://www.astro.caltech.edu/~george/ay21/eaa/eaa-powspec.pdf
        Pk = k**self.ns * T**2
    	# Normalization obtained from: https://ned.ipac.caltech.edu/level5/March02/White/White3.html
        norm = 2927876.81638

        return Pk * norm * self.d()**2


    def E(self, Om = 0.2815):

        """
        David Hogg, Fist Light...
        """
        return np.sqrt(Om * (1 + self.z)**3 + 1.0 - Om)


    def Xi(self, Om0 = 0.2807):

        """
        Comoving radial distance
        """

        return (D_H / h) * quad(lambda z: 1/self.E(Om0), 0, self.z)[0]


    def DA(self):

        """
        Angular diameter distance
        """

        return (1 + self.z)**(-1) * self.Xi(self.z)


    def setMlim(self):

        """
        Equaton 1 from paper
        """
        Y = self.Y200rhoc            # arcmin^2
        Y = np.deg2rad(Y/(60.)**2)    # arcmin to degree and degree to radian
        self.Mlim = 1e14  #(self.DA()**2 * self.E()**(-2/3) * (Y/2.5e-4))**0.533 * 1e15   # M_solar


    def D_plus(self, Om = 0.2815):

        """
        HMF paper
        """
        integrand = lambda z_: (1 + z_)/(np.sqrt(self.Omegam0 * (1 + z_)**3 + 1.0 - self.Omegam0))**3
        I = quad(integrand, self.z, np.inf, limit = 1500)[0]

        return 5./2 * Om * self.E() * I


    def D_plus_0(self, Om = 0.2815):

        """
        HMF paper
        """
        integrand = lambda z_: (1 + z_)/(np.sqrt(self.Omegam0 * (1 + z_)**3 + 1.0 - self.Omegam0))**3
        I = quad(integrand, 0.0, np.inf, limit = 1500)[0]

        return 5./2 * Om * np.sqrt(self.Omegam0 * (1 + 0)**3 + 1.0 - self.Omegam0) * I


    def d(self):

        """
        Eq.11, HMF paper
        """

        return self.D_plus()/self.D_plus_0()


    def R(self, M):

        """
        Radius of virialized mass
        """

        return ((3.0 * M)/(4.0 * np.pi * rho_c0))**0.33


    def W(self, M, k = None, bSingle_k = False):

        """
        Top-hat Window function, M in units of M_sun
        """

        if bSingle_k == False:

            k = self.Ks()

            return 3.0 * (np.sin(k * self.R(M)) - (k*self.R(M)) * np.cos(k*self.R(M)))/(k*self.R(M))**3.0

        elif bSingle_k == True:

            return 3.0 * (np.sin(k * self.R(M)) - (k*self.R(M)) * np.cos(k*self.R(M)))/(k*self.R(M))**3.0


    def sigma_squared(self, M, k = None, j = 0, bSingle_k = False):

        """
        Mass variance. arXiv:0712.0034v4 eq.1
        """

        if bSingle_k == False:

            k = self.Ks()

            win = self.W(M, k)

            I = (k)**(2. + 2.*j) * self.P_th(k) * win**2

            return 1.0/(2.0*np.pi**2) * trapz(I, x = k)

        elif bSingle_k == True:

            win = self.W(M, k, bSingle_k = bSingle_k)

            I = (k)**(2. + 2.*j) * self.P_th(k) * win**2

            return 1.0/(2.0*np.pi**2) * trapz(I, x = k)


    def nu(self, M):

        """
        Peak height
        """

        return delta_crit/np.sqrt(self.sigma_squared(M))


    def M_R(self, k, M):

        """

        """

        return 2.0/3 * (self.T_we(k, bSingle_k = True) * k**2)/(((self.H0/3e5)**2) * self.Omegam0) * self.W(M, k, bSingle_k = True)


    def P_phi(self, k):

        """
        Newtonian potential
        """

        A = 2927876.81638
        #H0, Omega_m0, ns = self.H0, self.Omegam0, self.ns

        return 9./4 * A * (self.H0/3e5)**4 * self.Omegam0**2 * k**(self.ns - 4)       # Divided by c = 3e5 to make the units correct


    def F_R_k(self, k, M):

        """
        Based on equation Verde Mattarasee paper
        For a single k value
        """
        fnl = self.fnl
        M_lim = self.Mlim

        sR2 = self.sigma_squared(M_lim)
        #print(sR2)
        ksi = self.Ks()# - 1e-9
        mu = np.linspace(-1.0, 1.0, len(ksi))

        kksi, mmu = np.meshgrid(ksi, mu, indexing = 'ij')

        alpha = np.sqrt(kksi**2 + k**2 + 2 * mmu * k * kksi)

        I = kksi**2 * self.M_R(kksi, M) * self.P_phi(kksi) * self.M_R(alpha, M) * (2.0 + self.P_phi(alpha)/self.P_phi(k))

        return (2 * fnl/(8. * np.pi**2 * sR2)) * simps(simps(I, ksi), mu)

    def Delta_b_nL(self, M, k):

        """
        Non-gaussian scale dependent correction to bias
        """
        Dbnl = self.F_R_k(k, M)/self.M_R(k, M)   #  for k in np.linspace(0.0001, 0.2, 55)])#self.F_R_k(k, M)/(self.M_R(k, M))

        return Dbnl


    def b_L(self, M, k = None, bScale = False):

        """
        Linear bias
        """
        a, p = 0.75, 0.3

        nu_ = delta_crit/np.sqrt(self.sigma_squared(M, bSingle_k = False))

        bL = 1 + (a * nu_**2 - 1)/delta_crit + 2 * p/(delta_crit * (1 + (a * nu_**2)**p))


        if bScale == False:

            return bL

        elif bScale == True:

            return bL + (bL - 1) * delta_crit * self.Delta_b_nL(M, k)


    def dW2dM(self, M):

        """
        From HMF paper
        """
        k = self.Ks()
        return (np.sin(k * self.R(M)) - k * self.R(M) * np.cos(k*self.R(M))) * (np.sin(k*self.R(M)) * (1 - 3./(k*self.R(M))**2)
                + 3 * np.cos(k*self.R(M))/(k*self.R(M)))


    def dlnsdlnM(self, M):

        """
        dlnSigma/dlnM from HMF paper
        """
        K_ = self.Ks()
        I = self.Ps/(K_**2) * self.dW2dM(M)

        return (3./(2 * np.pi**2 * self.R(M)**4 * self.sigma_squared(M))) * simps(I, K_)


    def f_J(self, M):

        """
        Jenkins fitting function (There is a typo in HMF?)
        """
        s = np.sqrt(self.sigma_squared(M))

        return 0.315 * np.exp(-np.absolute(0.61 + np.log(s**(-1.)*self.d()))**(3.8))


    def f_T(self, M):

        """
        Tinker fitting function
        """

        return None


    def dndlnM(self, M):

        """
        n(M,z)   ADD REDSHIFT DEPENDENCE HERE ===================================================
        """

        return (rho_c0/M) * self.f_J(M) * np.absolute(self.dlnsdlnM(M))


    def dndM(self, M):

        """
        n(M,z)
        """
        return (rho_c0/M**2) * self.f_J(M) * np.absolute(self.dlnsdlnM(M))


    def b_eff(self, m, k = None, bScale = False):     #b_eff(self, M, m):

        """
        """
        if bScale == False:

            bL = np.array([self.b_L(M_) for M_ in self.Md()])
            nM = np.array([self.dndlnM(M_) for M_ in self.Md()])

            ERF = np.array([erfc(hp.xm(self.z, M_, self.Mlim, m)) - erfc(hp.xm(self.z, M_, self.Mlim, m + 1)) for M_ in self.Md()])
            bias = trapz(nM * bL * ERF, x = self.Md())/ trapz(nM * ERF, x = self.Md())

            return bias


        elif bScale == True:

            bL = np.array([self.b_L(M_, k, bScale = bScale) for M_ in self.Md()])
            nM = np.array([self.dndlnM(M_) for M_ in self.Md()])

            ERF = np.array([erfc(hp.xm(self.z, M_, self.Mlim, m)) - erfc(hp.xm(self.z, M_, self.Mlim, m + 1)) for M_ in self.Md()])
            bias = trapz(nM * bL * ERF, x = self.Md())/ trapz(nM * ERF, x = self.Md())

            return bias


    def Ph_mn(self, m, n, k = None, bScale = False):

        """
        Galaxy power spectrum
        """
        if bScale == False:

            bm = self.b_eff(m)
            bn = self.b_eff(n)

            return bm * bn * self.Ps

        elif bScale == True:

            bm = np.array([self.b_eff(m, k, bScale = bScale) for k in self.Ks() + 0.00001]) #np.linspace(1e-3, 0.1, 30)])
            bn = np.array([self.b_eff(n, k, bScale = bScale) for k in self.Ks() + 0.00001])

            return bm * bn * self.Ps


    def Veff_mn(self, m, n, k = None, bScale = False):

        """

        """

        if bScale == False:

            Pmn = self.b_eff(m) * self.b_eff(n) * self.Ps
            Pnm = self.b_eff(n) * self.b_eff(m) * self.Ps
            Pmm = self.b_eff(m) * self.b_eff(m) * self.Ps
            Pnn = self.b_eff(n) * self.b_eff(n) * self.Ps

            Mpl = self.Md()[self.Md() < 1e16]
            nm = self.dndlnM(Mpl[m])
            nn = self.dndlnM(Mpl[n])

            n_ = (Pmn)**2 * nm * nn

            delta_nm = 1 if n == m else 0
            d1 = (nm * Pmm + 1) * (nn * Pnn + 1)
            d2 = nm * nn * (Pnm + delta_nm/nm)**2

            return n_/(d1 + d2)


        elif bScale == True:

            bm = np.array([self.b_eff(m, k, bScale = bScale) for k in self.Ks() + 0.00001])
            bn = np.array([self.b_eff(n, k, bScale = bScale) for k in self.Ks() + 0.00001])

            Pmn = bm * bn * self.Ps
            Pnm = bn * bm * self.Ps
            Pmm = bm * bm * self.Ps
            Pnn = bn * bn * self.Ps

            Mpl = self.Md()[self.Md() < 1e16]
            nm = self.dndlnM(Mpl[m])
            nn = self.dndlnM(Mpl[n])

            n_ = (Pmn)**2 * nm * nn

            delta_nm = 1 if n == m else 0
            d1 = (nm * Pmm + 1) * (nn * Pnn + 1)
            d2 = nm * nn * (Pnm + delta_nm/nm)**2

            return n_/(d1 + d2)


    def VarCrossPower(self, m, n):

        """
        Variance of cross powerspectrum term, eq A1
        """

        Nmod = 1

        Pmn = self.b_eff(m) * self.b_eff(n) * self.Ps
        Pmm = self.b_eff(m) * self.b_eff(m) * self.Ps
        Pnn = self.b_eff(n) * self.b_eff(n) * self.Ps

        Mpl = self.Md()[self.Md() < 1e16]
        nm = self.dndlnM(Mpl[m])
        nn = self.dndlnM(Mpl[n])

        delta_nm = 1 if n == m else 0

        return 1/Nmod * ((Pmm + 1/nm) * (Pnn + 1/nm) + (Pmn + delta_nm/nm)**2)