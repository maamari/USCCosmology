from cosmoboost.lib import jeong
from scipy.ndimage import gaussian_filter as GF
import healpy as hp
import numpy as np
import random
import numpy.random as rand
import os
import matplotlib.pyplot as plt

plt.rcParams['font.size']=13
plt.rcParams['font.family']='stix'
plt.rcParams['text.usetex']=False
plt.rcParams['figure.figsize']= (6.5,4)
plt.rcParams['figure.dpi']=150

import cosmoboost as cb

#read the default parameters from cosmoboost
pars = cb.DEFAULT_PARS

lmax=pars['lmax']=6000
delta_ell = pars['delta_ell']=10
beta=pars['beta']
T_0 = pars["T_0"]

# initialize the kernel 
pars['d']=0
pars['method']='analytic'
kernel_a = cb.Kernel(pars,save_kernel=False, overwrite=True)

def func(seed):
    random.seed(seed)
    rand.seed(seed)    
    T_0 = 2.725
    ell=np.arange(lmax+1)

    # here's a sample power spectrum generated with CAMB
    lib_dir = os.path.join(cb.COSMOBOOST_DIR,"lib")

    Cl_camb = np.load(os.path.join(lib_dir,"Cl_Planck2018_camb.npz"))

    Cl_TT = Cl_camb["TT"][:lmax+1]
    Cl_EE = Cl_camb["EE"][:lmax+1]
    Cl_BB = Cl_camb["BB"][:lmax+1]
    Cl_TE = Cl_camb["TE"][:lmax+1]

    # let's use it to simulate a CMB map
    Cl = np.array([Cl_TT,Cl_EE,Cl_BB,Cl_TE])
    alm_T, alm_E, alm_B = hp.synalm(Cl,lmax=lmax,new=True,verbose=True)

    # this is our alm in the rest frame
    alm_r = np.array([alm_T, alm_E, alm_B])

    # this is the power spectrum of the simulation
    Cl_r = hp.alm2cl(alm_r)
    
    alm_T_r = alm_T
    alm_T_b_a = cb.boost_alm(alm_T_r,kernel_a)
    
    # calculate the temperature power spectrum in the rest and boosted frame
    Cl_TT_r = Cl_r[0] 
    Cl_TT_b_a = hp.alm2cl(alm_T_b_a)
    Cl_TT_jeong_fs = jeong.jeong_boost_Cl_2ndorder(ell, Cl_TT_b_a,beta= 0.00123)
    
    # calculate the relative change of the boosted Cl using the accurate formalism
    dCl_TT_b_a = (Cl_TT_b_a - Cl_TT_r)
    dCl_Cl_TT_b_a = dCl_TT_b_a/Cl_TT_r

    dCl_TT_b_jeong_fs = (Cl_TT_jeong_fs - Cl_TT_r)
    dCl_Cl_TT_b_jeong_fs = Cl_TT_jeong_fs/Cl_TT_r

    # Gaussian smooth the results with a window
    dL = 50

    dCl_Cl_TT_b_a_GF = GF(dCl_Cl_TT_b_a[:lmax-20], dL, mode="constant")
    dCl_Cl_TT_b_jeong_GF_fs = GF(dCl_Cl_TT_b_jeong_fs, dL, mode="constant")
        
    def mask_cutbelowlat(cut_angle, lat_pix):
        mask = np.ones_like(lat_pix)
        mask[(lat_pix > cut_angle)] = 0.

        #approximate f_sky using the number of pixels
        f_sky = len(mask[mask == 1.]) / len(mask)

        return mask, f_sky
    
    NSIDE=2048
    NPIX = hp.nside2npix(NSIDE)
    lon_pix, lat_pix = hp.pix2ang(NSIDE,np.arange(NPIX),lonlat=True)
    
    mask_60,f_sky = mask_cutbelowlat(0,lat_pix)
    
    hp.mollview(mask_60, sub=221,title="Mask")
    plt.savefig("mask.png")
    plt.clf()
    
    Cl_TT_jeong = jeong.jeong_boost_Cl_1storder(ell, Cl_TT_b_a,beta= 0.00123, cos_avg=-0.5)
    dCl_Cl_TT_b_jeong = Cl_TT_jeong/Cl_TT_r
    
    boost_ang = -(lat_pix-90.)
    
    cos_map = np.cos(np.deg2rad(boost_ang))*mask_60
    ones_map = np.ones_like(lat_pix)*mask_60
    cos_avg = np.sum(cos_map)/np.sum(ones_map)
    cos2_avg = np.sum(cos_map*cos_map)/np.sum(ones_map)    
    #jeong.jeong_boost_Cl()
    
    T_map_r = hp.alm2map(alm_T_r,NSIDE)
    T_map_b_a = hp.alm2map(alm_T_b_a,NSIDE)
    T_map_r_ma = mask_60*T_map_r
    T_map_b_a_ma = mask_60*T_map_b_a
    
    T_map_b_a = hp.alm2map(alm_T_b_a,NSIDE)
    T_map_b_a_ma = mask_60*T_map_b_a
    
    Cl_TT_r_ma =(1/f_sky)*hp.anafast(T_map_r_ma,lmax=lmax)
    Cl_TT_b_a_ma =(1/f_sky)*hp.anafast(T_map_b_a_ma,lmax=lmax)
    
    # calculate the relative change of the boosted Cl using the accurate formalism
    dCl_TT_b_a_ma = (Cl_TT_b_a_ma - Cl_TT_r_ma)
    dCl_Cl_TT_b_a_ma = dCl_TT_b_a_ma/Cl_TT_r_ma

    # Gaussian smooth the results with a window
    dL = 50

    dCl_Cl_TT_b_a_ma_GF = GF(dCl_Cl_TT_b_a_ma[:lmax-20], dL, mode="constant")
    
    Cl_TT_jeong = jeong.jeong_boost_Cl_1storder(ell, Cl_TT,beta= 0.00123, cos_avg=cos_avg)
    dCl_TT_b_jeong = (Cl_TT_jeong - Cl_TT_r)
    dCl_Cl_TT_b_jeong = Cl_TT_jeong/Cl_TT_r
    dCl_Cl_TT_b_jeong_GF = GF(dCl_Cl_TT_b_jeong[:lmax-20], dL, mode="constant")
    
    return ell, 100*dCl_Cl_TT_b_a_GF, 100*dCl_Cl_TT_b_jeong_GF_fs, 100*(dCl_Cl_TT_b_a_ma_GF), 100*(dCl_Cl_TT_b_jeong_GF)

numRuns = 2
ell = []
fullANL = []
fullJEONG = []
maskedANL = []
maskedJEONG = []

for i in range(numRuns):
    run = func(i)
    ell.append(run[0])
    fullANL.append(run[1])
    fullJEONG.append(run[2])
    maskedANL.append(run[3])
    maskedJEONG.append(run[4])
    
fullANL_avg = np.array(fullANL[0])
fullJEONG_avg = np.array(fullJEONG[0])
for i in range(numRuns-1):
    fullANL_avg += np.array(fullANL[i+1])
    fullJEONG_avg += np.array(fullJEONG[i+1])
fullANL_avg /= numRuns
fullJEONG_avg /= numRuns

maskedANL_avg = np.array(maskedANL[0])
maskedJEONG_avg = np.array(maskedJEONG[0])
for i in range(numRuns-1):
    maskedANL_avg += np.array(maskedANL[i+1])
    maskedJEONG_avg += np.array(maskedJEONG[i+1])
maskedANL_avg /= numRuns
maskedJEONG_avg /= numRuns

plt.plot(maskedANL_avg,linewidth=2,color="tab:blue",label='analytic $a_{\ell m}$ (Gaussian smoothed)')
plt.plot(maskedJEONG_avg,linewidth=2,color="tab:red",label='jeong')

err=[]
for i in range(lmax-20):
	err.append(np.var(np.concatenate(maskedANL).reshape(numRuns,lmax-20)[:,i]))
plt.errorbar(ell[0][:lmax-20],maskedANL_avg,err,lw=0.1)

errJ=[]
for i in range(lmax-20):
	errJ.append(np.var(np.concatenate(maskedJEONG).reshape(numRuns,lmax-20)[:,i]))
plt.errorbar(ell[0][:lmax-20],maskedJEONG_avg,errJ,lw=0.1)

plt.xlabel("$\ell$")
plt.ylabel("$\Delta C_\ell/C_\ell (\%)$")

plt.grid()
plt.legend()

plt.xlim(100,lmax)
plt.savefig("masked.png")
