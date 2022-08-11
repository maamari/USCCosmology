import numpy as np

# Using PVP 04

# Where do we find these terms?
nbar = # ?
bias = # ?
nx = # ?
nr = # I believe would come from catalog/simulation
alpha = # Comes from catalog?

# Some initial guess on the power spectrum
Pi = # Maybe call camb? Not entirely sure

def w(bias,Pi,nbar):
    """
    Weight
    Equation 28 in PVP
    """
    return ((bias**2)*Pi) / (1.0+nbar*Pi*(bias**2))

def N(nbar,w):
    """
    Normalization
    Equation 7 in PVP
    """
    return np.sqrt(np.sum((nbar**2)*(w**2)))

def pShot(nbar,w,bias,alpha,N):
    """
    Equation 16 in PVP
    """
    summationTerm = nbar*((w**2)/(bias**2))
    return np.sum(summationTerm)*((1+alpha)/(N**2))

def Fr(w,N,bias,ng,alpha,nr):
    """
    Fluctutation field
    Equation 6 in PVP
    """
    return (w/(N*bias)) * (ng-alpha*nr)

def Fk(F):
    """
    Fourier transform
    """
    return np.fft.rfftn(F)

def Pk()
    """
    TODO: Requires summation in k-space
    """
    # Correct for shot noise
    Pk -= Pshot
    return

def Vk(k,dk):
    """
    TODO: Not entirely sure about this one
    """
    return 4*pi*(k**2)*(delta k)
    
def error(Pk, Vk):
    """
    Error bars
    """
    var = np.zeros(len(Pk))/Vk
    sig = np.sqrt(Pk**2*var) # 1-sigma? Not positive heres
