import numpy as np
from scipy.special import erfc

def BM(z):

    """
    Calibration parameter
    """
    BM0 = 0.0
    alpha = 0.0

    return BM0 * (1 + alpha)


def sigmaLnM(z):

    """
    Calibration parameter
    """
    beta = 0.0
    sLnM0 = 0.2 # random value

    return sLnM0 * (1 + z)**beta


def xm(z, M, Mthr, m):

    """

    """
    DlnM = 0.3
    #Mthr = 1e14 # what is this

    #return (np.log(Mthr) + (m * np.log(DlnM)) - BM(z) - np.log(M))/np.sqrt(2*sigmaLnM(z)**2)
    return (np.log(Mthr) + m * DlnM - BM(z) - np.log(M))/np.sqrt(2*sigmaLnM(z)**2)
