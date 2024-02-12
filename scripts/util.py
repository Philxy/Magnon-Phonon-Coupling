import numpy as np

def Power(a,b):
    return np.power(a,b)

def Sqrt(x):
    return np.sqrt(x)



def gaussian(x, mu, sig):
    return 1./(np.sqrt(2.*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2)