"""
utils
======

Miscellaneous utilities used in the pyMRMS package.

    dt64toDT

"""
from datetime import datetime as dt

def np64toDT(dt64):
    """
    Convert a numpy datetime64 value into a datetime.
    
    
    Parameters
    ----------
    dt64 : numpy.datetime64
        Single datetime64 value to be converted to a datetime
        
    
    Returns
    -------
    dtOut : datetime


    """
    dtOut = dt.utcfromtimestamp(int((dt64.astype(dt))*1e-9))
    
    return dtOut
    