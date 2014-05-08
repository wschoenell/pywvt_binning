'''
Created on Aug 2, 2013

@author: william
'''

import numpy as np
from pywvt_binning.wvt_binning import wvt_binning


# def wvt_image(signal, noise, targetSN, weight, center=None, max_area=None,
#               mask=None, ctsimage=None, keepfixed=None, gersho=False):

def wvt_image(signal, noise, targetSN, center=None, max_area=None,
              mask=None, ctsimage=None, keepfixed=None, gersho=False):

#   ; Interface for using WVT_BINNING with images instead of pixel lists
#   ; Takes the images, converts them into pixel lists, stores all of their
#   ; properties (signal, noise, cts, ...) in the global variable P (structure) 
#   ; and passes everything on to WVT_BINNING. This was designed particularly 
#   ; with X-ray images in mind. The structure P will hold all information 
#   ; about the individual pixels (signal, noise, cts, ...)
    
    if mask is None:
        mask = np.ones_like(signal) #FIXME: Should this be dtype=np.bool as well?
    if ctsimage is None:
        ctsimage = mask.copy()
    
    ngood = np.sum(mask == 1)
    
    y, x = np.indices(signal.shape)
    x = x[mask != 0].ravel()
    y = y[mask != 0].ravel()
    
    
#     ; Check if the pixel should be included or not
    
    p_signal = np.ravel(signal[mask != 0])
    p_noise = np.ravel(noise[mask != 0])
    p_cts = np.ravel(ctsimage[mask != 0])
    
    dens = p_signal/p_noise # If noise is == 0 we should have problems!
    
    pixelSize=1.0
    
    if keepfixed is not None:
        raise NotImplemented("KeepFixed is not implemented yet!")
    if gersho:
        raise NotImplemented("GerSHO is not implemented yet!")
    
#       ; Start the main binning algorithm

#     binnumber, xNode, yNode, SNbin, area, binValue = wvt_binning(x, y, pixelSize, targetSN,
    binnumber = wvt_binning(x, y, pixelSize, targetSN, 
                            p_signal, p_noise, p_cts, center,
                            max_area, dens, keepfixed, gersho)
    
    binned_image = np.zeros(signal.shape[0] * signal.shape[1], dtype=np.int64)
    binned_image[np.ravel(mask)] = binnumber
    binned_image.reshape(signal.shape)
    

#     return binned_image, xNode, yNode, SNbin, binnumber, area, binValue    
    return binned_image # binnumber