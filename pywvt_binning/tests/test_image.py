'''
Created on Aug 2, 2013

@author: william
'''

import pyfits
import numpy as np
import matplotlib.pyplot as plt
from pywvt_binning.wvt_image import wvt_image

fits_path = '../../contrib/idl/adaptive_binning/test_wvt_image/'

ctsimage = pyfits.open('%s/cts_256.fits' % fits_path)[0].data # ; Raw Counts
image = pyfits.open('%s/img_256.fits' % fits_path)[0].data # ; Fluxed image
noise = np.sqrt(pyfits.open('%s/var_256.fits' % fits_path)[0].data) # ; Associated noise in fluxed image
mask = pyfits.open('%s/mask_256.fits' % fits_path)[0].data # ; The mask, "1" for good, "0" bad

center=[127.,127.]

targetSN = 10.0

max_area = 100**2

binned_image, xNode, yNode, SNbin, binnumber, area, binValue = wvt_image(image, noise, targetSN, center=center, max_area=max_area)

plt.figure(1)
plt.clf()
plt.imshow(binnumber)
plt.figure(2)
plt.clf()
plt.imshow(image)


# wvt_image(image,noise,targetSN,binnedimage,xnode,ynode,snbin=snbin,
#           mask=mask, ctsimage=ctsimage, binnumber=binnumber, binvalue=binvalue, 
#           center=center, save_all=save_all, max_area=max_area)

# wvt_bin_accretion(x, y, dens, targetSN, binnumber, pixelSize)