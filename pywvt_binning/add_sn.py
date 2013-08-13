'''
Created on Aug 5, 2013

@author: william
'''

import numpy as np


def add_signal(signal, cts=None):
#   ; Adds the signal inside the bin, specified by the indeces of the bin members
    if cts is None:
        r = np.sum(signal) / len(signal)
    elif np.sum(cts > 0):
        r = np.sum(signal) / len(signal)
    else:
        r = 0.0

    return r


def add_noise(noise):
#  ; Adds the noise inside the bin, specified by the indeces of the bin members
    return np.sqrt( np.sum(np.power(noise,2)) ) / len(noise)