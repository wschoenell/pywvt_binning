'''
Created on Aug 2, 2013

@author: william
'''
from __future__ import division

from pywvt_binning.add_sn import add_noise, add_signal


import numpy as np

def wvt_binning(x, y, pixelSize, targetSN, signal, noise, cts=None, center=None,
                max_area=None, dens=None, keepfixed=None, gersho=False):
#   ; This is the main program that has to be called from external programs. 
#   ; It simply calls in sequence the different steps of the algorithms
#   ; and optionally plots the results at the end of the calculation.
    
    x = np.array(x, dtype=np.float) #FIXME:
    y = np.array(y, dtype=np.float) #FIXME:
    
    if keepfixed is not None:
        raise NotImplemented("KeepFixed is not implemented yet!")
    if gersho:
        raise NotImplemented("GerSHO is not implemented yet!")
    
    if dens is None:
        print 'WARNING: dens was not defined on wvt_binning call.'
        dens = signal / np.sqrt( np.power(noise,2) )
   
    print 'Bin-accretion...'
    
    binnumber, SNBin = wvt_bin_accretion(x, y, signal, noise, dens, targetSN, pixelSize, cts, max_area, keepfixed, center)
    
    print 'Reassign bad bins...'
    xNode, yNode, binnumber = wvt_reassign_bad_bins(x, y, dens, targetSN, binnumber, SNBin)
    
    np.savetxt('reassign.out', binnumber, fmt='%i') 


    return True, True, True, True, True, True
      

def wvt_bin_accretion(x, y, signal, noise, dens, targetSN, pixelSize, cts=None, max_area=None,
                      keepfixed=None, center=None):
    
    
    x = np.array(x, dtype=np.float) #FIXME:
    y = np.array(y, dtype=np.float) #FIXME:
    
    if keepfixed is not None:
        raise NotImplemented("KeepFixed is not implemented yet!")
    
    print '...making neighbor list...'
    neighborlist = wvt_makeneighborlist(x, y, pixelSize)    
    
    nx = len(x)
    binnumber = np.ones(shape=nx, dtype=np.int) * -1 #; will contain the bin number of each given pixel
    good = np.zeros(shape=nx, dtype=np.bool) #; will contain 1 IF the bin has been accepted as good
    SNbin = np.zeros(shape=nx, dtype=np.float) #; will contain the S/N of the corresponding bin
    
    if center is not None:
        #; Start at the specified center, if this keyword is set
        currentBin = [np.argmin((x-center[0])**2+(y-center[1])**2)]
        print 'Bin accretion started at the specified center:'
    else:
        if cts is not None:
#         ; For X-ray images, you might have a special cases, where you can have 
#         ; an artificially high S/N for a pixel that actually has no signal in 
#         ; it, namely when one knows the background very accurately. Thus, start
#         ; with the highest S/N pixel that has at least one count in it.
            aux = dens.copy()
            aux[cts < 1] = 0
            currentBin = [np.argmax(aux)]
        else:
            currentBin = [np.argmax(dens)]
            
        print 'Bin accretion started at highest S/N value:'
        
    SN = dens[currentBin]
    
    startingx=x[currentBin]
    startingy=y[currentBin]
    print 'First bin starts at ', startingx, startingy
    
#   ; The first bin will be assigned BINNUMBER = 1
#   ; With N pixels there will be at most N bins

#   Here we are changing from first bin equals one to equals 0 
    
    unBinned = []
    
    goodold = 0
    areaold = np.float(nx)
    area = 1.0
    mass = 1.e-200
    totarea = 0.0
    xBar_binned = startingx
    yBar_binned = startingy
    
    #TODO: Implement keepfixed.
    if keepfixed is not None:
        nfixed = 1
        currentBin = [np.argmin((x-keepfixed[0,0])**2+(y-keepfixed[1,0])**2)]
    
#     ; Check if the max_area is set. If not, leave the area unrestricted.
    if max_area is None:
        max_area = np.float(nx)
        
    totdens = 0
    totgoodbins = 0
    
    for i in range(nx):
        totgoodbins += goodold
        
        if goodold:
            print('bin: %i | n_pixels:  %i | %3.2f percent done' % 
                  ( totgoodbins, areaold,  100*(1-len(unBinned)/nx) ) )
    

        areaold = 0
    
        binnumber[currentBin] = i #; Here currentBin is still made of one pixel
        xBar = x[currentBin]
        yBar = y[currentBin]    #; Centroid of one pixels


        # Had to define these two:
        exneighbors, newmember = [], []
        
        neighbors = []
        
        while True:
            
# TODO: Implement keepfixed            
#             ; If nfixed is set, then leave the first nfixed pixels what they are.
#             if i < nfixed:
#                 good[currentBin] = True
            
#                   ; Test if the bin is good enough already
#                   ; Better way to decrease the scatter around the target S/N later on:
#                   ; Use the average bin members' S/N to estimate the S/N when another 
#                   ; pixel is added to see if adding another pixels would actually 
#                   ; increase the scatter around the target S/N due to "overshooting". 
#                   ; Also stop accreting if the bin size reached "max_area"
        
            modtargetSN = targetSN-SN*(np.sqrt(1.+1./area)-1.)/2.
            
#             print 'ifbreak1', SN, modtargetSN, area, max_area #, neighbors
            if np.bitwise_or(SN >= modtargetSN, area >= max_area):
                good[currentBin] = True
                SNOld = SN
                SNbin[currentBin] = SN
#                 print currentBin #, neighbors
#                 raw_input('break 1')
                break
            
#       ; Find nearest neighbors of pixels in the current bin, i.e. 
#       ; pixels contiguous with the bin, that have a chance of being accreted
#       ; For speed, remember the neighbors from the last accretion step and 
#       ; simply add the new ones.
            neighbors, exneighbors = wvt_findneighbors(currentBin,neighborlist, exneighbors, newmember, neighbors)
            
#       ; Accretable neighbors are those that aren't already binned
            if len(neighbors) > 0:
                wh = np.where(binnumber[neighbors] == -1)[0]
            else:
                wh = []
            
#       ; Stop if there aren't any accretable neighbors
            if len(wh) == 0:
                if SN > 0.8*targetSN:
                    good[currentBin] = True
                    SNbin[currentBin]=SN
                break
            
#       ; Otherwise keep only the accretable ones
            neighbors=neighbors[wh]
            
#       ; Search only the neighbors to get the next pixel
            dist_ = (x[neighbors]-xBar)**2 + (y[neighbors]-yBar)**2
            k = np.argmin(dist_)
            minDist = dist_[k]
            
#       ; Remember the old verified neighbors and the new members of the bin
            newmember = neighbors[k]
            nextBin = np.append(currentBin,neighbors[k])
            roundness = wvt_bin_roundness(x[nextBin], y[nextBin], pixelSize)
            
#       ; Compute the S/N one would obtain by adding
#       ; the candidate pixel to the current bin
            SNOld = SN
            SN = add_signal(signal[nextBin])/add_noise(noise[nextBin])
                        
#       ; Test whether the candidate pixel is connected to the
#       ; current bin, whether the possible new bin is round enough
#       ; Relaxed constraint on roundness to accept more bins
            if roundness > 1.0:
                if SN > 0.8*targetSN:
                    good[currentBin] = True
                    SNbin[currentBin]=SN
                break
            
#       ; If all the above tests are negative then accept the candidate pixel,
#       ; add it to the current bin, and continue accreting pixels
            binnumber[neighbors[k]] = i
            currentBin = nextBin
            
#       ; Update the centroid of the current bin
#TODO: Implement KeepFixed
#             nfixed = len(keepfixed)
#             if i > nfixed:
            xBar, yBar, mass = wvt_addto_weighted_centroid(x[newmember], y[newmember], dens[newmember]**2, xBar, yBar, mass)
#             else:
#                 xBar=keepfixed[0,i-1] #TODO: Implement KeepFixed
#                 yBar=keepfixed[1,i-1] #TODO: Implement KeepFixed
            
#       ; Update the area of the bin
            area = area + 1
            
        goodold = good[currentBin[0]]
        areaold = area
        area = 1
        mass = 0
        
        unBinned = np.argwhere(binnumber == -1)
        Binned = np.argwhere(binnumber != -1)
        
        if len(unBinned) == 0:
            break #; Stop if all pixels are binned
        
        totarea = totarea + areaold
        
#     ; Find the closest unbinned pixel to the centroid of all
#     ; the binned pixels, and start a new bin from that pixel.
#         if i >= nfixed: #TODO: Implement KeepFixed
        dist_ = (x[unBinned]-xBar_binned)**2 + (y[unBinned]-yBar_binned)**2
        k = np.argmin(dist_)
#         print k, unBinned, #unBinned[k]
#             minDist = dist_[k]
#         else: #TODO: Implement KeepFixed
#             dist_ = (x[unBinned]-keepfixed[0,i])**2 + (y[unBinned]-keepfixed[1,i])**2 #; keep the fixed centers unaltered!
        
        currentBin = unBinned[k]    #; The bin is initially made of one pixel
        SN = add_signal(signal[currentBin])/add_noise(noise[currentBin])

    binnumber[~good] = -1
    
    return binnumber, SNbin


def wvt_makeneighborlist(x, y, pixelSize):
#       ; Create a list of neighbors for each pixel
    nx = len(x)
    if len(y) != nx: raise Exception('wvt_makeneighborlist: x and y dimensions are different!')

    uniquex, tmp_ix = np.unique(x, return_index=True)
    uniquey, tmp_iy = np.unique(y, return_index=True)
    
    nuniquex = len(uniquex)
    nuniquey = len(uniquey)
        
    mask = np.ones((nuniquex+2,nuniquey+2), dtype=np.int) * -1

    for i in range(nx):
        ix = np.argwhere(x[i] == uniquex)
        iy = np.argwhere(y[i] == uniquey)
        ix=ix[0]
        iy=iy[0]
        mask[ix,iy] = i
    
    neighborlist = np.zeros((nx,4), dtype=np.int)
    
    for i in range(nuniquex):
        for j in range(nuniquey):
            pix=mask[i,j]
            if pix >= 0:
                neighborlist[pix,:] = np.sort([mask[i-1,j],mask[i,j+1],mask[i+1,j],mask[i,j-1]])[::-1]
    
    return neighborlist


def wvt_assign_to_bin(x, y, xnode, ynode, SNnode):
    # ; Assigns each pixel to the S/N weighted closest pixel 
    # ; i.e. this constructs the weighted voronoi tesselation
    return np.argmin( ((x-xnode)**2+(y-ynode)**2)*SNnode )

def wvt_assign_neighbors(neighborlist, binnumber):
    # ; Finds the binnumber of all neighbors of a bin and returns a list of adjacent
    # ; bins. You have to take into account that there might be no neighbors (-1)
    return np.append([-1], binnumber)[neighborlist]

def wvt_find_binneighbors(binneighbors, neighborlist, neighborbinnumber, binnumber, 
                          nbins, area=None):
#   ; Produces the final list of *unique* bin neighbors.
#   ; binneighbors will have the same format as the the REVERSE_INDICES
#   ; output of the HISTOGRAM function

    if area is None:
        area = np.histogram(binnumber, nbins, range=(0, nbins-1))
        npix = len(binnumber)
        binneighbors = np.array([nbins+1])
    
def wvt_findneighbors(group, neighborlist, exneighbors, newmember, neighbors):
#     ; Given the neighborlist for each pixel in the image, and a group of test
#     ; pixels, return a list of neighbors of the group's members that aren't
#     ; already group members themselves.
    ngroup = len(group)
    
    if np.bitwise_or(ngroup <= 2, len(exneighbors) < 1):
        neighbors  = np.unique(neighborlist[group,:])
        subneighbor = np.where(np.in1d(neighbors, group))[0] # np.where(np.ind1d()) is the same as MATCH on IDL
        subgroup = np.where(np.in1d(group, neighbors))[0][::-1] #also for the matching to the MATCH IDL routine

    
        if len(subneighbor) != 0:
            neighbors[subneighbor] = -1
        wh = np.where(neighbors != -1)[0]
        
        if len(wh) == 0:
            neighbors = np.array([])
        else:
            neighbors = neighbors[wh]
    
        exneighbors = np.zeros(neighborlist.shape[0], dtype=np.int)
        
        if len(neighbors != 0):
            exneighbors[neighbors] = 1
            
    else:
        tmp = neighborlist[newmember,:]
        wh = np.where( tmp != -1 )[0]
        
        if len(wh) > 0:
            tmp = tmp[wh]
            subtmp = np.where(exneighbors[tmp] == 0)[0]
            if len(subtmp) > 0:
                tmp = tmp[subtmp]
                exneighbors[tmp] = 1
                neighbors = np.append(neighbors, tmp)

#     ; We already know that the given list of neighbors is unique, since
#     ; we use the output from the last run. So there is no need to match them
#     ; again with everything.
#     ; First check which neighbors of the new members are already part of 
#     ; the group. All group members have an index of -1.

        exneighbors[newmember]=-1
        subneighbors = np.where(neighbors != newmember)
        if len(subneighbors) > 0:
            neighbors = neighbors[subneighbors]
        else:
            neighbors = None
            
    return neighbors, exneighbors

def wvt_bin_roundness(x, y, pixelSize):
#   ; Returns the "roundness" of a bin, as defined by equation (5) of 
#   ; Cappellari & Copin (2003)

    n = np.float(len(x))
    equivalentRadius = np.sqrt(n/np.pi)*pixelSize
    xBar = np.sum(x)/n #; unweighted centroid here!
    yBar = np.sum(y)/n
    maxDistance = np.sqrt( np.max((x-xBar)**2 + (y-yBar)**2) )
    roundness = maxDistance/equivalentRadius - 1.0
    
    return roundness

def wvt_addto_weighted_centroid(x, y, xymass, xBold, yBold, massold):
#   ; For speed, this procedure computes the geometric center of a bin by adding 
#   ; a new list of xy values for an existing bin. Useful in bin accretion step.
    xymass = np.array([xymass,], dtype = np.float)
    xymass[xymass < 1e-200] = 1e-200 # Same as IDL maximum operator 
    
    mass = np.sum(xymass) + massold
    
    if mass > 0: #FIXME: Is this if useless? I think so...
        xBar = (np.sum(xymass*x) + massold*xBold)/mass
        yBar = (np.sum(xymass*y) + massold*yBold)/mass
    else:
        xBar = xBold
        yBar = yBold

    return xBar, yBar, mass

def wvt_reassign_bad_bins(x, y, dens, targetSN, binnumber, SNBin):
#   ; Find pixels that the bin accretion step wasn't able to assign to a bin and
#   ; reassign them to the next closest bin
#   ; Implements steps (vi)-(vii) in section 5.1 of Cappellari & Copin (2003)

    binnumber = wvt_renumber_binnumber(binnumber)
    
#     I've omitted these lines of code here beacuse this checking is done on wvt_renumber_binnumber
#     area = histogram(binnumber, REVERSE_INDICES=r, MIN=1)
#     good = where(area gt 0, nnodes) ; Obtain the index of the good bins

    nnodes = np.max(binnumber) + 1
    
    xNode, yNode = np.transpose([wvt_unweighted_centroid(x[index], y[index]) for index in 
                                 [np.where(binnumber == i)[0] for i in range(nnodes)]])
    
#   ; Reassign pixels to the closest centroid of a good bin
    bad = np.where(binnumber == -1)[0]
    good = range(len(xNode))
    
    print 'binnumber', binnumber
    
    for i in range(len(bad)):
        binnumber[bad[i]] = good[wvt_assign_to_bin(x[bad[i]],y[bad[i]],xNode, yNode, 1)]
    
    print binnumber[bad]
    
#   ; Recompute all centroids of the reassigned bins.
#   ; These will be used as starting point for the WVT.
    xNode, yNode = np.transpose([wvt_unweighted_centroid(x[index], y[index]) for index in 
                                 [np.where(binnumber == i)[0] for i in range(nnodes)]])
    
    
    return xNode, yNode, binnumber
    
def wvt_unweighted_centroid(x, y):
#   ; Computes the geometric center of one bin
  
    mass=len(x)
    xBar=np.sum(x)/mass
    yBar=np.sum(y)/mass    
    
    return xBar, yBar
    
    
def wvt_renumber_binnumber(binnumber, area = None):
#   ; Kicks out zero pixel bins and renumbers the rest continuously, starting at
#   ; the value minbinnumber
    
    if area is not None:
        raise NotImplemented("Area can not be anything but None on wvt_renumber_binnumber yet!")

    uniq = np.unique(binnumber[binnumber >= 0]) # Because binnumber = -1 is a null bin.
    
    for i in range(len(uniq)):
        binnumber[binnumber == uniq[i]] = i
    
    return binnumber