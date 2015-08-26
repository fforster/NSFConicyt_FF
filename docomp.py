# F.Forster
# compare LSST and HiTS subtractions
# needs HiTS difference and inverse variance maps, catalogues, and astrometric solution files

import numpy as np
import os
import subprocess
import re
from pylab import *
import pyfits as fits
from scipy import stats



doplothisto = False
doplotstamps = False

basedir = '/home/fforster/Work/LSST/COMP'
lsstdir = '%s/LSST' % basedir
hitsdir = '%s/HiTS' % basedir 
outdir = '%s/plots' % basedir#/nocovariance'

field = 'Blind14A_10'
epoch1 = 2
EXP1 = 288976
epoch2 = 12
EXP2 = 289493
epoch2 = 19
EXP2 = 289820

lsstfits = sort(os.listdir(lsstdir))
lsstfiles = {}
hitsfits = sort(os.listdir(hitsdir))

# read lsst files and get file corresponding to given CCD
print "LSST files (file, EXPNUM1, EXPNUM2, CCD):"
for filei in lsstfits:
    if filei[-4:] != "fits":
        continue
    (EXPNUM1, EXPNUM2, nCCD) = re.findall("sub\_(.*?)\_(.*?)\_(.*?)\.fits", filei)[0]
    
    if int(EXPNUM1) != EXP1 or int(EXPNUM2) != EXP2:
        continue

    # check CCD name
    keys = ['EXTNAME']
    for key in keys:
        command = "listhead %s/%s | grep -i %s" % (lsstdir, filei, key)
        output = subprocess.check_output(command, shell=True)
        CCD = re.findall("EXTNAME\s=\s\'(.*?)\s+\'", output)[0]
        print filei, EXPNUM1, EXPNUM2, CCD
    lsstfiles[CCD] = filei

print "HiTS files (file, epoch1, epoch2, CCD):"
for filei in hitsfits:

    if filei[:4] != "Diff":
        continue

    CCDread = re.findall("Diff_%s_(.*?)_%02it-%02i_grid02.lanczos2.fits" % (field, epoch2, epoch1), filei)
    if CCDread == []:
        CCDread = re.findall("Diff_%s_(.*?)_%02i-%02it_grid02.lanczos2.fits" % (field, epoch2, epoch1), filei)
        if CCDread == []:
            print "WARNING: Cannot find CCD or inconsistent epoch in HiTS file name %s" % filei
            continue#            sys.exit()
    CCD = CCDread[0]
#    if CCD != "S25":
#        continue
    print filei, epoch1, epoch2, CCD

    
    lsstfilename = "%s/%s" % (lsstdir, lsstfiles[CCD])
    hitsdifffilename = "%s/%s" % (hitsdir, filei)
    hitsinvVARfilename = "%s/%s" % (hitsdir, filei.replace("Diff", "invVAR"))

    try:
        hitsdiff = fits.open(hitsdifffilename)
        hitsinvVAR = fits.open(hitsinvVARfilename)
        lsstfits = fits.open(lsstfilename)
    except:
        print "WARNING: cannot open fits files"
        sys.exit()

    datahitsdiff = hitsdiff[0].data
    e_datahitsdiff = np.sqrt(1. / hitsinvVAR[0].data)
    datalsstdiff = lsstfits[1].data
    e_datalsstdiff = np.sqrt(lsstfits[3].data)

    (nxhits, nyhits) = np.shape(datahitsdiff)
    (nxlsst, nylsst) = np.shape(datalsstdiff)

    hits = datahitsdiff#[400:1000, 400:1000] 
    lsst = datalsstdiff#[400:1000, 400:1000] 

    e_hits = e_datahitsdiff#[400:1000, 400:1000] 
    e_lsst = e_datalsstdiff#[400:1000, 400:1000] 

    hitsSNR = hits / e_hits
    lsstSNR = lsst / e_lsst

    if doplothisto:
    
        # differences
        fig, ax = plt.subplots(ncols = 2, figsize = (12, 6))
        ax[0].hist(hits[(hits > -150) & (hits < 150)], log = True, alpha = 0.5, color = 'b', bins = 50, label = 'HiTS', lw = 0)
        ax[0].hist(lsst[(lsst > -150) & (lsst < 150)], log = True, alpha = 0.5, color = 'r', bins = 50, label = 'LSST', lw = 0)
        ax[0].legend(loc = 1)
        ax[0].set_xlabel("SNR")
        ax[0].set_ylabel("N")
        bins = np.logspace(-1, 6, 100)
        ax[1].hist(np.abs(hits[np.isfinite(hits)]), log = True, alpha = 0.5, color = 'b', bins = bins, label = 'HiTS', lw = 0)
        ax[1].hist(np.abs(lsst[np.isfinite(lsst)]), log = True, alpha = 0.5, color = 'r', bins = bins, label = 'LSST', lw = 0)
        ax[1].set_xscale('log')
        ax[1].legend(loc = 1)
        ax[1].set_xlabel("SNR")
        ax[1].set_ylabel("N")
        plt.savefig("%s/histodiff_%s_%s.png" % (outdir, field, CCD))
        
        
        # errors
        del fig, ax
        fig, ax = plt.subplots(ncols = 2, figsize = (12, 6))
        ax[0].hist(e_hits[(e_hits > 0) & (e_hits < 70)], log = True, alpha = 0.5, color = 'b', bins = 50, label = 'HiTS', lw = 0)
        ax[0].hist(e_lsst[(e_lsst > 0) & (e_lsst < 70)], log = True, alpha = 0.5, color = 'r', bins = 50, label = 'LSST', lw = 0)
        ax[0].legend(loc = 1)
        ax[0].set_xlabel("SNR")
        ax[0].set_xlabel("N")
        bins = np.logspace(-1, 6, 100)
        ax[1].hist(np.abs(e_hits[np.isfinite(e_hits)]), log = True, alpha = 0.5, color = 'b', bins = bins, label = 'HiTS', lw = 0)
        ax[1].hist(np.abs(e_lsst[np.isfinite(e_lsst)]), log = True, alpha = 0.5, color = 'r', bins = bins, label = 'LSST', lw = 0)
        ax[1].set_xscale('log')
        ax[1].legend(loc = 1)
        ax[1].set_xlabel("SNR")
        ax[1].set_ylabel("N")
        plt.savefig("%s/histoerror_%s_%s.png" % (outdir, field, CCD))
        
        # SNR
        del fig, ax
        fig, ax = plt.subplots(ncols = 2, figsize = (12, 6))
        ax[0].hist(hitsSNR[(hitsSNR > -10) & (hitsSNR < 10)], log = True, alpha = 0.5, color = 'b', bins = 50, label = 'HiTS', normed = True, lw = 0)
        ax[0].hist(lsstSNR[(lsstSNR > -10) & (lsstSNR < 10)], log = True, alpha = 0.5, color = 'r', bins = 50, label = 'LSST', normed = True, lw = 0)
        bins = np.linspace(-10, 10, 100)
        ax[0].set_ylim(1e-8, 1)
        ax[0].plot(bins, stats.norm.pdf(bins), c = 'k', label = "G(0, 1)")
        ax[0].legend(loc = 1)
        ax[0].set_xlabel("SNR")
        ax[0].set_ylabel("pdf")
        bins = np.logspace(-1, 6, 100)
        ax[1].hist(np.abs(hitsSNR[np.isfinite(hitsSNR)]), log = True, alpha = 0.5, color = 'b', bins = bins, label = 'HiTS', lw = 0)
        ax[1].hist(np.abs(lsstSNR[np.isfinite(lsstSNR)]), log = True, alpha = 0.5, color = 'r', bins = bins, label = 'LSST', lw = 0)
        ax[1].set_xscale('log')
        ax[1].legend(loc = 1)
        ax[1].set_xlabel("SNR")
        ax[1].set_ylabel("pdf")
        plt.savefig("%s/histoSNR_%s_%s.png" % (outdir, field, CCD))
    

    # open catalogues
    
    cat1 = np.loadtxt("%s/%s_%s_%02i_image_crblaster.fits-catalogue_wtmap_backsize64.dat" % (hitsdir, field, CCD, epoch1)).transpose()
    cat2 = np.loadtxt("%s/%s_%s_%02i_image_crblaster.fits-catalogue_wtmap_backsize64.dat" % (hitsdir, field, CCD, epoch2)).transpose()

    # open transformation
    match = np.load("%s/match_%s_%s_%02i-%02i.npy" % (hitsdir, field, CCD, epoch2, epoch1))

    aflux, e_aflux, rms, order = match[0:4]
    sol_astrometry = match[4:]

    # find stars in first catalogue
    x1 = cat1[1]
    y1 = cat1[2]
    z1 = cat1[5]
    e_z1 = cat1[6]

    # find stars in 2nd catalogue
    x2 = cat2[1]
    y2 = cat2[2]
    z2 = cat2[5]
    e_z2 = cat2[6]

    def applytransformation(order, x1, y1, sol):

        # this is slow, but I prefer fewer bugs than speed at the moment...                                                                               
        x1t = sol[0] + sol[2] * x1 + sol[3] * y1
        y1t = sol[1] + sol[4] * x1 + sol[5] * y1
        if order > 1:
            x1t = x1t + sol[6] * x1 * x1 + sol[7] * x1 * y1 + sol[8] * y1 * y1
            y1t = y1t + sol[9] * x1 * x1 + sol[10] * x1 * y1 + sol[11] * y1 * y1
        if order > 2:
            x1t = x1t + sol[12] * x1 * x1 * x1 + sol[13] * x1 * x1 * y1 + sol[14] * x1 * y1 * y1 + sol[15] * y1 * y1 * y1
            y1t = y1t + sol[16] * x1 * x1 * x1 + sol[17] * x1 * x1 * y1 + sol[18] * x1 * y1 * y1 + sol[19] * y1 * y1 * y1

        return (x1t, y1t)

    (x1t, y1t) = applytransformation(order, x1, y1, sol_astrometry)

    # compute residuals
    x2m = np.zeros(len(x1))
    y2m = np.zeros(len(x1))
    for i in range(len(x1)):
        dist = np.sqrt((x1t[i] - x2)**2 + (y1t[i] - y2)**2)
        idx = np.argmin(dist)
        x2m[i] = x2[idx]
        y2m[i] = y2[idx]

    # extract limits for stars in 1st catalogue
    nstamp = 10
    nstamph = nstamp / 2
    stampshits = np.zeros((len(x1), nstamp, nstamp))
    stampslsst = np.zeros((len(x1), nstamp, nstamp))
    stampshitsSNR = np.zeros((len(x1), nstamp, nstamp))
    stampslsstSNR = np.zeros((len(x1), nstamp, nstamp))
    maskstars = np.ones(len(x1), dtype = bool)

    # save stamps
    for i in range(len(x1)):
        # star must be well inside image
        if x1[i] > nstamp and x1[i] < nxhits - 2 * nstamp and y1[i] > nstamp and y1[i] < nxhits - 2 * nstamp and x1t[i] > 2 * nstamp and x1t[i] < nylsst - 2 * nstamp and y1t[i] > 2 * nstamp and y1t[i] < nxlsst - 2 * nstamp:
            
            #print x1[i], y1[i], x1t[1], y1t[i]
            
            yloc = np.round(y1[i])
            xloc = np.round(x1[i])
            yloct = np.round(y1t[i]) - 1 # for LSST comparison (not sure why)
            xloct = np.round(x1t[i]) - 1 # "
            
            # read the arrays
            stamphits = hits[yloc - nstamph: yloc + nstamph, xloc - nstamph: xloc + nstamph]
            stamplsst = lsst[yloct - nstamph: yloct + nstamph, xloct - nstamph: xloct + nstamph]
            stamphitsSNR = hitsSNR[yloc - nstamph: yloc + nstamph, xloc - nstamph: xloc + nstamph]
            stamplsstSNR = lsstSNR[yloct - nstamph: yloct + nstamph, xloct - nstamph: xloct + nstamph]

            # check that the median is not exactly zero and the fraction of negative pixels is between 30 and 70% of the total
            if np.median(stamphits) != 0 and np.median(stamplsst) != 0 and np.isfinite(np.sum(stamphits)) and np.isfinite(np.sum(stamplsst)) \
               and np.sum(stamphits < 0) > 0.3 * nstamp**2 and np.sum(stamphits < 0) < 0.7 * nstamp**2 \
               and np.sum(stamplsst < 0) > 0.3 * nstamp**2 and np.sum(stamplsst < 0) < 0.7 * nstamp**2:
                stampshits[i, :, :] = stamphits
                stampslsst[i, :, :] = stamplsst
                stampshitsSNR[i, :, :] = stamphitsSNR
                stampslsstSNR[i, :, :] = stamplsstSNR
            else:
                #print "ZERO", x1[i], y1[i]
                maskstars[i] = False
        else:
            #print "OUTSIDE", x1[i], y1[i], x1t[i], y1t[i]
            maskstars[i] = False

    x1 = x1[maskstars]
    y1 = y1[maskstars]
    x2m = x2m[maskstars]
    y2m = y2m[maskstars]
    x1t = x1t[maskstars]
    y1t = y1t[maskstars]
    z1 = z1[maskstars]
    stampshits = stampshits[maskstars]
    stampslsst = stampslsst[maskstars]
    stampshitsSNR = stampshitsSNR[maskstars]
    stampslsstSNR = stampslsstSNR[maskstars]

    idxsort = np.argsort(z1)
    x1 = x1[idxsort[::-1]]
    y1 = y1[idxsort[::-1]]
    x2m = x2m[idxsort[::-1]]
    y2m = y2m[idxsort[::-1]]
    x1t = x1t[idxsort[::-1]]
    y1t = y1t[idxsort[::-1]]
    z1 = z1[idxsort[::-1]]
    stampshits = stampshits[idxsort[::-1]]
    stampslsst = stampslsst[idxsort[::-1]]
    stampshitsSNR = stampshitsSNR[idxsort[::-1]]
    stampslsstSNR = stampslsstSNR[idxsort[::-1]]

    # check residual
    fig, ax = plt.subplots()
    match = np.array((np.sqrt((x2m - x1t)**2 + (y2m - y1t)**2) < 2) & (z1 > 200) & (z1 < 1e5) & (x1 > 100) & (x1 < nyhits - 100) & (y1 > 100) & (y1 < nxhits - 100))
    ax.scatter(x2m[match] - x1t[match], y2m[match] - y1t[match], marker = '.', lw = 0, alpha = 0.5)
    plt.savefig("%s/residuals.png" % outdir)

    # mesh of x and y values for center of mass computation
    X, Y = np.meshgrid(np.array(range(nstamp)), np.array(range(nstamp)))
    R = np.sqrt((X - nstamph)**2 + (Y - nstamph)**2)
    xneghits = np.ones(len(x1), dtype = float)
    xposhits = np.ones(len(x1), dtype = float)
    xneglsst = np.ones(len(x1), dtype = float)
    xposlsst = np.ones(len(x1), dtype = float)
    yneghits = np.ones(len(x1), dtype = float)
    yposhits = np.ones(len(x1), dtype = float)
    yneglsst = np.ones(len(x1), dtype = float)
    yposlsst = np.ones(len(x1), dtype = float)

    # plot stamps
    ncols = 15
    nrows = 7
    if doplotstamps:
        del fig, ax
        fig, ax = plt.subplots(ncols = ncols, nrows = 2 * nrows, figsize = (13.5, 12))
    for i in range(len(x1)):

        # compute centers of mass
        pos = (stampshits[i] > 0) & (R <= nstamph)
        neg = (stampshits[i] < 0) & (R <= nstamph)
        xposhits[i] = np.sum(X[pos] * stampshits[i][pos]) / np.sum(stampshits[i][pos])
        yposhits[i] = np.sum(Y[pos] * stampshits[i][pos]) / np.sum(stampshits[i][pos])
        xneghits[i] = np.sum(X[neg] * stampshits[i][neg]) / np.sum(stampshits[i][neg])
        yneghits[i] = np.sum(Y[neg] * stampshits[i][neg]) / np.sum(stampshits[i][neg])
        pos = (stampslsst[i] > 0)
        neg = (stampslsst[i] < 0)
        xposlsst[i] = np.sum(X[pos] * stampslsst[i][pos]) / np.sum(stampslsst[i][pos])
        yposlsst[i] = np.sum(Y[pos] * stampslsst[i][pos]) / np.sum(stampslsst[i][pos])
        xneglsst[i] = np.sum(X[neg] * stampslsst[i][neg]) / np.sum(stampslsst[i][neg])
        yneglsst[i] = np.sum(Y[neg] * stampslsst[i][neg]) / np.sum(stampslsst[i][neg])

        if i >= ncols * nrows:
            continue
        if doplotstamps:
            j = i / ncols
            k = mod(i, ncols)

            # plot stamps
            ax[2 * j, k].imshow(stampshits[i], interpolation = 'nearest', cmap = 'gray')
            ax[2 * j + 1, k].imshow(stampslsst[i], interpolation = 'nearest', cmap = 'gray')
            
            # plot center of mass
            ax[2 * j, k].scatter(xposhits[i], yposhits[i], marker = 'x', c = 'r', lw = 2)
            ax[2 * j, k].scatter(xneghits[i], yneghits[i], marker = 'x', c = 'b', lw = 2)
            ax[2 * j, k].plot([xneghits[i], xposhits[i]], [yneghits[i], yposhits[i]], c = 'y', lw = 2)
            ax[2 * j + 1, k].scatter(xposlsst[i], yposlsst[i], marker = 'x', c = 'r', lw = 2)
            ax[2 * j + 1, k].scatter(xneglsst[i], yneglsst[i], marker = 'x', c = 'b', lw = 2)
            ax[2 * j + 1, k].plot([xneglsst[i], xposlsst[i]], [yneglsst[i], yposlsst[i]], c = 'y', lw = 2)
        
            # fix axis
            ax[2 * j, k].set_xlim(-0.5, nstamp - 0.5)
            ax[2 * j, k].set_ylim(-0.5, nstamp - 0.5)
            ax[2 * j + 1, k].set_xlim(-0.5, nstamp - 0.5)
            ax[2 * j + 1, k].set_ylim(-0.5, nstamp - 0.5)
            ax[2 * j, k].spines['top'].set_color('red')
            ax[2 * j, k].axes.get_xaxis().set_visible(False)
            ax[2 * j, k].axes.get_yaxis().set_visible(False)
            ax[2 * j + 1, k].axes.get_xaxis().set_visible(False)
            ax[2 * j + 1, k].axes.get_yaxis().set_visible(False)
            ax[2 * j + 1, k].spines['bottom'].set_color('red')


    if doplotstamps:
        fig.subplots_adjust(wspace = 0, hspace = 0)
        plt.savefig("%s/stamps_%s_%s_%02i-%02i_hits.png" % (outdir, field, CCD, epoch2, epoch1), pad_inches = 0, bbox_inches = 'tight', hspace = 0, vspace = 0)

        
    dipolelim = 1.5
    fluxlim = 1e3

    # plot source flux vs center of masses distance
    del fig, ax
    fig, ax = plt.subplots()
    dipolestrengthhits = np.sqrt(np.sum(stampshits[:, R <= nstamph]**2, axis = 1))
    dipolestrengthlsst = np.sqrt(np.sum(stampslsst[:, R <= nstamph]**2, axis = 1))
    dipolelengthhits = np.sqrt((xneghits - xposhits)**2 + (yneghits - yposhits)**2)
    dipolelengthlsst = np.sqrt((xneglsst - xposlsst)**2 + (yneglsst - yposlsst)**2)
    dipoleanglehits = 180. / np.pi * np.arctan2(yposhits - yneghits, xposhits - xneghits)
    dipoleanglelsst = 180. / np.pi * np.arctan2(yposlsst - yneglsst, xposlsst - xneglsst)
    residualanglehits = 180. / np.pi * np.arctan2(y1t - y2m, x1t - x2m)
    residual = np.sqrt((y2m - y1t)**2 + (x2m - x1t)**2)
    ax.scatter(dipolelengthhits[z1 > 0], z1[z1 > 0], marker = '.', c = 'b', alpha = 0.5, lw = 0)
    ax.scatter(dipolelengthlsst[z1 > 0], z1[z1 > 0], marker = '.', c = 'r', alpha = 0.5, lw = 0)
    ax.axvline(dipolelim)
    ax.axhline(fluxlim)
    ax.set_yscale('log')
    ax.set_xlabel("Dipole length [pix]")
    ax.set_ylabel("sextractor flux [ADU]")
    plt.savefig("%s/dipoleflux_%s_%s_%02i-%02i.png" % (outdir, field, CCD, epoch2, epoch1))

    # dipole field
    del fig, ax
    fig, ax = plt.subplots(figsize = (15, 7))
    mask = np.array(((dipolelengthhits > dipolelim) | (dipolelengthlsst > dipolelim)) & (z1 > fluxlim) & match & (residual > 0.2))

    # check residual
    ax.scatter(y1[match & mask], x1[match & mask], marker = '.', alpha = 0.5, c = 'k', s = 200)
    ax.quiver(y1[mask & match], x1[mask & match],  y2m[mask & match] - y1t[mask & match], x2m[mask & match] - x1t[mask & match], pivot='middle', headwidth=4, headlength=6, lw = 0.1, alpha = 0.5, color = 'k', label = 'astrom. res.')
    ax.quiver(y1[mask], x1[mask], yposhits[mask] - yneghits[mask], xposhits[mask] - xneghits[mask], pivot='middle', headwidth=4, headlength=6, lw = 0, alpha = 0.5, color = 'b', label = 'HiTS')
    ax.quiver(y1[mask], x1[mask], yposlsst[mask] - yneglsst[mask], xposlsst[mask] - xneglsst[mask], pivot='middle', headwidth=4, headlength=6, lw = 0, alpha = 0.5, color = 'r', label = 'LSST')
    ax.legend(loc = 1)
    ax.set_xlabel("j [pix]")
    ax.set_ylabel("i [pix]")
    ax.set_xlim(0, nxhits)
    ax.set_ylim(0, nyhits)
    ax.set_title("Dipole field (blue: HiTS, red: LSST, gray: HiTS astrometric residual, length: distance between positive and negative dipoles)", fontsize = 10)
    plt.savefig("%s/dipolefield_%s_%s_%02i-%02i.png" % (outdir, field, CCD, epoch2, epoch1), bbox_inches = 'tight')

    # dipole sqrt(squared sum of flux)
    del fig, ax
    fig, ax = plt.subplots()
    mask = np.array(((dipolelengthhits > dipolelim) | (dipolelengthlsst > dipolelim)))
    ax.scatter(z1[mask], dipolestrengthhits[mask], c = 'b', label = 'HiTS', lw = 0)
    ax.scatter(z1[mask], dipolestrengthlsst[mask], c = 'r', label = 'LSST', lw = 0)
    ax.legend(loc = 2)
    ax.set_xlabel("Star flux [ADU]")
    ax.set_ylabel("sqrt(sum(dipole_(r<=5)^2))")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(5e1, 2e4)
    plt.savefig("%s/dipolestrength_%s_%s_%02i-%02i.png" % (outdir, field, CCD, epoch2, epoch1), bbox_inches = 'tight')

    # dipole sqrt(squared sum of flux)
    del fig, ax
    fig, ax = plt.subplots()
    mask = np.array(z1 > 0)
    ax.scatter(dipolestrengthhits[mask], dipolelengthhits[mask], c = 'b', label = 'HiTS', lw = 0, alpha = 0.5)
    ax.scatter(dipolestrengthlsst[mask], dipolelengthlsst[mask], c = 'r', label = 'LSST', lw = 0, alpha = 0.5)
    ax.legend(loc = 2)
    ax.set_xlabel("sqrt(sum(dipole_(r<=5)^2))")
    ax.set_ylabel("Dipole length [pix]")
    ax.set_xscale('log')
    plt.savefig("%s/dipolestrengthsize_%s_%s_%02i-%02i.png" % (outdir, field, CCD, epoch2, epoch1), bbox_inches = 'tight')

    # dipole strength histograms
    del fig, ax
    fig, ax = plt.subplots()
    mask = np.array((z1 > 0) & np.isfinite(dipolestrengthlsst))
    ax.hist(dipolestrengthhits[mask], color = 'b', label = 'HiTS', alpha = 0.5, log = True, bins = np.logspace(1, 4, 100))
    ax.hist(dipolestrengthlsst[mask], color = 'r', label = 'LSST', alpha = 0.5, log = True, bins = np.logspace(1, 4, 100))
    ax.legend(loc = 1)
    ax.set_xlabel("sqrt(sum(dipole_(r<=5)^2))")
    ax.set_ylabel("N")
    ax.set_xscale('log')
    ax.set_xlim(2e1, 1e4)
    plt.savefig("%s/dipolestrength_hist_%s_%s_%02i-%02i.png" % (outdir, field, CCD, epoch2, epoch1), bbox_inches = 'tight')

    # dipole angles
    del fig, ax
    fig, ax = plt.subplots()
    mask = np.array((z1 > fluxlim))
    ax.scatter(dipolelengthhits[mask], dipoleanglehits[mask], marker = '.', c = 'r', alpha = 0.5)
    ax.scatter(dipolelengthlsst[mask], dipoleanglelsst[mask], marker = '.', c = 'r', alpha = 0.5)
    ax.set_xlabel("Dipole length [pix]")
    ax.set_ylabel("Dipole angle [deg]")
    plt.savefig("%s/dipoleangle_%s_%s_%02i-%02i.png" % (outdir, field, CCD, epoch2, epoch1), bbox_inches = 'tight')

    # dipole angle differences
    del fig, ax
    fig, ax = plt.subplots()
    mask = np.array(((dipolelengthhits > dipolelim) | (dipolelengthlsst > dipolelim)) & (z1 > fluxlim) & (residual > 0.05))
    ax.hist(np.mod(dipoleanglehits[mask] - dipoleanglelsst[mask] + 360, 360), lw = 0, normed = True, label = 'HiTS - LSST', alpha = 0.5)
    ax.hist(np.mod(dipoleanglehits[mask & match] - residualanglehits[mask & match] + 360, 360), lw = 0, normed = True, label = 'HiTS - astrom.res.(HiTS)', alpha = 0.5)
    ax.set_xlabel("Dipole angle differences [deg]")
    ax.set_ylabel("N")
    ax.legend(loc = 1, fontsize = 9)
    ax.set_xlim(0, 360)
    plt.savefig("%s/dipoleanglediff_%s_%s_%02i-%02i.png" % (outdir, field, CCD, epoch2, epoch1), bbox_inches = 'tight')

    # dipole angle histogram
    del fig, ax
    fig, ax = plt.subplots()
    mask = np.array(((dipolelengthhits > dipolelim) | (dipolelengthlsst > dipolelim)) & (z1 > fluxlim))
    ax.hist(np.mod(dipoleanglehits[mask] + 180, 180), alpha = 0.5, color = 'b', label = 'HiTS')
    ax.hist(np.mod(dipoleanglelsst[mask] + 180, 180), alpha = 0.5, color = 'r', label = 'LSST')
    ax.legend(loc = 2)
    ax.set_xlabel("mod(Dipole angle + 180 deg, 180 deg) [deg]")
    ax.set_ylabel("N")
    plt.savefig("%s/dipoleangle_hist_%s_%s_%02i-%02i.png" % (outdir, field, CCD, epoch2, epoch1), bbox_inches = 'tight')

    # dipole lengths vs residual length
    del fig, ax
    fig, ax = plt.subplots()
    mask = np.array((z1 > fluxlim))
    ax.scatter(dipolelengthhits[match & mask], residual[match & mask], alpha = 0.5, color = 'b', label = 'HiTS')
    ax.legend(loc = 2)
    ax.set_xlabel("Dipole length [pix]")
    ax.set_ylabel("Residual length [pix]")
    plt.savefig("%s/dipoleresiduallengths_%s_%s_%02i-%02i.png" % (outdir, field, CCD, epoch2, epoch1), bbox_inches = 'tight')

    # dipole angle vs residual angle
    del fig, ax
    fig, ax = plt.subplots()
    mask = np.array((z1 > fluxlim))
    ax.scatter(dipoleanglehits[match & mask], residualanglehits[match & mask], alpha = 0.5, color = 'b', label = 'HiTS')
    ax.legend(loc = 2)
    ax.set_xlabel("Dipole angle [deg]")
    ax.set_ylabel("Residual angle [deg]")
    plt.savefig("%s/dipoleresidualangles_%s_%s_%02i-%02i.png" % (outdir, field, CCD, epoch2, epoch1), bbox_inches = 'tight')

    # histogram of raw dipole values
    fluxlim = int(sys.argv[1])
    del fig, ax
    fig, ax = plt.subplots()
    mask = np.array(z1 < fluxlim)
    stampshits = stampshits[mask].flatten()
    stampslsst = stampslsst[mask].flatten()
    histlim = 500
    bins = np.linspace(-histlim, histlim, 150)
    ax.hist(stampshits[np.isfinite(stampshits) & (stampshits < histlim) & (stampshits > -histlim)], color = 'b', label = 'HiTS', bins = bins, alpha = 0.5, log = True, linewidth = 0)
    ax.hist(stampslsst[np.isfinite(stampslsst) & (stampslsst < histlim) & (stampslsst > -histlim)], color = 'r', label = 'LSST', bins = bins, alpha = 0.5, log = True, linewidth = 0)
    ax.legend(loc = 1)
    ax.set_title("Star flux < %i ADU" % fluxlim)
    ax.set_xlim(-histlim, histlim)
    ax.set_xlabel("dipole signal [ADU]")
    ax.set_ylim(1., 1e5)
    plt.savefig("%s/dipole_%s_%s_%02i-%02i_%i.png" % (outdir, field, CCD, epoch2, epoch1, fluxlim), bbox_inches = 'tight')

    # histogram of dipole SNR values
    del fig, ax
    fig, ax = plt.subplots()
    mask = np.array(z1 < fluxlim)
    stampshitsSNR = stampshitsSNR[mask].flatten()
    stampslsstSNR = stampslsstSNR[mask].flatten()
    histlim = 50
    bins = np.linspace(-histlim, histlim, 150)
    ax.hist(stampshitsSNR[np.isfinite(stampshits) & (stampshits < histlim) & (stampshits > -histlim)], color = 'b', label = 'HiTS', bins = bins, alpha = 0.5, log = True, linewidth = 0, normed = True)
    ax.hist(stampslsstSNR[np.isfinite(stampslsst) & (stampslsst < histlim) & (stampslsst > -histlim)], color = 'r', label = 'LSST', bins = bins, alpha = 0.5, log = True, linewidth = 0, normed = True)
    ax.plot(bins, stats.norm.pdf(bins), c = 'k', label = "G(0, 1)")
    ax.legend(loc = 1)
    ax.set_title("Star flux < %i ADU" % fluxlim)
    ax.set_xlim(-histlim, histlim)
    ax.set_xlabel("dipole SNR")
    ax.set_ylim(1e-8, 1)
    ax.set_xlim(-10, 10)
    plt.savefig("%s/dipoleSNR_%s_%s_%02i-%02i_%i.png" % (outdir, field, CCD, epoch2, epoch1, fluxlim), bbox_inches = 'tight')

