#!/usr/bin/env python

"""
rmselfcal.py

Perform RM-selfcal as developed by M Brentjens.
Needs Q and U images for single time and frequency steps.

Copyright 2017 by Andreas Horneffer

pyrmsynth is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
"""

import time, calendar, os
from astropy.io import fits
import astropy.time as astrotime
import numpy as np
import rm_tools
from scipy.optimize import curve_fit

def progress(width, percent):
    """
    Print a progress-bar
    """
    import sys, math
    marks = math.floor(width * (percent / 100.0))
    spaces = math.floor(width - marks)
    loader = '[' + ('=' * int(marks)) + (' ' * int(spaces)) + ']'
    sys.stdout.write("%s %d%%\r" % (loader, percent))
    if percent >= 100:
        sys.stdout.write("\n")
    sys.stdout.flush()

def create_memmap_file_and_array(fn, SHAPE, DTYPE):
    """
    Creates an empty numpy memmap array and the associated file on disk.
    """
    if np.isscalar(SHAPE):
        npix = SHAPE
    else:
        npix = 1
        for i in range(len(SHAPE)):
            npix = npix * SHAPE[i]
    with open(fn, 'wb') as f:
        # OPEN THE FILE, SKIP TO THE END, AND WRITE A 0
        # Quickly creates an empty file on disk.
        f.seek(npix * DTYPE.itemsize - 1)
        f.write('\x00')
    m = np.memmap(fn, dtype=DTYPE, shape=SHAPE)
    return m

def sort_check_input_files(infiles, maskfile=None, debug=False):
    """
    Sort the input files by time and frequency

    Parameters:
    -----------
    infiles : list or array of str
        Names (pathes) of the input maps.
    maskfile : str (optional)
        Name (path) of a mask file that will be chcked for conformity 
        with the polarization files
    debug : bool (optional)
        Print out some more debug information
    """
    filesdict = {}

    print "Sorting all files by time and frequency"
    #read in first fits-file
    inheader = fits.getheader(infiles[0])    
    #datestring = inheader['DATE-OBS']
    #datestamp = calendar.timegm(time.strptime(datestring.split('.')[0]+' UTC',"%Y-%m-%dT%H:%M:%S %Z"))
    # LOFAR (or CASA MS?) style datestamp is MJD in seconds.
    t = astrotime.Time(inheader['DATE-OBS'],format='isot')
    datestamp = t.mjd*24.*3600.
    try:
        assert inheader['NAXIS']   == 4
        assert inheader['NAXIS3']  == 1
        assert inheader['NAXIS4']  == 1
        assert inheader['CTYPE3']  == 'FREQ'
        assert inheader['CUNIT3']  == 'Hz'
        assert inheader['CTYPE4']  == 'STOKES'
    except AssertionError:
        raise ValueError("Input file %s does not match required structure"%(infiles[0]))
    NAXIS1 = inheader['NAXIS1']
    NAXIS2 = inheader['NAXIS2']
    CRVAL1 = inheader['CRVAL1']
    CRVAL2 = inheader['CRVAL2']
    CDELT1 = inheader['CDELT1']
    CDELT2 = inheader['CDELT2']
    CDELT3 = inheader['CDELT3']
    #filesdict['dnu'] = CDELT3
    #filesdict['RA-len'] = NAXIS1
    #filesdict['DEC-len'] = NAXIS2
    filesdict[datestamp] = {}
    filesdict[datestamp]['frequencies'] = [inheader['CRVAL3']]
    if inheader['CRVAL4'] == 2.:
        filesdict[datestamp]['q-files'] = [infiles[0]]
        filesdict[datestamp]['u-files'] = [None]
    elif inheader['CRVAL4'] == 3.:
        filesdict[datestamp]['q-files'] = [None]
        filesdict[datestamp]['u-files'] = [infiles[0]]
    else:
        raise ValueError("Input file %s is not a Q- or U-map."%(infiles[0]))

    # Check the mask-file if needed
    # (If this failes, we don't need to wait until all the polarization maps are processed.)
    if maskfile:
        inheader = fits.getheader(maskfile) 
        try:
            assert inheader['NAXIS']   == 4
            assert inheader['NAXIS3']  == 1
            assert inheader['NAXIS4']  == 1
            assert NAXIS1 == inheader['NAXIS1']
            assert NAXIS2 == inheader['NAXIS2']
            assert CRVAL1 == inheader['CRVAL1']
            assert CRVAL2 == inheader['CRVAL2']
            assert CDELT1 == inheader['CDELT1']
            assert CDELT2 == inheader['CDELT2']
        except AssertionError:
            print "Mask file %s has different structure than %s"%(maskfile,infiles[0])
            raise 
        
    
    # read in all the other files and compare to first file
    for (num, filename) in enumerate(infiles[1:]):
        inheader = fits.getheader(filename) 
        try:
            assert inheader['NAXIS']   == 4
            assert inheader['NAXIS3']  == 1
            assert inheader['NAXIS4']  == 1
            assert inheader['CTYPE3'].upper() == 'FREQ'
            assert inheader['CUNIT3'].upper() == 'HZ'
            assert inheader['CTYPE4'].upper() == 'STOKES'
            assert NAXIS1 == inheader['NAXIS1']
            assert NAXIS2 == inheader['NAXIS2']
            assert CRVAL1 == inheader['CRVAL1']
            assert CRVAL2 == inheader['CRVAL2']
            assert CDELT1 == inheader['CDELT1']
            assert CDELT2 == inheader['CDELT2']
            assert CDELT3 == inheader['CDELT3']
        except AssertionError:
            print "Input file %s has different structure than %s"%(filename,infiles[0])
            raise
        #datestring = inheader['DATE-OBS']
        #datestamp = calendar.timegm(time.strptime(datestring.split('.')[0]+' UTC',"%Y-%m-%dT%H:%M:%S %Z"))
        t = astrotime.Time(inheader['DATE-OBS'],format='isot')
        datestamp = t.mjd*24.*3600.
        freq = inheader['CRVAL3']
        if datestamp not in filesdict:
            filesdict[datestamp] = {}
            filesdict[datestamp]['frequencies'] = []
            filesdict[datestamp]['q-files'] = []
            filesdict[datestamp]['u-files'] = []
        if freq not in filesdict[datestamp]['frequencies']:
            filesdict[datestamp]['frequencies'].append(freq)
            if inheader['CRVAL4'] == 2.:
                filesdict[datestamp]['q-files'].append(filename)
                filesdict[datestamp]['u-files'].append(None)
            elif inheader['CRVAL4'] == 3.:
                filesdict[datestamp]['q-files'].append(None)
                filesdict[datestamp]['u-files'].append(filename)
            else:
                raise ValueError("Input file %s is not a Q- or U-map."%(filename))
        else:
            findex = filesdict[datestamp]['frequencies'].index(freq)
            if inheader['CRVAL4'] == 2.:
                filesdict[datestamp]['q-files'][findex] = filename
            elif inheader['CRVAL4'] == 3.:
                filesdict[datestamp]['u-files'][findex] = filename
            else:
                raise ValueError("Input file %s is not a Q- or U-map."%(filename))
        if (num % 100) == 0:
            progress(70,100.*num/len(infiles))
    progress(70,100.)
    # clean up the cube (remove entries with only Q or U maps)
    for datestamp in filesdict.keys():
        findex = 0
        while findex < len(filesdict[datestamp]['frequencies']):
            if ( (not filesdict[datestamp]['q-files'][findex]) or
                 (not filesdict[datestamp]['u-files'][findex]) ):
                if debug:
                    print "removing entry:",datestamp, filesdict[datestamp]['frequencies'][findex], filesdict[datestamp]['q-files'][findex], filesdict[datestamp]['u-files'][findex]
                filesdict[datestamp]['frequencies'].pop(findex)
                filesdict[datestamp]['q-files'].pop(findex)
                filesdict[datestamp]['u-files'].pop(findex)
            else:
                findex += 1

    return (NAXIS1, NAXIS2, CDELT3, filesdict)


def get_datacube_from_files(q_files, u_files, RAlen, DEClen, temfilename):
    """
    Create the (dirty) RM-cube for a given set of files

    Parameters:
    -----------
    q_files : list or array of str
        Names (pathes) of the input Q maps.
    u_files : list or array of str
        Names (pathes) of the input U maps.
    RAlen : int
        size of the input maps in RA
    DEClen : int
        size of the input maps in DEC
    temfilename : str
        Name (path) of the temporary file in which to store the input data.

    Notes:
    ------
    - q_files, u_files,  have to have the same length
    - the size of the output cube will be:
      len(x_files) x  DEClen x RAlen
    """
    nchan = len(q_files) 
    assert len(u_files) == nchan

    incube = create_memmap_file_and_array(temfilename,
                                          (nchan, DEClen, RAlen),
                                          np.dtype('complex128'))
    for idx in xrange(nchan):
        if q_files[idx] != None and u_files[idx] != None :
            try:
                tdata_q = fits.getdata(q_files[idx])
                tdata_u = fits.getdata(u_files[idx])
                incube.real[idx] = tdata_q[0, 0, : , : ]
                incube.imag[idx] = tdata_u[0, 0, : , : ]
            except IOError:
                print "IO-Exception! idx: %d, Q-file: %s, U-file: %s"%(idx,q_files[idx],u_files[idx])
                raise
    return incube

def correlate_cubes(reference_cube, reference_frequencies, data_cube, data_frequencies,
                    dnu, phi_values, mask=None, numrows=None, storage_file=None, debug=False):
    """
    correlate the two data-cubes

    Parameters:
    -----------
    reference_cube : 3-d numpy array of complex
        Array with the image data from the reference time-step 
        with dimensions (nchan, DEClen, RAlen)
    reference-frequencies : 1-d numpy array of float
        Array with the frequencies of the reference data, 
        length must be equal to nchan from reference cube
    data_cube : 3-d numpy array of complex
        Array with the image data from the given time-step 
        with dimensions (nchan, DEClen, RAlen)
        DEClen and RAlen must be equal to reference_cube
    data-frequencies : 1-d numpy array of float
        Array with the of the data, length must be equal to nchan from data cube    
    dnu : float
        Bandwidth of a data in a single frequency channel
    phi_values : 1-d array or list of float
        (Relative-)FR values for which the correlation should be computed
    mask : 2-d numpy array
        Mask array with dimensions (DEClen, RAlen). DEClen, RAlen have to be the 
        same as in input cubes. Only lines-of-sight in which mask is _not_ set will
        be used for the correlation. (I.e. mask can be a stokes-I clean mask.)
    numrows : int
        Process only this many rows in DEC.
    storage_file : str
        Name (path) of a file in which to store the corrlation data
    debug : bool
        Print debug output
        
    """
    # check data consitency
    assert reference_cube.ndim == 3
    assert data_cube.ndim == 3
    assert reference_cube.shape[0] == len(reference_frequencies)
    assert data_cube.shape[0] == len(data_frequencies)
    assert reference_cube.shape[1] == data_cube.shape[1] 
    assert reference_cube.shape[2] == data_cube.shape[2]
    if type(mask) is np.ndarray:
        assert mask.shape[0] == reference_cube.shape[1]
        assert mask.shape[1] == reference_cube.shape[2]
        not_use_mask = False
    else:
        not_use_mask = True
    # make list in indices of common frequencies
    ref_indices = []
    data_indices = []
    ref_idx = 0
    data_idx = 0
    while ref_idx < len(reference_frequencies) and data_idx < len(data_frequencies):
        if reference_frequencies[ref_idx] == data_frequencies[data_idx]:
            ref_indices.append(ref_idx)
            data_indices.append(data_idx)
            ref_idx += 1
            data_idx += 1
        elif reference_frequencies[ref_idx] < data_frequencies[data_idx]:
            ref_idx += 1
        else:
            data_idx += 1
    if debug:
        print("correlate_cubes: Num-Frequencies: reference: %d, data: %d, matched: %d"%(len(reference_frequencies), len(data_frequencies), len(ref_indices)))

    # make RM-synthesis object    
    weights = np.ones(len(ref_indices))
    frequency_array = np.array(reference_frequencies)
    rms = rm_tools.RMSynth(frequency_array[ref_indices], dnu, phi_values, weights)
    
    if storage_file:
        corr_cube = np.zeros((len(phi_values),reference_cube.shape[1],reference_cube.shape[2]), dtype="float32")

    sum_corr = np.zeros(len(phi_values),dtype="float32")
    if numrows and int(numrows)>0:
        numrows=min(int(numrows),reference_cube.shape[1])
        if debug:
            print "Processing",numrows,"rows."
    else:
        numrows = reference_cube.shape[1]      
    progress(40,0.)
    for DEC_idx in xrange(numrows):
        for RA_idx in xrange(reference_cube.shape[2]):
            if not_use_mask or not mask[DEC_idx,RA_idx]:
                los_data = data_cube[data_indices,DEC_idx,RA_idx] * \
                           reference_cube[ref_indices,DEC_idx,RA_idx].conj()
                los_corr = np.abs(rms.compute_dirty_image(los_data))
                if storage_file:
                    corr_cube[:,DEC_idx,RA_idx] = los_corr
                sum_corr += los_corr
        progress(40,100.*(DEC_idx+1.001)/reference_cube.shape[1])
    if storage_file:
        np.save(storage_file, corr_cube)
    print ""
    
    return sum_corr


def gauss_function(x, a, x0, sigma):
    """
    Compute a simple Gauss-distribution.
    """
    return a*np.exp(-(x-x0)**2/(2*sigma**2))
          
def do_RMselfcal(infiles,maskfile=False,invert_mask=False,numrows=None,
                 debug=False,debug_plots=False):
    """
    Perform RM-selfcal

    Parameters:
    -----------
    infiles : list of str
        Names (pathes) of the input Q- and U-maps.
    maskfile : str (optional)
        Name (path) of a maskfile. Only positions not in the mask 
        will be used for the correlation. (I.e. a clean-mask for Stokes-I
        cleaning will reduce the influence of instrumental polarization.)
    invert_mask : bool (optional)
        Invert the logic of the mask, i.e. use th position where the mask is set.
    numrows : int (optional)
        Number of linese in DEC that should be processed. Intended to reduce 
        processing while debugging the code.
    debug : bool (optional)
        Print out debug information.
    debug_plots : bool (optional)
        Create plots (png files) with the summed correlation for each time-step.

    Results:
    --------
    timestamps : list of float
        List with the time-stamps (in Unix time) of the time-steps
    dFR_values : 1-d array of float with length of timestamps
        Faraday rotation difference to reference time of the time-steps
    sigma_values : 1-d array of float with length of timestamps
        Width of the gaussian fit for the time-steps
    """
    # setup
    reffilename = 'rmselfcal_refcube.dat'
    datafilename = 'rmselfcal_tmpcube.dat'
    phi_values = np.arange(-6.,6.,0.1)
    if invert_mask and debug:
        print "Going to invert the mask logic!"
    
    # sort the input files
    (RAlen, DEClen, dnu, filesdict) = sort_check_input_files(infiles, maskfile=maskfile, debug=True)
    
    # get the reference cube
    print "Loading the reference data-cube."
    timestamps = filesdict.keys()
    timestamps.sort()
    refidx = len(timestamps)/2
    reftime = timestamps[refidx]
    refcube = get_datacube_from_files(filesdict[reftime]['q-files'], filesdict[reftime]['u-files'],
                                      RAlen, DEClen, reffilename)

    if maskfile:
        maskdata = fits.getdata(maskfile)
        if maskdata.ndim == 4:
            mask = maskdata[0,0,:,:].astype(np.bool)
        elif maskdata.ndim == 3:
            mask = maskdata[0,:,:].astype(np.bool)
        elif maskdata.ndim == 2:
            mask = maskdata[:,:].astype(np.bool)
        else:
            raise ValueError("Mask %s has unsupported number of dimensions %d"%(maskfile, maskdata.ndim))
        if invert_mask:
            mask = np.logical_not(mask)
    else:
        mask = None
            
    # correlate time-steps with reference cube
    print "Performing the correlarions."
    dFR_values = np.zeros(len(timestamps))
    sigma_values = np.zeros(len(timestamps))    
    peak_values = np.zeros(len(timestamps))    
    for (dateidx, datestamp) in enumerate(timestamps):
        print "Starting timeindex %d of %d"%(dateidx+1,len(timestamps))
        if datestamp == reftime:
            if debug: print "This is the reference time, continuing."
            dFR_values[dateidx] = 0.
            sigma_values[dateidx] = 0.
            continue
        if debug: print "Loading the comparison data-cube."
        datacube = get_datacube_from_files(filesdict[datestamp]['q-files'],
                                           filesdict[datestamp]['u-files'],
                                           RAlen, DEClen, datafilename)
        if debug: print "Running the correlations."
        FR_corr = correlate_cubes(refcube, filesdict[reftime]['frequencies'],
                                  datacube, filesdict[datestamp]['frequencies'],
                                  dnu, phi_values, 
                                  mask=mask, numrows=numrows,
                                  storage_file="corr-cube.npy",
                                  debug=debug)
        # fit a Gaussian to get the position of the maximum
        maxidx = np.argmax(FR_corr)
        startvals = [FR_corr[maxidx],(phi_values[maxidx]+.01),1.]
        fstart = max(0, maxidx-5)
        fend = min(len(FR_corr), maxidx+5)
        (popt, pcov) = curve_fit(gauss_function, phi_values[fstart:fend], FR_corr[fstart:fend],
                                 startvals)
        if debug:
            print "startvals",startvals,"popt",popt,"pcov",pcov
        if debug_plots:
            import matplotlib as mpl
            mpl.use('Agg')
            import matplotlib.pyplot as plt
            plt.plot(phi_values, FR_corr)
            plt.plot(phi_values[fstart:fend], gauss_function(phi_values[fstart:fend], popt[0], popt[1],  popt[2]))
            figurename = "corr-plot-%+05.0f.png"%((datestamp-reftime)/60.)
            plt.savefig(figurename,dpi=100)
            plt.clf()
        peak_values[dateidx] = popt[0]
        dFR_values[dateidx] = popt[1]
        sigma_values[dateidx] = popt[2]

    # return list of time-stamps and FR differences
    return timestamps, dFR_values, sigma_values, peak_values


if __name__ == '__main__':
    import argparse
    import glob

    descriptiontext = "Calculates differential Faraday-Rotation values from polarization images.\n"
    parser = argparse.ArgumentParser(description=descriptiontext)
    parser.add_argument('im_file_pattern', help='Glob-able filename-pattern of input images. '
                        '(Usually needs to be put in quotation marks: \" or \')')
    parser.add_argument("outfile", help='Name of (gnuplot-style) output file to be written.')
    parser.add_argument("-m", "--mask", dest="Mask",
                        help="Name (path) of a mask file. Only positions not in the mask "
                        "will be used for the correlation. (I.e. a clean-mask for Stokes-I "
                        "cleaning will reduce the influence of instrumental polarization.)")
    parser.add_argument("--invert_mask", action="store_true",
                        help="Invert the logic of the mask, i.e. use the position in which "
                        "the mask is set.")
    parser.add_argument("-n", "--num-rows", dest="Rows",
                        help="Process only this many rows in DEC when calucating the correlations. "
                        "Mainly for debugging purposes.")
    parser.add_argument("-p", "--plots", dest="Plots", action="store_true",
                        help="Write out png plots with the correlation for each time-step.")
    parser.add_argument("-d", "--debug", dest="Debug", action="store_true",
                        help="Generate Debug output.")
    

    
    args = parser.parse_args()

    imlist = glob.glob(args.im_file_pattern)

    (timestamps, dFR_values, sigma_values, peak_values) = do_RMselfcal(imlist, maskfile=args.Mask,
                                                                       numrows=args.Rows,
                                                                       invert_mask=args.invert_mask,
                                                                       debug=args.Debug, debug_plots=args.Plots)
    fd = open(args.outfile,"w")
    fd.write("# timestamp dFR_value sigma_value\n")
    for (timestamp, dFR, sigma, peak) in zip(timestamps, dFR_values, sigma_values, peak_values):
        fd.write("%f %f %f %f\n"%(timestamp, dFR, sigma, peak))
    fd.close()
