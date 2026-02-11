"""
calc_flow contains the functions necessary to calculate optical flow in 2D and
3D, as well as a function for parsing files and their metadata.
"""

import math
import numpy as np
from scipy.ndimage import correlate1d
from pathlib import Path
import os
import re
import tifffile as tf
import pandas as pd
from datetime import datetime
import sys
from natsort import natsorted

def calc_flow2D(images,xySig=3,tSig=1,wSig=4):
    """
    Calculate two-dimensional optical flow fields from input images.

    The calc_flow2D function calculates optical flow velocities for a single
    z-slice. Surrounding images in time are necessary to perform the 
    calculations. To peform calculations on an entire timelapse, see the
    function parse_flow.

    This script uses the convention that (0,0) is located in the upper-left
    corner of an image. This is inline with conventions used in other programs
    (e.g., ImageJ/FIJI), but note that it means that positive y-velocities point
    down, which can be non-intuitive in some cases.

    ARGS:
    images: 3D numpy array with dimensions N_T, N_Y, N_X
             N_T should be odd as only the central timepoint will be analyzed.
             N_T must be greater than or equal to 3*tSig+1.
    xySig:  sigma value for smoothing in all spatial dimensions. Default 3.
             Larger values remove noise but remove spatial detail.
    tSig:   sigma value for smoothing in the temporal dimension. Default 1.
             Larger values remove noise but remove temporal detail.
    wSig:   sigma value for Lucas-Kanade neighborhood. Default is 4.
             Larger values include a larger neighboorhood in the
             Lucas-Kanade constraint and will smooth over small features.

    RETURNS:
    vx:    Velocity in the x direction, reported as pixels/frame
    vy:    Velocity in the y direction, reported as pixels/frame
    rel:   Reliability, the smallest eigenvalue of (A'wA)
            This is a measure of confidence in the linear algebra solution
            and can be used to mask the velocities for downstream analysis.
    """

    ### Check the function inputs ##############################################
    # Check that the images are 2D + time
    if not(len(images.shape)==3): 
        sys.exit('ERROR: Input image must be a 3D matrix with dimensions N_T, N_Y, N_X')
    # Check image size against tSig
    Nt = images.shape[0] 
    if Nt < 6*tSig+1:
        # The kernel size is 3*tSig in time (or 6*tSig total)
        # There is also a central pixel, so need at least 6*tSig+1 images in t
        sys.exit('ERROR: Input images will lead to edge effects. N_T must be >= 6*tSig+1')
    # Check for an odd number of frames
    if not(Nt % 2):        
        sys.exit('ERROR: Input images must have an odd number of timepoints. Only the central time point is analyzed')
    NtSlice = math.ceil(Nt/2)-1 # -1 because python indexing starts from 0
    # Use float for calculations
    images = images.astype(np.float64)


    ### Set up filters #########################################################
    # Common terms
    x = np.arange(-math.ceil(3*xySig),math.ceil(3*xySig)+1)
    xySig2 = xySig/4
    y = np.arange(-math.ceil(3*xySig2),math.ceil(3*xySig2)+1)
    fderiv = np.exp(-x*x/2/xySig/xySig)/math.sqrt(2*math.pi)/xySig
    fsmooth = np.exp(-y*y/2/xySig2/xySig2)/math.sqrt(2*math.pi)/xySig2
    gderiv = x/xySig/xySig
    gsmooth = 1

    # Build y-gradient filter kernels (along first spatial dimension)
    yFil1 = (fderiv*gderiv)
    xFil1 = (fsmooth*gsmooth)
    # Build x-gradient filter kernels (along second spatial dimension)
    yFil2 = (fsmooth*gsmooth)
    xFil2 = (fderiv*gderiv)

    # Build t-gradient filter kernels (t = third dimension)
    t = np.arange(-math.ceil(3*tSig),math.ceil(3*tSig)+1)
    fx = np.exp(-x*x/2/xySig/xySig)/math.sqrt(2*math.pi)/xySig
    ft = np.exp(-t*t/2/tSig/tSig)/math.sqrt(2*math.pi)/tSig
    gx = 1
    gt = t/tSig/tSig
    yFil3 = (fx*gx)
    xFil3 = yFil3
    tFil3 = ft*gt

    # Structure tensor -- Lucas Kanade neighborhood filter
    wRange = np.arange(-math.ceil(3*wSig),math.ceil(3*wSig)+1)
    gw = np.exp(-wRange*wRange/2/wSig/wSig)/math.sqrt(2*math.pi)/wSig
    yFil4 = gw
    xFil4 = gw

    # Throughout will use del to keep the memory clear as this processing is memory intensive
    del gderiv, gsmooth, gt, gw, gx, ft, fx, fsmooth, fderiv, x, y, t, wRange


    ### Spatial and Temporal Gradients #########################################
    # Spatial gradients require only the frame of interest, while  the temporal
    # gradient requires N_T >= 2*3*tSig+1. Keep only the relevant slice of dtI
    # after it is calculated.

    # dtI is split into two steps to save memory and processing time
    dtI = correlate1d(images, tFil3, axis=0, mode='nearest')
    dtI = dtI[NtSlice,:]
    images = images[NtSlice,:]
    dtI = correlate1d(correlate1d(dtI, yFil3, axis=0, mode='nearest'), xFil3, axis=1, mode='nearest')
    del xFil3, yFil3, tFil3
    
    dyI = correlate1d(correlate1d(images, yFil1, axis=0, mode='nearest'), xFil1, axis=1, mode='nearest')
    del xFil1, yFil1

    dxI = correlate1d(correlate1d(images, yFil2, axis=0, mode='nearest'), xFil2, axis=1, mode='nearest')
    del xFil2, yFil2
    del images


    ### Structure Tensor Inputs ################################################
    # The following calculations are for the individual elements of the
    # matrices required for the optical flow calculation, incorporating
    # Gaussian weighting into the Lucas-Kanade constraint.

    # Time components
    wdtx = correlate1d(correlate1d(dxI*dtI, yFil4, axis=0, mode='nearest'), xFil4, axis=1, mode='nearest')
    wdty = correlate1d(correlate1d(dyI*dtI, yFil4, axis=0, mode='nearest'), xFil4, axis=1, mode='nearest')
    del dtI

    # Spatial Components
    wdxy = correlate1d(correlate1d(dxI*dyI, yFil4, axis=0, mode='nearest'), xFil4, axis=1, mode='nearest')
    wdx2 = correlate1d(correlate1d(dxI*dxI, yFil4, axis=0, mode='nearest'), xFil4, axis=1, mode='nearest')
    del dxI
    wdy2 = correlate1d(correlate1d(dyI*dyI, yFil4, axis=0, mode='nearest'), xFil4, axis=1, mode='nearest')
    del dyI
    del xFil4, yFil4


    ### Optical Flow Calculations ##############################################
    # Equation is v = (A' w A)^-1 A' w b
    # A = -[dxI dyI]
    # b = [dtI]
    # w multiplication is incorporated in the structure tensor inputs above
    # A' w b = -[wdtx wdty]  (minus sign because of negative sign on A)
    # (A' w A) = [a=wdx2 b=wdxy ; c=wdxy d=wdy2]
    # A^-1 = [a b ; c d]^-1 = (1/det(A))[d - b; -c a]
    determinant = (wdx2*wdy2) - (wdxy*wdxy)
    vx = ((determinant+np.finfo(float).eps)**-1)*((wdy2*-wdtx)+(-wdxy*-wdty))
    vy = ((determinant+np.finfo(float).eps)**-1)*((-wdxy*-wdtx)+(wdx2*-wdty))
    del wdtx, wdty, wdxy


    ### Eigenvalues for Reliability ############################################
    # solve det(A^T w A - lamda I) = 0
    # (A' w A) = [a=wdx2 b=wdxy ; c=wdxy d=wdy2]
    trace = wdx2 + wdy2
    del wdx2, wdy2

    L1 = (trace + np.sqrt(trace**2 - 4*determinant))/2
    L2 = (trace - np.sqrt(trace**2 - 4*determinant))/2
    rel = np.real(np.minimum(L1,L2))
    del L1, L2


    ### Return Outputs #########################################################
    return vx, vy, rel

def calc_flow3D(images,xyzSig=3,tSig=1,wSig=4):
    """
    Calculate three-dimensional optical flow fields from input z-stacks.

    The calc_flow3D function calculates optical flow velocities for a single
    z-stack of images. Surrounding z-stacks in time are necessary to perform
    the calculations. To peform calculations on an entire timelapse, see
    the function parse_flow.

    This script uses the convention that (0,0) is located in the upper-left
    corner of an image. This is inline with conventions used in other
    programs (e.g., ImageJ/FIJI), but note that it means that positive
    y-velocities point down, which can be non-intuitive in some cases.

    ARGS:
    images: 4D array with dimensions N_T, N_Z, N_Y, N_X
             N_T should be odd as only the central timepoint will be analyzed.
             N_T must be greater than or equal to 6*tSig+1.
    xyzSig:  sigma value for smoothing in all spatial dimensions. Default 3.
             Larger values remove noise but remove spatial detail.
    tSig:   sigma value for smoothing in the temporal dimension. Default 1.
             Larger values remove noise but remove temporal detail.
    wSig:   sigma value for Lucas-Kanade neighborhood. Default is 4.
             Larger values include a larger neighboorhood in the
             Lucas-Kanade constraint and will smooth over small features.

    RETURNS:
    vx:    Velocity in the x direction, reported as pixels/frame
    vy:    Velocity in the y direction, reported as pixels/frame
    vz:    Velocity in the z direction, reported as pixels/frame
    rel:   Reliability, the smallest eigenvalue of (A'wA)
            This is a measure of confidence in the linear algebra solution
            and can be used to mask the velocities for downstream analysis.
    """

    ### Check the function inputs ##############################################
    # Check that the images are 2D + time
    if not(len(images.shape)==4): 
        sys.exit('ERROR: Input image must be a 3D matrix with dimensions N_T, N_Z, N_Y, N_X')
    # Check image size against tSig
    Nt = images.shape[0] 
    if Nt < 6*tSig+1:
        # The kernel size is 3*tSig in time (or 6*tSig total)
        # There is also a central pixel, so need at least 6*tSig+1 images in t
        sys.exit('ERROR: Input images will lead to edge effects. N_T must be >= 6*tSig+1')
    # Check for an odd number of frames
    if not(Nt % 2):
        sys.exit('ERROR: Input images must have an odd number of timepoints. Only the central time point is analyzed')
    NtSlice = math.ceil(Nt/2)-1 # -1 because python indexing starts from 0
    # Use float for calculations
    images = images.astype(np.float64)


    ### Set up filters #########################################################
    # Common terms
    x = np.arange(-math.ceil(3*xyzSig),math.ceil(3*xyzSig)+1)
    xyzSig2 = xyzSig/4
    y = np.arange(-math.ceil(3*xyzSig2),math.ceil(3*xyzSig2)+1)
    fderiv = np.exp(-x*x/2/xyzSig/xyzSig)/math.sqrt(2*math.pi)/xyzSig
    fsmooth = np.exp(-y*y/2/xyzSig2/xyzSig2)/math.sqrt(2*math.pi)/xyzSig2
    gderiv = x/xyzSig/xyzSig
    gsmooth = 1

    # Build y-gradient filter kernels
    yFil1 = (fderiv*gderiv)
    xFil1 = (fsmooth*gsmooth)
    zFil1 = (fsmooth*gsmooth)
    # Build x-gradient filter kernels
    yFil2 = (fsmooth*gsmooth)
    xFil2 = (fderiv*gderiv)
    zFil2 = (fsmooth*gsmooth)
    # Build z-gradient filter kernels
    yFil3 = (fsmooth*gsmooth)
    xFil3 = (fsmooth*gsmooth)
    zFil3 = (fderiv*gderiv)

    # Build t-gradient filter kernels (t = third dimension)
    t = np.arange(-math.ceil(3*tSig),math.ceil(3*tSig)+1)
    fx = np.exp(-x*x/2/xyzSig/xyzSig)/math.sqrt(2*math.pi)/xyzSig
    ft = np.exp(-t*t/2/tSig/tSig)/math.sqrt(2*math.pi)/tSig
    gx = 1
    gt = t/tSig/tSig
    yFil4 = (fx*gx)
    xFil4 = (fx*gx)
    zFil4 = (fx*gx)
    tFil4 = ft*gt

    # Structure tensor -- Lucas Kanade neighborhood filter
    wRange = np.arange(-math.ceil(3*wSig),math.ceil(3*wSig)+1)
    gw = np.exp(-wRange*wRange/2/wSig/wSig)/math.sqrt(2*math.pi)/wSig
    yFil5 = gw
    xFil5 = gw
    zFil5 = gw

    # Throughout will use del to keep the memory clear as this processing is memory intensive
    del gderiv, gsmooth, gt, gw, gx, ft, fx, fsmooth, fderiv, x, y, t, wRange


    ### Spatial and Temporal Gradients #########################################
    # Spatial gradients require at least N_T = 3 to avoid edge effets, while 
    # the temporal gradient requires N_T >= 6*tSig+1.
    dtI = correlate1d(images, tFil4, axis=0, mode='nearest')
    dtI = dtI[NtSlice,:]
    images = images[NtSlice,:]
    dtI = correlate1d(correlate1d(correlate1d(dtI, yFil4, axis=1, mode='nearest'), xFil4, axis=2, mode='nearest'), zFil4, axis=0, mode='nearest')
    del xFil4, yFil4, zFil4, tFil4

    dyI = correlate1d(correlate1d(correlate1d(images, yFil1, axis=1, mode='nearest'), xFil1, axis=2, mode='nearest'), zFil1, axis=0, mode='nearest')
    del xFil1, yFil1, zFil1

    dxI = correlate1d(correlate1d(correlate1d(images, yFil2, axis=1, mode='nearest'), xFil2, axis=2, mode='nearest'), zFil2, axis=0, mode='nearest')
    del xFil2, yFil2, zFil2

    dzI = correlate1d(correlate1d(correlate1d(images, yFil3, axis=1, mode='nearest'), xFil3, axis=2, mode='nearest'), zFil3, axis=0, mode='nearest')
    del xFil3, yFil3, zFil3
    
    del images


    ### Structure Tensor Inputs ################################################
    # The following calculations are for the individual elements of the
    # matrices required for the optical flow calculation, incorporating
    # Gaussian weighting into the Lucas-Kanade constraint.

    # Time components
    wdtx = correlate1d(correlate1d(correlate1d(dxI*dtI, yFil5, axis=1, mode='nearest'), xFil5, axis=2, mode='nearest'), zFil5, axis=0, mode='nearest')
    wdty = correlate1d(correlate1d(correlate1d(dyI*dtI, yFil5, axis=1, mode='nearest'), xFil5, axis=2, mode='nearest'), zFil5, axis=0, mode='nearest')
    wdtz = correlate1d(correlate1d(correlate1d(dzI*dtI, yFil5, axis=1, mode='nearest'), xFil5, axis=2, mode='nearest'), zFil5, axis=0, mode='nearest')
    del dtI

    # Spatial Components
    wdxy = correlate1d(correlate1d(correlate1d(dxI*dyI, yFil5, axis=1, mode='nearest'), xFil5, axis=2, mode='nearest'), zFil5, axis=0, mode='nearest')
    wdxz = correlate1d(correlate1d(correlate1d(dxI*dzI, yFil5, axis=1, mode='nearest'), xFil5, axis=2, mode='nearest'), zFil5, axis=0, mode='nearest')
    wdx2 = correlate1d(correlate1d(correlate1d(dxI*dxI, yFil5, axis=1, mode='nearest'), xFil5, axis=2, mode='nearest'), zFil5, axis=0, mode='nearest')
    del dxI
    wdyz = correlate1d(correlate1d(correlate1d(dyI*dzI, yFil5, axis=1, mode='nearest'), xFil5, axis=2, mode='nearest'), zFil5, axis=0, mode='nearest')
    wdy2 = correlate1d(correlate1d(correlate1d(dyI*dyI, yFil5, axis=1, mode='nearest'), xFil5, axis=2, mode='nearest'), zFil5, axis=0, mode='nearest')
    del dyI
    wdz2 = correlate1d(correlate1d(correlate1d(dzI*dzI, yFil5, axis=1, mode='nearest'), xFil5, axis=2, mode='nearest'), zFil5, axis=0, mode='nearest')
    del dzI
    del xFil5, yFil5, zFil5


    ### Optical Flow Calculations ##############################################
    # Equation is v = (A' w A)^-1 A' w b
    # A = -[dxI dyI dzI]
    # b = [dtI]
    # w multiplication is incorporated in the structure tensor inputs above
    # A' w b = -[wdtx wdty wdtz]  (minus sign because of negative sign on A)
    # (A' w A) = [a=wdx2 b=wdxy c=wdxz ; d=wdxy e=wdy2 f=wdyz ; g=wdxz h=wdyz i=wdz2]
    # A^-1 = 1/determinant * [A D G ; B E H ; C F I];
    # using inverse notation from here: https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_3_%C3%97_3_matrices
    #   A = wdy2*wdz2 - wdyz*wdyz; %ei-fh
    #   B = wdyz*wdxz - wdxy*wdz2; %fg-di
    #   C = wdxy*wdyz - wdy2*wdxz; %dh-eg
    #   D = wdxz*wdyz - wdxy*wdz2; %ch-bi
    #   E = wdx2*wdz2 - wdxz*wdxz; %ai-cg
    #   F = wdxy*wdxz - wdx2*wdyz; %bg-ah
    #   G = wdxy*wdyz - wdxz*wdy2; %bf-ce
    #   H = wdxz*wdxy - wdx2*wdyz; %cd-af
    #   I = wdx2*wdy2 - wdxy*wdxy; %ae-bd

    determinant = (wdx2*wdy2*wdz2) + (2*wdxy*wdxz*wdyz) - (wdy2*wdxz**2) - (wdz2*wdxy**2) - (wdx2*wdyz**2)
    vx = -((determinant + np.finfo(float).eps)**-1)*((wdy2*wdz2 - wdyz*wdyz)*wdtx + (wdxz*wdyz - wdxy*wdz2)*wdty + (wdxy*wdyz - wdxz*wdy2)*wdtz)
    vy = -((determinant + np.finfo(float).eps)**-1)*((wdyz*wdxz - wdxy*wdz2)*wdtx + (wdx2*wdz2 - wdxz*wdxz)*wdty + (wdxz*wdxy - wdx2*wdyz)*wdtz)
    vz = -((determinant + np.finfo(float).eps)**-1)*((wdxy*wdyz - wdy2*wdxz)*wdtx + (wdxy*wdxz - wdx2*wdyz)*wdty + (wdx2*wdy2 - wdxy*wdxy)*wdtz)

    del wdtx, wdty, wdtz
    del determinant


    ### Eigenvalues for Reliability ############################################
    # solve det(A^T w A - lamda I) = 0
    # (A' w A) = [a=wdx2 b=wdxy c=wdxz ; d=wdxy e=wdy2 f=wdyz ; g=wdxz h=wdyz i=wdz2]
    # det(A^T w A - lambda I) = (a-lambda)(e-lambda)(i-lambda) + 2*b*c*f - c^2(e-lambda) - b^2(i-lamda) - f^2(a-lambda)

    # Solve for the eigenvalues
    w = np.array([[wdx2, wdxy, wdxz],[wdxy, wdy2, wdyz],[wdxz, wdyz, wdz2]])
    del wdx2, wdxy, wdxz, wdy2, wdyz, wdz2    
    w = np.moveaxis(w,[0,1],[-1,-2])
    w = w.astype(np.complex64) # Allow for complex eignenvalues
    rel = np.linalg.eigvals(w)
    rel = np.real(np.amin(rel,axis=-1))

    ### Return Outputs #########################################################
    return vx, vy, vz, rel

def process_flow(imDir,imName,fileType="SequenceT",spatialDimensions=3,xyzSig=3,tSig=1,wSig=4):
    """
    Parse images for input into calc_flow2D or calc_flow3D.

    Function to organize commands to calc_flow given a single image or image
    sequence for processing. The script can parse two tif formats. If your
    data is another format, this is the function to change to adapt the code
    to your uses.

    OneTif files are assumed to be created using ImageJ when reading metadata.

    ARGS:
    imDir:              Full path to the folder of image(s) to process
    imName:             File name for the image(s) to process, excluding .tif
    fileType:           Either 'OneTif' or 'SequenceT.' In the case 'OneTif',
                        the entire timelapse and all z-slices are assumed to 
                        be saved in one single .tif file, and in the format 
                        generated by ImageJ. In the case of 'SequenceT', the
                        images are assumed to be saved as a sequence, with one
                        tif per timepoint (but all z-slices saved as one tif).
                        In this case, imName must be specified with a wildcard 
                        (.*) for the time label in the file names. For a list of
                        files with names:
                            myexperiment_t000_ch0.tif
                            myexperiment_t001_ch0.tif
                            myexperiment_t002_ch0.tif
                        the correct imName would be 'myexperiment_t.*_ch0'.
    spatialDimensions:  Either 2 (2D) or 3 (3D). Default 3.
    xyzSig:             sigma value for smoothing in all spatial dimensions. Default 3.
                          Larger values remove noise but remove spatial detail.
    tSig:               sigma value for smoothing in the temporal dimension. Default 1.
                          Larger values remove noise but remove temporal detail.
    wSig:               sigma value for Lucas-Kanade neighborhood. Default is 4.
                          Larger values include a larger neighboorhood in the
                          Lucas-Kanade constraint and will smooth over small features.

    RETURNS:
    This function does not have explict returns, but saves several output files.
    Output files are saved in subfolders of the input imDir.
    Output files are prefaced with the input imName.
    *_parameters.csv:   Input parameters as a comma seprated file
    *_vx.tif:           Velocity in the x direction, reported as pixels/frame
    *_vy.tif:           Velocity in the y direction, reported as pixels/frame
    *_vz.tif:           Velocity in the z direction, reported as pixels/frame
    *_rel.tif:          Reliability, the smallest eigenvalue of (A'wA)
                          This is a measure of confidence in the linear algebra solution
                          and can be used to mask the velocities for downstream analysis.
    """

    ### Check Inputs and Set Up Paths ##########################################
    # Check that the directory exists
    imDir = Path(imDir)
    if not imDir.is_dir():
        sys.exit('ERROR: image path \'%s\' does not exist' % imDir)    
    # Check that the image files exist and that the type is correct
    # First get the list of relevant images
    imNamePattern = re.compile(imName + '.tif')
    fileList = []
    files = os.listdir(imDir)
    for f in files:
        m = imNamePattern.fullmatch(f)
        if m:
            fileList.append(f)
    # Now check that the number of files makes sense
    if len(fileList) == 0:
        sys.exit('ERROR: No image files found. imName: ' + imName + ' imDir: ' + str(imDir))
    if fileType=='OneTif':
        if len(fileList) > 1:
            sys.exit('ERROR: Type is OneTif but more than one file was found for imName: ' + imName)
    elif fileType=='SequenceT':
        if len(fileList) < 6*tSig+1: # Minimum requirment for calc_flow
            sys.exit('ERROR: Image sequence found for file name ' + imName + ' only contains ' + str(len(fileList)) + ' files. Minimum 6*tsig+1 ('+ str(6*tSig+1) + ') files required.')
    else:
        sys.exit('ERROR: fileType must be either OneTif or SequenceT.')

    # Make sure files are sorted in numerical order not necessarily ASCII order
    fileList = natsorted(fileList)

    # Must be either 2D or 3D processing
    if spatialDimensions < 2 or spatialDimensions > 3:
        sys.exit('ERROR: Number of spatial dimensions must be either 2 or 3.')

    ### Metadata parsing and parameter saving ##################################
    meta = tf.TiffFile(imDir/ fileList[0]) # Assume first file is representative of whole set
    Ny = meta.pages[0].shape[0]
    Nx = meta.pages[0].shape[1]
    imj = meta.imagej_metadata
    if fileType=='OneTif':
        if not(imj):
            sys.exit('ERROR: fileType is OneTif, but no ImageJ metadata was detected')
        Nt = imj["frames"]
        if spatialDimensions==3:
            Nz = imj["slices"]
        elif spatialDimensions==2:
            Nz = 1
    elif fileType=='SequenceT':
        Nz = len(meta.pages)
        Nt = len(fileList)
        if spatialDimensions==2:
            if Nz != 1:
                sys.exit('ERROR: More than one z-slice detected for 2D processing')
        elif spatialDimensions==3:
            if Nz <= 1:
                sys.exit('ERROR: 3D processing requested but Nz = ' + str(Nz))
    
    NtChunk = 6*tSig+1
    if not(NtChunk%2):
        NtChunk = NtChunk+1
    NtSlice = math.ceil(NtChunk/2)-1

    # Set up the saving folder    
    # Main folder is inside the image directory. 
    if spatialDimensions==3:
        savedir = imDir / 'OpticalFlow3D'
    elif spatialDimensions==2:
        savedir = imDir / 'OpticalFlow2D'
    savedir.mkdir(exist_ok=True)
    # Subfolder is imNameSave.
    imNameSave = imName.replace('.*','')
    savedir = savedir / imNameSave
    savedir.mkdir(exist_ok=True)

    # Save parameters
    param = {'xyzSig': [xyzSig],
            'tiSig': [tSig],
            'wSig': [wSig],
            'Nx': [Nx],
            'Ny': [Ny],
            'Nz': [Nz],
            'Nt': [Nt]}
    param = pd.DataFrame(param)
    savename = imNameSave + '_parameters.csv'
    param.to_csv(savedir / savename, index=False)

    ### Processing Loop ########################################################
    print('Note: regardless of input filenames, the first image = frame 0.')
    print('If your file names start from 0, adjust indexing accordingly for reading the output files.')
    print(' ')

    # Because >6*tSig time frames are necessary for processing, some frames at
    # the start and the end of the timelapse will be ignored.
    for hh in range(0,NtSlice):
        print(str(datetime.now()) + ' - No data will be saved for frame ' + str(hh) + ' to avoid edge effects')
    
    # Loop through the files to be processed
    if fileType=='OneTif': # assuming a tif made with ImageJ containing all z and all t
        imName = imName + '.tif'
        allImages = tf.memmap(imDir / imName) # Does not load images into memory until called later

        if spatialDimensions == 3:    
            for hh in range(0,Nt-NtChunk+1): # If you have enough memory, this could become a parfor loop.
            
                loopStart = datetime.now()
                print(str(datetime.now()) + ' - Processing frame ' + str(hh+NtSlice) + '...')
            
                # Load images
                images = allImages[hh:hh+NtChunk,:]

                # Run the optical flow
                vx,vy,vz,rel = calc_flow3D(images ,xyzSig, tSig, wSig)
            
                # Save this frame
                tstr = str(hh+NtSlice)
                tstr = tstr.zfill(4)
                tf.imwrite(str(savedir / imNameSave) + '_vx_t' + tstr + '.tiff',vx, photometric='minisblack')
                tf.imwrite(str(savedir / imNameSave) + '_vy_t' + tstr + '.tiff',vy, photometric='minisblack')
                tf.imwrite(str(savedir / imNameSave) + '_vz_t' + tstr + '.tiff',vz, photometric='minisblack')
                tf.imwrite(str(savedir / imNameSave) + '_rel_t' + tstr + '.tiff',rel, photometric='minisblack')
            
                del rel, vx, vy, vz, images
            
                framestime = datetime.now()
                print(str(datetime.now()) + ' - Frame ' + str(hh+NtSlice) + ' saved.  Duration: ' + str(framestime-loopStart))
    
        elif spatialDimensions == 2:    
            for hh in range(0,Nt-NtChunk+1):
        
                loopStart = datetime.now()
                print(str(datetime.now()) + ' - Processing frame ' + str(hh+NtSlice) + '...')
            
                # Load images
                images = allImages[hh:hh+NtChunk,:]
            
                # Run the optical flow
                vx,vy,rel = calc_flow2D(images ,xyzSig, tSig, wSig)
            
                # Save this frame
                tstr = str(hh+NtSlice)
                tstr = tstr.zfill(4)
                tf.imwrite(str(savedir / imNameSave) + '_vx_t' + tstr + '.tiff',vx, photometric='minisblack')
                tf.imwrite(str(savedir / imNameSave) + '_vy_t' + tstr + '.tiff',vy, photometric='minisblack')
                tf.imwrite(str(savedir / imNameSave) + '_rel_t' + tstr + '.tiff',rel, photometric='minisblack')

                del rel, vx, vy, images
            
                framestime = datetime.now()
                print(str(datetime.now()) + ' - Frame ' + str(hh+NtSlice) + ' saved.  Duration: ' + str(framestime-loopStart))
    
        else:
            sys.exit('ERROR: Spatial Dimension must be 2 or 3.')
    
    elif fileType=='SequenceT': # assuming 1 tif per timepoint that is a z-stack (or slice for 2D)    
        if spatialDimensions == 3:    
            for hh in range(0,Nt-NtChunk+1):
        
                loopStart = datetime.now()
                print(str(datetime.now()) + ' - Processing frame ' + str(hh+NtSlice) + '...')
            
                # Load images
                images = tf.imread(imDir / fileList[hh])
                images = np.append([images],[tf.imread(imDir / fileList[hh+1])],axis=0)
                for jj in range(2,NtChunk):
                    images = np.append(images,[tf.imread(imDir / fileList[hh+jj])],axis=0)
            
                # Run the optical flow
                vx,vy,vz,rel = calc_flow3D(images ,xyzSig, tSig, wSig)
            
                # Save this frame
                tstr = str(hh+NtSlice)
                tstr = tstr.zfill(4)
                tf.imwrite(str(savedir / imNameSave) + '_vx_t' + tstr + '.tiff',vx, photometric='minisblack')
                tf.imwrite(str(savedir / imNameSave) + '_vy_t' + tstr + '.tiff',vy, photometric='minisblack')
                tf.imwrite(str(savedir / imNameSave) + '_vz_t' + tstr + '.tiff',vz, photometric='minisblack')
                tf.imwrite(str(savedir / imNameSave) + '_rel_t' + tstr + '.tiff',rel, photometric='minisblack')
            
                del rel, vx, vy, vz, images
            
                framestime = datetime.now()
                print(str(datetime.now()) + ' - Frame ' + str(hh+NtSlice) + ' saved.  Duration: ' + str(framestime-loopStart))

        elif spatialDimensions == 2:    
            for hh in range(0,Nt-NtChunk+1):
        
                loopStart = datetime.now()
                print(str(datetime.now()) + ' - Processing frame ' + str(hh+NtSlice) + '...')
            
                # Load images
                images = tf.imread(imDir / fileList[hh])
                images = np.append([images],[tf.imread(imDir / fileList[hh+1])],axis=0)
                for jj in range(2,NtChunk):
                    images = np.append(images,[tf.imread(imDir / fileList[hh+jj])],axis=0)
            
                # Run the optical flow
                vx,vy,rel = calc_flow2D(images ,xyzSig, tSig, wSig)
            
                # Save this frame
                tstr = str(hh+NtSlice)
                tstr = tstr.zfill(4)
                tf.imwrite(str(savedir / imNameSave) + '_vx_t' + tstr + '.tiff',vx, photometric='minisblack')
                tf.imwrite(str(savedir / imNameSave) + '_vy_t' + tstr + '.tiff',vy, photometric='minisblack')
                tf.imwrite(str(savedir / imNameSave) + '_rel_t' + tstr + '.tiff',rel, photometric='minisblack')
            
                del rel, vx, vy, images
            
                framestime = datetime.now()
                print(str(datetime.now()) + ' - Frame ' + str(hh+NtSlice) + ' saved.  Duration: ' + str(framestime-loopStart))
    
        else:
            sys.exit('ERROR: Spatial Dimension must be 2 or 3')
    
    # Because >6*tSig time frames are necessary for processing, some frames at
    # the start and the end of the timelapse will be ignored.
    for hh in range(Nt-NtSlice,Nt):
        print(str(datetime.now()) + ' - No data will be saved for frame ' + str(hh) + ' to avoid edge effects')

