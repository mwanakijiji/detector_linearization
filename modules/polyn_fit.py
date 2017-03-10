# These are functions concerned with the actual polynomial fitting to pixel responses

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.modeling import models, fitting
import ipdb

# fit polynomial to pixel responses
def fit_polyn_each_pix(intTimesDarksPass,
                       stemDarkSubtedFlats,
                       patternFilenamesPass,
                       numberStartPass,
                       numberStopPass,
                       ySizePass,
                       xSizePass,
                       orderPolynPass,
                       correction=True
):

    pfit = fitting.LinearLSQFitter()
    
    # initialize a cube to hold dark-subtracted flats at different DITs
    dataCube = np.zeros((numberStopPass-numberStartPass+1,ySizePass,xSizePass))

    # initialize a cube to hold the polynomial coefficients of the pixel responses
    coeffCube = np.zeros((orderPolynPass+1,ySizePass,xSizePass))

    # initialize a cube to hold the polynomial coefficients of empirical --> linear mapping
    correctionCube = np.zeros((orderPolynPass+1,ySizePass,xSizePass))  

    for s in range(numberStartPass,numberStopPass+1): # loop over frames with the same DIT
        image, header = fits.getdata(stemDarkSubtedFlats+
                                     patternFilenamesPass+
                                     str("{:0>2d}".format(s))+
                                     '.fits',
                                     0,
                                     header=True)

        dataCube[s-numberStartPass,:,:] = image # put image into cube



    ############## START GENERATNG PIX RESPONSE CORRECTION  ############
    # loop over columns
    if (correction==True):
        for currentCol in range(0+5,xSizePass-5): # boundary to avoid screwy pixels on the edge
            print('Finding polynomial correction to pixels in col '+str(currentCol)+'...')

            # establish ideal 'linear' response
            modelsInit_linearIdeal = models.Polynomial1D(1, # linear fit
                                         n_models=ySizePass,
                                         model_set_axis=1)
            newModelsCurrentCol_linearIdeal = pfit(modelsInit_linearIdeal,
                                   intTimesDarksPass[2][0:3], # DIT axis
                                   dataCube[0:3,:,currentCol])  # along one column

            allIntTimes = np.array([intTimesDarksPass[2][numberStartPass:numberStopPass+1]])
            # find ideal 'linear' points
            idealYvals = newModelsCurrentCol_linearIdeal.c0+np.multiply(newModelsCurrentCol_linearIdeal.c1,
                                                                               np.tile(np.transpose(allIntTimes),(1,np.ma.size(newModelsCurrentCol_linearIdeal.c1))))

            for currentRow in range(0,ySizePass):
                result_rev = np.polyfit(dataCube[:,currentRow,currentCol], idealYvals[:,currentRow], orderPolynPass)
                correctionCube[:,currentRow,currentCol] = np.flipud(result_rev)
            
            returnCube = correctionCube

    ############## END GENERATNG PIX RESPONSE CORRECTION  ############

    ############## START CURVE-FITTING TO PIX RESPONSES ############
    # loop over columns
    if (correction==False): # i.e., if I just want polynomials to describe the pixel responses
        for currentCol in range(0,xSizePass):
            print('Fitting polynomials to pixel responses in col '+str(currentCol)+'...')

            # initialize set of models for an entire column    
            modelsInit = models.Polynomial1D(orderPolynPass,
                                     n_models=ySizePass,
                                     model_set_axis=1)
            newModelsCurrentCol = pfit(modelsInit,
                                   intTimesDarksPass[2][numberStartPass:numberStopPass+1], # DIT axis
                                   dataCube[:,:,currentCol])  # along one column

            # loop over coeff numbers
            for coeffNum in range(0,orderPolynPass+1):
                coeffCube[coeffNum,:,currentCol] = newModelsCurrentCol.param_sets[coeffNum,0,:]

        returnCube = coeffCube
    ############## STOP CURVE-FITTING TO PIX RESPONSES ############

    return returnCube, dataCube

