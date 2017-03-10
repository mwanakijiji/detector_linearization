# fcn to return all values of a header card in a series of FITS files

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# determine header card values for a series of FITS files
# input:
# [0]: directory stem of files
# [1]: date string (i.e., 20170310)
# [2]: camera ('n' or 'lm')
# [3]: number of first frame
# [4]: number of last frame
# [5]: string corresponding to the card in the header
# output:
# [0]: array of all frame numbers
# [1]: array of card values corresponding to all frames
# [2]: array of unique card values 
def return_header_card_vals(stemPass, 
                            dateStringPass,
                            camPass,
                            firstFramePass,
                            lastFramePass,
                            headerCardString
                           ):
    
    cardValArrayPass = np.zeros(lastFramePass-firstFramePass+1)
    frameNumArrayPass = np.copy(cardValArrayPass)
    
    for p in range(firstFramePass,lastFramePass+1):
        hdul = fits.open(stemPass+camPass+'_'+dateStringPass+'_'+str("{:0>6d}".format(p))+'.fits')
        cardValArrayPass[p-firstFramePass] = hdul[0].header[headerCardString]
        frameNumArrayPass[p-firstFramePass] = p
        
    cardValArrayPassUnique = np.unique(cardValArrayPass) # unique card values
    
    return frameNumArrayPass, cardValArrayPass, cardValArrayPassUnique


# make a median dark and flat for each DIT
# input:
# [0]: directory stem of files to read in
# [1]: " " to write files out in
# [2]: array of integration times returned by return_header_card_vals
# [3]: string denoting nature of the frame (i.e., 'flat' or 'dark')
# [4]: string ref to thing being done (i.e.,'median')
# [5]: date string (i.e., 20170310)
# [6]: camera ('n' or 'lm')
# [7]: number of rows in frame
# [8]: number of cols
# output:
# none; just writes out frames to file
def write_median_frames(stemReadPass,
                        stemWritePass,
                        intTimesPass,
                        descriptorStringPass,
                        typeStringPass,
                        dateStringPass,
                        camPass,
                        ySizePass,
                        xSizePass,):

    print('Warning! This is wired to clobber frames. Type c to continue.')

    print('Making median of '+descriptorStringPass+'s for integration time')
    # note: 2017 Feb 09 NOMIC data has frames
    # 282580 --> 282879: darks with increasing DITs (10 frames per DIT)
    # 282880 --> 283179: dome flats with increasing DITs ( "  " )

    for l in range(0,np.ma.size(intTimesPass[2])): # loop over unique DITs
        thisDIT = intTimesPass[2][l]
        frameNumsOfInterest = intTimesPass[0][np.where(intTimesPass[1] == thisDIT)]
    
        print(str(thisDIT)+' sec\nwhich includes frame numbers')
        print(frameNumsOfInterest)
    
        # initialize a cube
        dataCubeForMedian = np.zeros((np.ma.size(frameNumsOfInterest),ySizePass,xSizePass))

        for s in range(0,np.ma.size(frameNumsOfInterest)): # loop over frames with the same DIT
            image, header = fits.getdata(stemReadPass+
                                         camPass+'_'+
                                         dateStringPass+'_'+
                                         str("{:0>6d}".format(frameNumsOfInterest[s].astype(int)))+
                                         '.fits',
                                         0,
                                         header=True)
            #image = np.squeeze(hdul[0].data.copy())                                   
            dataCubeForMedian[s,:,:] = image # put image into cube
            
        # take median
        thisMedian = np.median(dataCubeForMedian, axis=0)

        # clear cube
        del dataCubeForMedian
        
        # write out median
        fits.writeto(stemWritePass+
                     camPass+'_'+
                     descriptorStringPass+'_'+
                     typeStringPass+'_'+
                     dateStringPass+'_'+
                     str("{:0>2d}".format(l))+
                     '.fits',
                     thisMedian, clobber=True)
    
        print('---')


# dark-subtract flat frames
# inputs:
# [0]: integration times array for the darks
# [1]: " " flats
# [2]: stem for the darks directory
# [3]: dark filename pattern before the '[number].fits'
# [4]: " " flats " "
# [5]: flats " "
# [6]: starting number in the median frames
# [7]: ending number " "
# [8]: camera ('n' or 'lm')
# [9]: date string
# [10]: descriptor to put in written out filenames
# [11]: stem for the writeouts directory
def dark_subtract(intTimesDarksPass,
                  intTimesFlatsPass,
                  stemDarksPass,
                  patternDarksPass,
                  stemFlatsPass,
                  patternFlatsPass,
                  numberStartPass,
                  numberStopPass,
                  camPass,
                  dateStringPass,
                  descriptorStringPass,
                  stemWritePass):

    print('Dark-subtracting flats...')
    print('Warning! Clobber is set to true')
    
    # sanity check: are DITs the same between the flats and darks?
    timeDiffs = np.subtract(intTimesFlatsPass[2][:],intTimesDarksPass[2][:])

    if max(np.abs(timeDiffs)!=0):
        print("Error: inconsistent DITs between flats, darks.")
        exit()
    
    for t in range(numberStartPass,numberStopPass+1): # loop over frames with the same DIT
        imageDark, header = fits.getdata(stemDarksPass+
                                         patternDarksPass+
                                         str("{:0>2d}".format(t))+
                                         '.fits',
                                         0,
                                         header=True)
      
        imageFlat, header = fits.getdata(stemFlatsPass+
                                     patternFlatsPass+
                                     str("{:0>2d}".format(t))+
                                     '.fits',
                                     0,
                                     header=True)

        # dark-subtraction
        imageDiff = np.subtract(imageFlat,imageDark)

        # write out
        fits.writeto(stemWritePass+
                 camPass+'_'+
                 descriptorStringPass+'_'+
                 dateStringPass+'_'+
                 str("{:0>2d}".format(t))+
                 '.fits',
                 imageDiff, clobber=True)
