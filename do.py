#from modules import header_info
#from header_info
import numpy as np
import matplotlib.pyplot as plt
import ipdb
from modules.constants import *
from modules import header_info
from modules import polyn_fit
#from modules.header_info import return_header_card_vals

#import ipdb; ipdb.set_trace()

# get integration times for the darks and flats
intTimesDarks = header_info.return_header_card_vals(stem+'copies_raw_data/',
                                                    dateString,
                                                    cam,
                                                    firstFrameDarks,
                                                    lastFrameDarks,
                                                    'EXPTIME')
intTimesFlats = header_info.return_header_card_vals(stem+'copies_raw_data/',
                                                    dateString,
                                                    cam,
                                                    firstFrameFlats,
                                                    lastFrameFlats,
                                                    'EXPTIME')

'''
# write out median darks and flats for each DIT
header_info.write_median_frames(stem+'copies_raw_data/',
                                stem+'step01_dark_flat_medians/',
                                intTimesDarks,'dark','median',dateString,cam,ySize,xSize)
header_info.write_median_frames(stem+'copies_raw_data/',
                                stem+'step01_dark_flat_medians/',
                                intTimesFlats,'flat','median',dateString,cam,ySize,xSize)
'''

# subtract the darks from the flats
'''
header_info.dark_subtract(intTimesDarks,
                          intTimesFlats,
                          stem+'step01_dark_flat_medians/',
                          'n_dark_median_'+dateString+'_',
                          stem+'step01_dark_flat_medians/',
                          'n_flat_median_'+dateString+'_',
                          0,
                          27,
                          cam,
                          dateString,
                          'darkSubtractedFlat',
                          stem+'step02_dark_subtracted_flats/')
'''

# fit a 3rd-order polynomial to pixel responses
coeffCube3, flatCube3 = polyn_fit.fit_polyn_each_pix(intTimesDarks,
                             stem+'step02_dark_subtracted_flats/',
                             cam+'_darkSubtractedFlat_'+dateString+'_',
                             0,
                             5,
                             ySize,
                             xSize,
                                                     3,correction=False)

ipdb.set_trace()
'''
# fit a 2nd-order polynomial to pixel responses
coeffCube2, flatCube2 = polyn_fit.fit_polyn_each_pix(intTimesDarks,
                             stem+'step02_dark_subtracted_flats/',
                             cam+'_darkSubtractedFlat_'+dateString+'_',
                             0,
                             5,
                             ySize,
                             xSize,
2,correction=False)
'''

# generate a 2nd-order correction

######################################################################
# PLOTTING

xRandArray = [30,30,98,98,99,61,61,224]
yRandArray = [30,98,30,98,63,64,190,61]
strings = ['Test1','Test2','Test3','Test4','POS1','StarNOD1','StarNOD2','POS2']

'''
ipdb.set_trace()
# plot response corrections
for p in range(0,8):
    xRand = xRandArray[p] #int(np.random.rand()*256)
    yRand = yRandArray[p] #int(np.random.rand()*256)
    #plt.scatter(intTimesDarks[2][0:6],flatCube3[:,yRand,xRand],s=50,color='k')
    # intTimesDarks[2][0:6] # x-vals from empirical data
    empirical = flatCube3[0:6,yRand,xRand]
    ideal = coeffCube3[0,yRand,xRand]+coeffCube3[1,yRand,xRand]*empirical+coeffCube3[2,yRand,xRand]*np.power(empirical,2)+coeffCube3[3,yRand,xRand]*np.power(empirical,3)
    plt.plot(empirical,ideal,color='b',marker='o')
    plt.plot(np.arange(np.min(empirical),np.max(empirical)),np.arange(np.min(empirical),np.max(empirical)),'--k')
    plt.xlabel('Empirical (ADU)')
    plt.ylabel('Linearized (ADU)')
    #plt.xlim([0,0.05])
    plt.title('3rd-order Correction to Pixel ('+str(strings[p])+'; x='+str(xRand)+', y='+str(yRand)+')')
    plt.show()
    plt.plot(empirical,np.divide(np.subtract(ideal,empirical),empirical))
    plt.title('3rd-order Correction Error to Pixel ('+str(strings[p])+'; x='+str(xRand)+', y='+str(yRand)+')')
    plt.xlabel('Empirical (ADU)')
    plt.ylabel('Error')
    plt.show()
'''

'''
# plot pixel responses
for p in range(0,8):
    xRand = xRandArray[p] #int(np.random.rand()*256)
    yRand = yRandArray[p] #int(np.random.rand()*256)
    plt.scatter(intTimesDarks[2][0:6],flatCube2[:,yRand,xRand],s=50,color='k')
    times = np.arange(0,50)*0.001
    # intTimesDarks[2][0:6] # x-vals from empirical data
    plt.plot(times,coeffCube2[0,yRand,xRand]+
             coeffCube2[1,yRand,xRand]*times+
             coeffCube2[2,yRand,xRand]*np.power(times,2), #coeffCube[3,yRand,xRand]*np.power(times,3),
             color='r')
    plt.plot(times,coeffCube3[0,yRand,xRand]+
             coeffCube3[1,yRand,xRand]*times+
             coeffCube3[2,yRand,xRand]*np.power(times,2)+coeffCube3[3,yRand,xRand]*np.power(times,3),
             color='b')
    plt.xlabel('Integration Times (sec)')
    plt.ylabel('ADU')
    plt.xlim([0,0.05])
    plt.title('Polyn Fit to Pixel ('+str(strings[p])+'; x='+str(xRand)+', y='+str(yRand)+') Response\n(black: empirical; red: 2nd order polyn; blue: 3rd order)')
    plt.show()
'''

# coefficient plots, with dots to indicate sample pixels

fig, ax = plt.subplots()
coeffCubeToPlot = coeffCube3
titleString = '3rd-order correction, coeff $c_{3}$'
cax = ax.imshow(coeffCubeToPlot[3,:,:], origin='lower')
ax.scatter(xRandArray,yRandArray,s=20,color='r')
for i, txt in enumerate(strings):
    ax.annotate(txt, (xRandArray[i],yRandArray[i]))
cbar = fig.colorbar(cax)
ax.set_title(titleString)
plt.show()


