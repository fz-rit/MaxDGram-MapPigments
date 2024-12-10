'''
+
=======================================================================

 NAME:
      spectral_tools

 DESCRIPTION:
simple spectral tools, starting with spectral angle calculation

 USES:
spectral
spectral.io.envi
numpy
matplotlib
openENVI

 PARAMETERS:


 KEYWORDS:


 RETURNS:

 NOTES:
Tools include:
- vector magnitude
- spectral angle between two vectors
- spectral angle map
- normalized band difference
- gets a band index from wavelength value in an image band array
- makes a three color image for display; can be true RGB or false color

 HISTORY:

04/20/2023: D. Messinger - created


=======================================================================
-
'''

import numpy as np
import math

#print ('in spectral tools')

###
# function definition to compute magnitude of the vector
###
def magnitude(vector):
    return math.sqrt(sum(pow(element, 2) for element in vector))
 

###
# Spectral angle between two pixels
###
def spec_ang(pix1,pix2): 
    #compute the spectral angle between the two pixels passed to this code, result is in radians.
    SA = np.arccos((np.dot(pix1,pix2))/(magnitude(pix1) * magnitude(pix2)))
    return(SA)



###
# Spectral Angle Map
###
def spectral_angle_map(image,target):
    ###
    # image - format of (nrows, ncols, nbands)
    # target - format of (nbands)
    # returns a single band image of the spectral angle between each pixel in image and the target

    print('....> in spectral_angle_map')

    # create the new image
    nrows = image.nrows
    ncols = image.ncols
    nbands = image.nbands

    print('nrows, ncols, nbands', nrows,ncols,nbands)

    sam_image = np.ndarray([nrows,ncols])
    print(sam_image.shape)
#pix1 = np.reshape(myimg[75,75,:],nbands)
    for irow in range(nrows):
        for icol in range(ncols):
           sam_value = np.arccos((np.dot(np.reshape(image[irow,icol,:],nbands),target))/(magnitude(np.reshape(image[irow,icol,:],nbands)) * magnitude(target)))
           sam_image[irow,icol] = sam_value
    else:
     print('done with matrix')#
     
     # debug
#     print(sam_image.shape)
     print('<.... done')
     print('')

     return(sam_image)


###
# Normalized Band Difference (NBD)
# computes a new image that is (Band 1 - Band 2) / (Band 1 + Band 2)
###
def NBD(image, band1, band2):
    ###
    # image - format of (nrows, ncols, nbands)
    # band1 - first band of formula
    # band2 - second band of formula
    # returns nbd_image, a single band image of the normalized band difference (NBD)

    print('....> in NBD_image')
#    print('in NBD_image')
    print('Band 1: ', band1)
    print('Band 2: ', band2)

# create the new image
    nbd_image = np.ndarray([image.nrows,image.ncols])

# compute the new image
    nbd_image = (image[:,:,band1] - image[:,:,band2])/(image[:,:,band1] + image[:,:,band2])
# final image should only be nrows by ncols so reshape it
    nbd_image = np.reshape(nbd_image,(image.nrows,image.ncols))
# debug
#    print(nbd_image.shape)
    print('<..... done')
    #print('<....')
    print('')
    return (nbd_image)

###
# get_band_index
# returns the band array index of the band array at a desired wavelength
# NOTE: assumes bandarray and WL are in the same units, e.g., nanometers
#
def get_band_index(bandarray, WL):
    ###
    # bandarray - array of wavelength information
    # WL - wavelength of interest
    # returns band_index, the index of the array at that wavelength
    # NOTE: ASSUMES bandarray and WL are in the same units!!!!

    print('....> in get_band_index')
    #print('WL: ', WL)

    #set up array to find minimum
    nbands = np.size(bandarray)
    temp_array= np.ones(nbands) * WL

    # find the minimum value
    min_array = np.abs(bandarray - temp_array)
    index_array = np.where(min_array == np.amin(min_array))
    band_index = int(index_array[0])

    #print('Associated band index: ', band_index)
    print('<..... done')
    #print('<....')
    print('')
    return (band_index)

###
# Make a three color image
# [b1, b2, b3] will be indices 0, 1, 2 for rendering as a color image
# note that in python displaying an RGB image, imshow assumes indices [0, 1, 2] are [red, green, blue]
###
def make_color_image(image, b1, b2, b3):

    nrows = image.shape[0]
    ncols = image.shape[1]

    redimg = image[:, :, b1]
    redimg = redimg.reshape(nrows, ncols)
    #
    greenimg = image[:, :, b2]
    greenimg = greenimg.reshape(nrows, ncols)
    #
    blueimg = image[:, :, b3]
    blueimg = blueimg.reshape(nrows, ncols)
    #
    # NOTE: displaying an RGB image imshow assumes indices [0,1,2] are [red, green, blue]
    #
    rgbimg = np.ndarray((nrows, ncols, 3))
    rgbimg[:, :, 2] = blueimg
    rgbimg[:, :, 1] = greenimg
    rgbimg[:, :, 0] = redimg
    return (rgbimg)