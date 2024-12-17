'''
+
=======================================================================

 NAME:
      pigment_map_SAM_EM

 DESCRIPTION:
Map pigments in a HSI using endmembers and spectral angle mapping; endmembers
are extracted from the image using the MaxD and Gram Matrix approach

 USES:
spectral
spectral.io.envi
numpy
matplotlib
openENVI
spectral_tools - my spectral tools in an external file
MaxD_Gram - external code

 PARAMETERS:

 KEYWORDS:
needs to know the location of the input file and how many EM's to keep

 RETURNS:
saves the colored class map as a pdf with colors related to EM's; spectra
are also saved as a pdf

 NOTES:


 HISTORY:
08/30/2023: D. Messinger - created, based on pigment_SAM.py



=======================================================================
-
'''

from spectral import *
import spectral.io.envi as envi
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from src.MaxD_Gram import *

import time


#define function to open up an ENVI file and return the image
#from openENVI import *
# my spectral tools
from src.spectral_tools import *

###
# starting stuff
###
ifigure = 1
start_time = time.time()
start_hour = time.gmtime().tm_hour
start_min = time.gmtime().tm_min
start_sec = time.gmtime().tm_sec
print('Starting time [GMT]: ', start_hour,':', start_min,':',start_sec)

###
# input image file location and name
###

infolder = r'C:\Users\fzhcis\Documents\projects\from_DavidM\for_Fei\Data\Symeon\\'
image_file = infolder + 'Symeon_VNIR_cropped'
outfolder = infolder
outfile = outfolder+'Symeon_VNIR_MaxD_class_map_v2.pdf'
specfolder = outfolder
specfile = outfolder+'junk.pdf'
saveimages = False

###
# enter number of endmembers to start off with (this should be more than the expected endmembers)
# and the number to keep in the end
###
num_EM = 15
num_to_keep = 6

###
# open the image file
###
#image_file = infolder +'refl_img'
#image_file = infolder +'GM_HSI_3R_test_chip_4_1600x2000'
print ('opening image file: ', image_file)
image = open_image(image_file+'.hdr').load()
nrows = image.nrows
ncols = image.ncols
nbands = image.nbands
print('IMAGE rows, cols, bands: ', image.nrows, image.ncols, image.nbands)
print('')
#
# make an RGB picture
#
# extract the red green and blue bands to make a picture
#

bands = np.array(image.nbands)
bands = image.bands.centers
# iblue = get_band_index(bands,475.0)
# igreen = get_band_index(bands,535.0)
# ired = get_band_index(bands,650.0)
# for SWIR
iblue = get_band_index(bands,1034.0)
igreen = get_band_index(bands,1195.0)
ired = get_band_index(bands,1600.0)

if str(image.bands.centers) == 'None' :
    print('creating band centers array')
    bands = np.arange(nbands) + 1

###
# for MSI, max normalize the image
###
#image[:,:,:] = image[:,:,:]/np.amax(image[:,:,:])

# ####
# #For MSI, blue green and red are 2, 4, 6
# ####
# iblue = 2
# igreen = 4
# ired = 6

###
# make the RGB image
###
rgbimg = make_color_image(image, ired, igreen, iblue)

plt.figure(ifigure)
picture = plt.imshow((rgbimg/np.amax(rgbimg)))
plt.title('RGB')
#plt.show()
ifigure += 1

###
# extract the endmembers with MaxD / Gram matrix
###
# reshape the image to pass to MaxD
image2D = np.reshape(image, [nrows * ncols, nbands], order="F")

# choose between using general gram matrix or local gram matrix
gram = 'general'
# gram = 'local'

# if you have mnf data of the image, set mnf_data to that, else code will use image data as mnf_data
mnf_data = 0

# number of EM's to extract
print(f'Extracting {num_EM} EMs: ')

endmembers, endmembers_index, volume = maximumDistance(image2D, num_EM, mnf_data, gram)

# normalize volume
volume_norm = volume / sum(volume)

# plot volume function
x = np.arange(3, num_EM + 1)
plt.figure(ifigure)
plt.plot(x, volume_norm[2:])
plt.xlabel('Number of endmembers')
plt.ylabel('Normalized estimated volume')
plt.title('Grammian Volume Function')
#plt.show()
ifigure +=1


#choose number of endmembers to keep based on graph
# value = input("How many EM's to keep? [integer]:\n")
# value = int(value)
# print(f'Keeping {value} endmembers....')
#num_to_keep = value

print(f'---> Keeping {num_to_keep} endmembers ....')
endmembers[:,0:num_to_keep-1]
endmembers_index[:,0:num_to_keep-1]

###
# build the class_spec array from the EMs
###
# class_spec = np.zeros((library.spectra.shape[0],library.spectra.shape[1]), float)
# for i,c in enumerate(library.spectra) :
#     class_spec[i] = library.spectra[i]
lib_nspec = num_to_keep
class_spec = np.zeros((num_to_keep,nbands), float)
for i in range(num_to_keep):
    class_spec[i,:] = endmembers[:,i]

###
# plot the spectra used for classification
#
# first, set up the color map for plotting and making the class map so they match up
#
#-> import the color map
mycmap = plt.get_cmap('jet', lib_nspec)
#-> generate color array
newcolors = mycmap(np.arange(0,mycmap.N))
#-> cast it as a ListedColormap object with attributes
newcmp = ListedColormap(newcolors)

###
# plot the class spectra
###
plt.figure(ifigure)
#plt.plot(library.bands.centers,class_spec[5,:],color='black',label = library.names[5])
for j in range(lib_nspec) :
    EM_names = 'EM: '+str(j+1)
    plt.plot(bands,class_spec[j,:],color=newcmp.colors[j], label = EM_names)

plt.xlabel('Wavelength (nm)')
plt.ylabel('Reflectance')
plt.title('Endmember Spectra')
plt.legend(loc='best')
if saveimages:
    plt.savefig(specfile,dpi = 600)
#plt.show()
ifigure += 1

###
# run spectral angle classification with the library
###
print('')
print('--> making the class map...')
angles = spectral_angles(image, class_spec)

###
# identify the smallest angle for each pixel
###
class_map = np.argmin(angles,2)
print('<--- done')

print('')
print("RUN TIME --- %s seconds ---" % (time.time() - start_time))
###
# show the class map image with the same color map as used to plot the spectra
###
plt.figure(ifigure)
picture = plt.imshow(np.flip(np.flip(class_map,axis=1),axis=0),cmap = newcmp)
plt.title('Class Map from EMs: HSI')
#plt.show()
ifigure +=1
if saveimages:
    #plt.imsave(outfolder+'WS_FF_classmap_4classes.pdf',class_map)
    plt.savefig(outfile,dpi=600)

print("--- %s seconds ---" % (time.time() - start_time))
plt.show()
