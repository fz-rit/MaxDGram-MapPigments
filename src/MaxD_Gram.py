'''
+
=======================================================================

 NAME:
      MaxD_Gram

 DESCRIPTION:
Extracts end members from a spectral image using the MaxD approach.  Also estimates the
Material diversity in a scene using the Gram matrix approach, both global and local.

TO CALL:
   endmembers, endmembers_index, volume = maximumDistance(data, num, mnf_data, gram)

 USES:
numpy

 PARAMETERS:
    data = 2D data [npixels, nbands]
    num = number of endmembers to be calculated (choose more than expected to find)
    mnf_data = [MNF data cube or 0]; if MNF data is not available (i.e., mnf_data = 0),
        code will assign img as mnf_data

 KEYWORDS:


 RETURNS:
maximumDistance -  endmembers, endmembers_index, volume
calcGram; calcGramLocal - gram (Grammian, the volume of the parallelotope enclosed by the end members)

 NOTES:
- original python code developed by Tania Kleynhans based on papers by Amanda Zieman and Dave Messinger

 HISTORY:
08/30/2023 - created

=======================================================================
-
'''

import numpy as np
from numpy import matlib as mb


def maximumDistance(data, num, mnf_data, gram):
    # data = 2D data [npixels, nbands]
    # num = number of endmembers to be calculated (choose more than expected to find)
    # if MNF data is not available, code will assign img as mnf_data
    print('---> In MaxD extracting endmembers and Grammian ...')
    if mnf_data == 0:
        mnf_data = data

    data = np.transpose(data)
    data2 = np.transpose(mnf_data)

    # find data size
    num_bands = data.shape[0]
    num_pix = data.shape[1]

    # calculate magnitude of all vectors to find min and max
    magnitude = np.sum(np.square(data), axis=0)
    idx1 = np.int_(np.where(magnitude == np.max(magnitude))[0])
    idx2 = np.int_(np.where(magnitude == np.min(magnitude))[0])

    # create empty output arrays for endmembers
    endmembers = np.zeros([num_bands, num])
    endmembers_index = np.zeros([1, num])

    # assign largest and smallest vector as first and second endmembers
    endmembers[:, 0] = np.transpose(data[:, idx1])
    endmembers[:, 1] = np.transpose(data[:, idx2])
    endmembers_index[0, 0] = idx1
    endmembers_index[0, 1] = idx2

    data_proj = np.matrix(data2)
    identity_matrix = np.identity(num_bands)

    # create array for volume of determinant of Gram matrix
    volume = np.zeros([num])

    loop = np.arange(3, num + 1)
    for i in loop:
        diff = []
        pseudo = []
        # calc difference between endmembers
        diff = np.matrix(data_proj[:, idx2] - data_proj[:, idx1])
        # caclualte pseudo inverse of difference vector
        pseudo = np.linalg.pinv(diff)
        data_proj = np.matmul((identity_matrix - np.matmul(diff, pseudo)), data_proj)

        idx1 = idx2
        diff_new = np.sum(np.square((np.matmul(data_proj[:, idx2], np.ones([1, num_pix])) - data_proj)), axis=0)

        # find ne maximum distance for next endmember
        idx2 = np.int_(np.where(diff_new == np.max(diff_new))[1])

        #print('DEBUG: idx2: ', idx2,np.size(idx2))
        ###
        # DWM: looks like there may be cases where idx2 has more than one element, i.e., there are
        # two elements of diff_new that are equal to the max.  In that case, just grab the first one
        ###
        if np.size(idx2) > 1:
            idx2 = idx2[0]

        # assign to endmember file
        endmembers[:, i - 1] = np.transpose(data[:, idx2])
        endmembers_index[0, i - 1] = idx2

        if gram == 'local':
            # calculate local gram matrix
            loc_gram = calcGramLocal(endmembers, i)
            volume[i - 1] = np.sqrt(np.abs(np.linalg.det(loc_gram)))

        elif gram == 'general':
            # calculate general gram matrix
            gen_gram = calcGramGeneral(endmembers[:, 0:i])
            volume[i - 1] = np.sqrt(np.abs(np.linalg.det(gen_gram)))

    return endmembers, endmembers_index, volume


def calcGramGeneral(data_endmembers):
    # calculate gram matrix = V^T * V
    gram = np.matmul(np.transpose(data_endmembers), data_endmembers)

    return gram


def calcGramLocal(data_endmembers, iteration):
    # use only endmembers already calculated
    data_endmembers = data_endmembers[:, 0:iteration]

    # calculate the Gram matrix based on local information (points nearest to mean)
    # num_bands = data_endmembers.shape[0]
    num_pix = data_endmembers.shape[1]

    # create mean vector
    mean_spec = np.mean(data_endmembers, axis=1)

    # calculate normalized difference between mean vector and endmembers and find closest vector to mean vector
    diffdist = np.linalg.norm(np.transpose(mb.repmat(mean_spec, num_pix, 1)) - data_endmembers, axis=0)
    min_idx = np.int(np.where(diffdist == np.min(diffdist))[0])

    # create index of rows to keep
    index = np.ones([num_pix])
    # keep all but min distance one
    index[min_idx] = 0
    # find index of all nonzero entires
    keep_idx = np.squeeze(np.where(index == 1))
    nearpix = data_endmembers[:, keep_idx]

    # calculate local Gram
    num_neighbors = nearpix.shape[1]
    # gram = np.zeros([num_neighbors, num_neighbors])
    diff_matrix = nearpix - np.transpose(mb.repmat(mean_spec, num_neighbors, 1))

    gram = np.matmul(np.transpose(diff_matrix), diff_matrix)
    print('<--- done')

    return gram

