from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
import numpy as np
import matplotlib.pyplot as plt
import glob, random, fitsio

from astropy.table import Table, join, vstack
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u
from astropy.cosmology import LambdaCDM, z_at_value
from astropy.io import fits

import pandas as pd

from functions import *

save_directory = '/pscratch/sd/c/clamman/abacus/all_halos_info/'

halo_sims_paths = np.sort(glob.glob(
        '/global/cfs/projectdirs/desi/cosmosim/Abacus/AbacusSummit_base_c000_ph*'))


# going through abacus halos and taking the largest halos, matching 2x the LRG density

def get_halo_info(halo_paths):
    
    # going through each slab and getting the right number of halos for each
    all_info=[]
    n_in_slab = 76800/2  # rough number to match LRG density x2 /2

    for h in halo_paths:
        #print('Working on', h.split('info_')[-1].strip('.asdf'), '/',len(halo_pahts))
        cat = CompaSOHaloCatalog(str(h), fields=['x_L2com', 'v_L2com', 'sigman_L2com', \
                                                 'sigman_eigenvecsMin_L2com', 'sigman_eigenvecsMid_L2com', \
                                                 'sigman_eigenvecsMaj_L2com', 'N'])
        halos = cat.halos
        
        halost = Table()
        halost['x_L2com'] = halos['x_L2com']
        halost['v_L2com'] = halos['v_L2com']
        halost['sigman_L2com'] = halos['sigman_L2com']
        halost['sigman_eigenvecsMin_L2com'] = halos['sigman_eigenvecsMin_L2com']
        halost['sigman_eigenvecsMid_L2com'] = halos['sigman_eigenvecsMid_L2com']
        halost['sigman_eigenvecsMaj_L2com'] = halos['sigman_eigenvecsMaj_L2com']
        halost['N'] = halos['N']
        
        # matching LRG density and n(z)
        halost.sort('N')   # to get the largest halos
        big_halos = halost[-int(n_in_slab):]  # matching LRG density and n(z)

        all_info.append(big_halos)    
    final_sample = vstack(all_info)
    
    # removing 20Mpc from either side of x direction
    x_min = np.min(final_sample['x_L2com'][:,0])+20; x_max = np.max(final_sample['x_L2com'][:,0])-20
    final_sample = final_sample[((final_sample['x_L2com'][:,0]>x_min)&(final_sample['x_L2com'][:,0]<x_max))]
    return final_sample


for s in halo_sims_paths[::-1]:
    print('Working on', s)
    halo_paths = np.sort(glob.glob(s+'/halos/z0.800/halo_info/halo_info_0*.asdf'))
    sample = get_halo_info(halo_paths)
    print('N: ', len(sample))
    sample.write(save_directory+'halo_info_ph'+s[-3:]+'_z0.80.fits', overwrite=True)