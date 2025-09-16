import numpy as np
import matplotlib.pyplot as plt
import glob, random, fitsio

from astropy.table import Table, join, vstack
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u
from astropy.cosmology import LambdaCDM, z_at_value
from astropy.io import fits

import pandas as pd

from functions.lightProfile_functions import *
from functions.ellipsoid_projection import *

def get_ab_selection(sample, save_path, ph_num='000', ap_rad = 1.5/2, N=100000, light_profile='SER'): 
    '''
    creates and saves file of halo 'mock' catalog which includes input columns plus:
    shape_r: angular half-light radius in arcminutes
    r_half_model: estimated half-light radius of model in arcseconds
    ap_sum: fraction of light within aperture
    e_los: projected shape along line of sight
    E1, E2: projected shapes on sky
    
    sample: table of halos containing at least sigman_L2com, sigman_eigenvecsMin_L2com, sigman_eigenvecsMid_L2com, sigman_eigenvecsMaj_L2com, Z
    ph_num: string, len of 3
    ap_rad: aperture radius in arcseconds, default is DESI
    N: number of points to sample light profile, default is 100000
    light_profile: 'SER' for Sersic, 'EXP' for exponential, 'DEV' for de Vaucouleurs
    '''
    
    rad_image_sq = ap_rad**2
    h_points = get_light_profile(light_profile)
    gaussian_displacement = np.random.normal(0, 1/2.355, N*2).reshape(2,N)
    
    ##### fiber mags ##### 
    
    # half light radii
    p_sizes = np.random.lognormal(mean=6.9, sigma=.7, size=len(sample)) # physical half-light radii
    sample['shape_r'] = (p_sizes * u.Mpc / LambdaCDM(H0=70.4 * u.km / u.s / u.Mpc, Om0=0.272, Ode0=0.73).\
             angular_diameter_distance(z=np.asarray(sample['Z']))).value
    
    k=0
    r_half_models = []
    ab_apSum = []
    ab_e1_LOS = []
    sky_e1 = []
    sky_e2 = []
    
    for ab in sample[k:]:
        if k%20000==0:
            print('working on ', k+1, ' / ', len(sample))

        # estimate model half-light from many orientations
        # this part is annoying, but necessary for good scale estimation
        ax = 1 / (ab['sigman_L2com']**2)  # axis lengths
        a3d = 1; b3d = ax[1]/ax[2]; c3d = ax[0]/ax[2]
        r_halfs = []
        for r_matrix in R.random(10):
            # scale, rotate, and project ellipsoid points along z axis
            M0 = shape_transformation_2D(a3d, b3d, c3d, r_matrix, h_points, p_axis='z')
            r_halfs.append(find_r_half(M0))
        r_half_model = np.mean(r_halfs) 
        r_half_models.append(r_half_model)
        

        # make light profile based on 3D shape
        # X is LOS in Abacus mock
        # scale, rotate, and project ellipsoid points along x axis
        scale_matrix = np.array([a3d*ab['sigman_eigenvecsMin_L2com'], b3d*ab['sigman_eigenvecsMid_L2com'], c3d*ab['sigman_eigenvecsMaj_L2com']]).transpose()
        M3D = np.matmul(scale_matrix, h_points)
        M0 = M3D[:2]  # retain just the y and z direction

        # find scale (arcesconds / grid unit) based on the galaxy's half light radius
        image_scale = (ab['shape_r']*0.7) / r_half_model
        M1 = M0 * image_scale

        # gaussian displacement matrix
        M = M1 + gaussian_displacement

        # estimate fraction within aperture
        a_sum = app_sum(M, rad_image_sq = rad_image_sq)
        ab_apSum.append(a_sum)
        
        
        ##### projected shapes ##### 
        
        # CALCULATING PROPERTIES OF 2D ELLIPSOID as viewed from transverse direction
        # want to project along y, to get theta relative to LOS (x), in direction of +z

        # eigen vectors
        evc = np.array([ab['sigman_eigenvecsMaj_L2com'], ab['sigman_eigenvecsMid_L2com'], ab['sigman_eigenvecsMin_L2com']])
        evl = ab['sigman_L2com']**2 # eigen values
        
        # for p_axis='y', theta is the angle relative to x, in the direciton of +z
        K = np.sum(evc[:,1][:,None] * (evc / evl[:,None]), axis=0)
        r = evc[:,0] - evc[:,1]*K[0]/K[1]
        s = evc[:,2] - evc[:,1]*K[2]/K[1]

        A = np.sum(r**2 / evl, axis=0)
        B = np.sum(2*r*s / evl, axis=0)
        C = np.sum(s**2 / evl, axis=0)
        # for p_axis='y', theta is the angle relative to x, in the direciton of +z
        theta = np.pi/2 + np.arctan2(B, A-C) / 2
        a_p = 1 / np.sqrt((((A+C)/2) + ((A-C)/(2*np.cos(2*theta)))))
        b_p = 1 / np.sqrt(A + C - (1/a_p**2))
        e1, e2 = e_complex(1, b_p/a_p, theta)
        ab_e1_LOS.append(e1)
    
    
        # PROPERTIES OF 2D ELLIPSOID ON PLANE OF SKY (y - z, curved because sky)
        
        # for p_axis='x', theta is the angle relative to z, in the direciton of +y
        K = np.sum(evc[:,0][:,None] * (evc / evl[:,None]), axis=0)
        r = evc[:,2] - evc[:,0]*K[2]/K[0]
        s = evc[:,1] - evc[:,0]*K[1]/K[0]

        A = np.sum(r**2 / evl, axis=0)
        B = np.sum(2*r*s / evl, axis=0)
        C = np.sum(s**2 / evl, axis=0)
        # for p_axis='x', theta is the angle relative to z, in the direciton of +y
        theta = np.pi/2 + np.arctan2(B, A-C) / 2
        a_p = 1 / np.sqrt((((A+C)/2) + ((A-C)/(2*np.cos(2*theta)))))
        b_p = 1 / np.sqrt(A + C - (1/a_p**2))
        e1, e2 = e_complex(1, b_p/a_p, theta)
        sky_e1.append(e1); sky_e2.append(e2)
        
        k+=1
        
    sample['r_half_model'] = r_half_models
    sample['ap_sum'] = ab_apSum
    sample['e_los'] = ab_e1_LOS
    sample['E1'] = sky_e1
    sample['E2'] = sky_e2
    
    print('saving')
    sample.write(save_path, overwrite=True)
    
catalog_paths = glob.glob('/pscratch/sd/c/clamman/abacus/halos_withRSD_shapes/halo_selection_ph*_z0.80.fits')

n_batches = 20
for cat in catalog_paths:
    ph_num = cat.split('selection_ph')[1][:3]
    save_directory = '/pscratch/sd/c/clamman/abacus/halo_selections/halo_selection_ph'+ph_num+'_z0.80_'
    print('Working on ph', ph_num)
    halo_sample = Table.read(cat)
    halo_sample['Z'] = halo_sample['Z_withRSD']
    halo_sample = halo_sample[((np.abs(halo_sample['RA']) < 20) & (np.abs(halo_sample['DEC']) < 20) & (halo_sample['Z'] < 1.4) & (halo_sample['Z'] > 0.6))]

    # run catalog in chunks to periodically save
    batch_size = len(halo_sample) // n_batches
    for i in range(n_batches):
        save_path = save_directory + str(i) + '.fits'
        # see if saved file already exists
        if glob.glob(save_path):
            print('batch', i, ' done')
            continue
        batch = halo_sample[i*batch_size:(i+1)*batch_size]
        get_ab_selection(batch, ph_num=ph_num, save_path=save_path)