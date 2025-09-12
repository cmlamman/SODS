from nbodykit import transform
import fitsio, glob
import numpy as np
from numpy.lib.recfunctions import append_fields
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.table import Table
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=69.6, Om0=0.286, Ode0=0.714)


import sys
sys.path.append('/global/homes/c/clamman/IA/spec-IA')
from geometry_functions.ellipsoid_projection import *

def e_complex(a, b, theta):
    '''complex ellipticity, theta must be in rad'''
    abs_e = (1 - (b/a)) / (1 + (b/a))
    e1 = abs_e * np.cos(2*theta)
    e2 = abs_e * np.sin(2*theta)
    return e1, e2

halo_paths = np.sort(glob.glob('/pscratch/sd/c/clamman/abacus/all_halos_info/halo_info_ph*_z0.80.fits'))

def e_complex(a, b, theta):
    '''complex ellipticity, theta must be in rad'''
    abs_e = (1 - (b/a)) / (1 + (b/a))
    e1 = abs_e * np.cos(2*theta)
    e2 = abs_e * np.sin(2*theta)
    return e1, e2

#halo_paths = ['/pscratch/sd/c/clamman/abacus/all_halos_info/z0.8_04.fits']
for path in halo_paths:
    ph_n = path.split('ph')[1].split('_z')[0]
    ##ph_n = path.split('z0.8_')[1].split('.fits')[0]
    #if float(ph_n) > 10:
    #    continue
    print("working on ", ph_n)
    sample = Table.read(path)
    positions = sample['x_L2com'] / .7 # positions, must be in Mpc
    velocities = sample['v_L2com']
    dist_observer = cosmo.comoving_distance(0.8).value # in Mpc
    pos_observer = np.asarray([-dist_observer, 0, 0]) # puts Earth outside box and must be in Mpc
    print(type(positions))
    ras, decs, zes_noRSD = transform.CartesianToSky(positions, cosmo, observer=pos_observer) 
    ras, decs, zes_withRSD = transform.CartesianToSky(positions, cosmo, observer=pos_observer, velocity=velocities) 
    sample['Z_noRSD'] = zes_noRSD
    sample['Z_withRSD'] = zes_withRSD
    sample['RA'] = ras 
    sample['DEC'] = decs
    
    print('getting projected shape')
    a_ps=[]
    b_ps=[]
    theta_ps = []
    for h in sample:
        ai, bi, thetai = get_projected_shape_sphere(h)
        a_ps.append(ai); b_ps.append(bi); theta_ps.append(thetai)
    a_ps = np.asarray(a_ps); b_ps = np.asarray(b_ps); theta_ps = np.asarray(theta_ps)
    e1s, e2s = e_complex(a_ps, b_ps, theta_ps)
    sample['E1'] = e1s
    sample['E2'] = e2s
    
    print('saving')
    sample.write('/pscratch/sd/c/clamman/abacus/halos_withRSD_shapes/halo_selection_ph'+ph_n+'_z0.80.fits', overwrite=True)
    
    
    
    ## add randoms??