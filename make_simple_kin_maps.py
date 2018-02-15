import numpy as np
import pyfits
import yt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import PchipInterpolator
from scipy.interpolate import interp1d
import astropy
from astropy.cosmology import WMAP9 as cosmo
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
plt.ioff()
import astropy
import rebin
from astropy import units as u
import numpy as np
import math
from astropy.io import fits
import cPickle
import scipy
import scipy.ndimage
import scipy.stats as ss
import scipy as sp
from scipy.optimize import curve_fit
import pyfits
import os
import math
import string
import astropy
import rebin
from astropy import units as u
import sys
import struct
import matplotlib
import matplotlib.pyplot as pyplot
import astropy
import make_color_image
from matplotlib.backends.backend_pdf import PdfPages
from astropy import cosmology
from astropy.convolution import Gaussian2DKernel, Gaussian1DKernel, convolve_fft, convolve
from astropy.cosmology import WMAP9 as cosmo
import time
import warnings
import bottleneck
from joblib import Parallel, delayed

plt.close('all')


def write_fits(fits_filename, gal, kmaps):
    print '\tGenerating fits for %s...'%fits_filename
    master_hdulist = []
    prihdr = fits.Header()
    prihdr['COMMENT'] = "Storing kinematic maps."
    prihdr['simname'] = gal

    prihdu = fits.PrimaryHDU(header=prihdr)    
    master_hdulist.append(prihdu)

    colhdr = fits.Header()

    for i in arange(len(kmaps)):
        kmap, name = kmaps[i]

        master_hdulist.append(fits.ImageHDU(data = kmap, header = colhdr, name = name))

    print '\tSaving to ' + fits_filename
    thdulist = fits.HDUList(master_hdulist)
    thdulist.writeto(fits_filename, clobber = True)

    return



def gauss(x, *p):
    A, mu, sigma= p
    return A*numpy.exp(-(x-mu)**2/(2.*sigma**2))

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (average, math.sqrt(variance))


def rebin_function(im, new_shape):

    im_new=zeros((new_x, new_y))

    im_expand


    for i in arange(new_x):
        good_old_x=linspace(shape(im)[0]/new_x*i, shape(im)[0]/new_x*(i+1), shape(im)[0]/new_x)

        for j in arange(new_y):


            im_new[new_x,new_y]=rebin_function()


    return im_new

def build_vel_maps(min_x, max_x, im_size, fact, vel_los, Lbol, weights_all, y1_1, y2_1):

    # Create image parameters
    bin_size=(abs(max_x-min_x))/(im_size)
    # Create empty velocity and dispersion maps
    vel_map=zeros((im_size,im_size))
    disp_map=zeros((im_size,im_size))
    vel_map_blurred=zeros((im_size,im_size))
    disp_map_blurred=zeros((im_size,im_size))
    L_map = zeros((im_size, im_size))

    L_map[:,:]=nan
    vel_map[:,:]=nan
    disp_map[:,:]=nan
    vel_map_blurred[:,:]=nan
    disp_map_blurred[:,:]=nan

    for i in arange(im_size):
        for j in arange(im_size):
            bin_x = [min_x+i*bin_size,min_x+(i+1)*bin_size]
            bin_y = [min_x+j*bin_size,min_x+(j+1)*bin_size]
            good=where((bin_x[0] < y1_1) & (y1_1 < bin_x[1]) & (bin_y[0] < y2_1) & (y2_1 < bin_y[1]))[0]
            Lbol_i=sum(Lbol[good])
            L_map[i,j]=Lbol_i
            if len(good) > 5:
                vel_map[i,j], disp_map[i,j] = weighted_avg_and_std(vel_los[good], weights_all[good]) #weight by Lbol
    return vel_map, disp_map, L_map





class camera:
    ###This class has been adopted from Greg Snyder

    # The camera class performs the basic operations.
    # It takes as input 10 parameters from the CAMERAX-PARAMETERS HDUs created by Sunrise
    # The position and FOV units are in KPC
    # It returns an object containing these data plus methods for converting generic
    #   x,y,z coordinates from the simulation frame (in Physical kpc!!) into a camera-based coordinate system.
    # The camera coordinates are defined such that the axis ranges are [-1,1].
    #     The return coordinates can be modified to use a pixel-based grid instead, but this more generic function can be used for both the 
    #     CANDELized and perfect images (i.e., on the same axis extent, in matplotlib terms)
    # There is one remaining uncertainty -- the sense of the rows and columns in the stored images (and how they are read into a given language).
    #     In Python:pyfits/astropy, given simulation coordinates, the returned pixel values correspond to the location on the image given the following assumptions:
    #     The "imshow" command was run with origin='lower' and extent=(-1,1,-1,1)
    #     The images returned by pyfits from the broadband.fits or _candelized_noise.fits must be **TRANSPOSED** first
    #     Presumably there are other iterations of these two settings (or other ways to manipulate the images) that will be satisfactory (or better, if the reason they work is known).
    #      -Greg Snyder, 8/21/2014

    def __init__(self,x,y,z,dirx,diry,dirz,upx,upy,upz,fov):
        self.x=x
        self.y=y
        self.z=z
        self.dirx=dirx  #x3 unit vector w/ ref to lab frame
        self.diry=diry
        self.dirz=dirz
        self.upx=upx  #x2 unit vector
        self.upy=upy
        self.upz=upz
        self.fov = fov
        #These vectors are defined following the convention at http://en.wikipedia.org/wiki/Pinhole_camera_model
        self.x3vector = np.asarray([self.dirx,self.diry,self.dirz])
        self.x2vector = np.asarray([self.upx,self.upy,self.upz])
        self.x1vector = np.cross(self.x2vector,self.x3vector)
        self.x1vector = self.x1vector/np.linalg.norm(self.x1vector)

        # This is the heart of the matter.  The idea is to express the object's coordinates in the frame of the camera model defined in __init__.
        # Let the object's position expressed in the original frame be A, and the unit vectors i1, i2, i3 be those along the simulation axes.
        # Let the object's position defined in the camera's reference (without shifting the origin yet) be A'.
        # Then, with linear operators, A' = M A, where M is constructed by taking dot products of the camera's unit vectors i' with the original unit vectors i.
        # When the original frame is standard cartesian coords, this equation reduces to the algebra below.
    def express_in_camcoords(self,x,y,z):
        new_x = x*self.x1vector[0] + y*self.x1vector[1] + z*self.x1vector[2]
        new_y = x*self.x2vector[0] + y*self.x2vector[1] + z*self.x2vector[2]
        new_z = x*self.x3vector[0] + y*self.x3vector[1] + z*self.x3vector[2]
        return np.asarray([new_x,new_y,new_z])

    # vel is the velocity vector in the oringinal frame
    def xyzvel_to_los(self, vel):
        camvec=self.express_in_camcoords(vel[:,0], vel[:,1], vel[:,2])
        return camvec[2]

    #Wrapper that reconstructs the Sunrise pinhole camera model, expresses the object's position in the camera frame, and computes its position in the image plane.
    def xyz_to_pixelvals(self,x,y,z):
        camdist = (self.x**2 + self.y**2 + self.z**2)**0.5
        camvec = self.express_in_camcoords(x,y,z)
        #define focal length such that image values span from -1 to 1.
        f = camdist/(0.5*self.fov)
        #See guidance at http://en.wikipedia.org/wiki/Pinhole_camera_model
        y1 = (f/camdist)*camvec[0]*1.0
        y2 = (f/camdist)*camvec[1]*1.0

        return y1,y2


def make_maps(mass, x_pos, y_pos, z_pos, vel_x, vel_y, vel_z, im_size, max_x, min_x, camobj_1, star_age = None, star_age_min = 0., star_age_max = 2.e7, gas_temp = None):

    rad_pos = sqrt(x_pos**2.+y_pos**2.+z_pos**2.)
    if gas_temp != None:
        good = where((gas_temp.value[()] < 1.e4) & (rad_pos < 100))[0]
    if star_age != None:
        good = where((star_age > star_age_min) & (star_age < star_age_max) & (rad_pos < 100))[0]
    vel_x = vel_x[good].value[()]
    vel_y = vel_y[good].value[()]
    vel_z = vel_z[good].value[()]

    vel = array([vel_x, vel_y, vel_z]).transpose()
    mass = mass[good]
    x_pos, y_pos, z_pos = x_pos[good], y_pos[good], z_pos[good]
    #This shows how the object works on a single entry of x,y,z (can be anything).
    y1_1,y2_1       = camobj_1.xyz_to_pixelvals(x_pos,y_pos,z_pos)
    # Calculate the line of sight velocity from the given camera position and the velocity vector
    # of the particle
    print nanmax(vel)
    vel_los         = camobj_1.xyzvel_to_los(vel)
    print nanmax(vel_los)
    kpcyr=u.kpc/u.yr
    kms=u.km/u.s
    weights_all=mass.value[()]/max(mass.value[()])

    #Want to build the reconstructed cube from the highest resolution

    fact = 1.
    vel_map, disp_map, L_map = build_vel_maps (min_x = min_x, max_x = max_x, im_size = im_size/fact, fact = fact,
                                                                      vel_los = vel_los, Lbol = mass, weights_all = weights_all, 
                                                                      y1_1 = y1_1, y2_1 = y2_1)
    return vel_los, vel_map, disp_map



snaps = [['/Volumes/wd/vela_v2/tracers/data/10MpcBox_csf512_a0.370.d', 'VELA_v2_20']]


define read_data(snap_file, gal):
    ds = yt.load(snap_file)
    dd = ds.all_data()
    nir_cat = np.loadtxt('/Volumes/wd/art_hydro/catalogs/GEN3/'+gal+'/galaxy_catalogue/Nir_simplified_disc_cat.txt', skiprows = 1, dtype='str')
    nir_disc_cat = np.loadtxt('/Volumes/wd/art_hydro/catalogs/GEN3/'+gal+'/galaxy_catalogue/Nir_disc_cat.txt', skiprows = 1, dtype='str')
    scale = round(1.0/(ds.current_redshift+1.0),3)
    dd = ds.all_data()
    good = where(nir_cat[:,0] == 'a%.3f'%scale)[0][0]
    nir_cat = nir_cat[good]
    nir_disc_cat = nir_disc_cat[good]
    L_disk_s = nir_cat[8:11].astype('float')
    id_cen_star = nir_cat[1].astype('int')
    cen_star_offset  = nir_cat[2:5].astype('float')
    cen_star_voffset = nir_cat[5:8].astype('float')
    cold_cen  = nir_disc_cat[1:4].astype('float')
    L_disk = nir_disc_cat[7:10].astype('float')

    #Load simulation data

    print 'Loading gas velocity...'
    gas_vx = dd['gas', 'velocity_x']
    gas_vy = dd['gas', 'velocity_y']
    gas_vz = dd['gas', 'velocity_z']

    print 'Loading gas temperature...'
    gas_temp = dd['gas', 'temperature']

    print 'Loading gas density...'
    gas_dens = dd['gas', 'density']
    
    print 'Loading gas cell mass...'
    gas_mass = dd['gas', 'cell_mass']

    print 'Loading gas cell position...'
    gas_x = dd['index', 'x'].in_units('kpc')
    gas_y = dd['index', 'y'].in_units('kpc')
    gas_z = dd['index', 'z'].in_units('kpc')

    print 'Loading star positions...'
    star_x = dd['stars', 'particle_position_x'].in_units('kpc')
    star_y = dd['stars', 'particle_position_y'].in_units('kpc')
    star_z = dd['stars', 'particle_position_z'].in_units('kpc')
    star_mass = dd['stars', 'particle_mass'].in_units('Msun')
    star_creation_time = dd['stars', 'particle_creation_time'].in_units('yr')
    star_age_all = ds.arr(cosmo.age(ds.current_redshift).value, 'Gyr').in_units('yr') - star_creation_time
    print 'Loading star velocities...'
    star_vx = dd['stars', 'particle_velocity_x'].in_units('km/s')
    star_vy = dd['stars', 'particle_velocity_y'].in_units('km/s')
    star_vz = dd['stars', 'particle_velocity_z'].in_units('km/s')


    #Recenter the gas cells to galaxy center

    cen_x       = star_x[id_cen_star-1] - ds.arr(cen_star_offset[0], 'kpc')
    cen_y       = star_y[id_cen_star-1] - ds.arr(cen_star_offset[1], 'kpc')
    cen_z       = star_z[id_cen_star-1] - ds.arr(cen_star_offset[2], 'kpc')
    cen_vx      = star_vx[id_cen_star-1] - ds.arr(cen_star_voffset[0], 'km/s')
    cen_vy      = star_vy[id_cen_star-1] - ds.arr(cen_star_voffset[1], 'km/s') 
    cen_vz      = star_vz[id_cen_star-1] -  ds.arr(cen_star_voffset[2], 'km/s')

    gas_x_cen   = gas_x - cen_x
    gas_y_cen   = gas_y - cen_y
    gas_z_cen   = gas_z - cen_z
    gas_pos     = array([gas_x_cen, gas_y_cen, gas_z_cen])
    gas_pos_mag = sqrt(gas_x_cen**2. + gas_y_cen**2. + gas_z_cen**2.)

    gas_vx_cen  = (gas_vx - cen_vx).in_units('km/s')
    gas_vy_cen  = (gas_vy - cen_vy).in_units('km/s')
    gas_vz_cen  = (gas_vz - cen_vz).in_units('km/s')
    
    star_vx_cen = (star_vx - cen_vx).in_units('km/s')
    star_vy_cen = (star_vy - cen_vy).in_units('km/s')
    star_vz_cen = (star_vz - cen_vz).in_units('km/s')

    star_x_cen  =   star_x - cen_x 
    star_y_cen  =   star_y - cen_y
    star_z_cen  =   star_z - cen_z







define make_simple_kmaps(scale, simname, scale, cam_n):
    nir_cat_name = simname[0:-2]+'_v2_'+simname[-2:]
    nir_cat = np.loadtxt('/nobackupp2/rcsimons/catalogs/nir_catalogs/GEN3/'+nir_cat_name+'/galaxy_catalogue/Nir_simplified_disc_cat.txt', skiprows = 1, dtype='str')
    nir_disc_cat = np.loadtxt('/nobackupp2/rcsimons/catalogs/nir_catalogs/GEN3/'+nir_cat_name+'/galaxy_catalogue/Nir_disc_cat.txt', skiprows = 1, dtype='str')
    nir_mstar_cat = np.loadtxt('/nobackupp2/rcsimons/catalogs/nir_catalogs/GEN3/'+nir_cat_name+'/galaxy_catalogue/Mstar.txt')

    print '\n\n\n\t\t Running on (%s, %.3f, %i)'%(gal, scale, cam_n)
    camera = mcrx_data['CAMERA%i'%(cam_n)]   
    camera_params =  mcrx_data['CAMERA%i-PARAMETERS'%(cam_n)].header
    if True:
        camera.header['z'] = (round(1./scale - 1,4), 'redshift (w/ %s)'%cosmo.name)
        camera.header['ascale'] = (scale, 'scale factor')
        camera.header['camera'] = (cam_n, 'camera')


        camera.header['cameradist'] = (camera_params['cameradist'], '[kpc] Distance from origin to camera')
        camera.header['theta'] = (camera_params['theta'], '[rad] Angular coordinate theta of camera positi')
        camera.header['phi']   = (camera_params['phi'], '[rad] Angular coordinate phi of camera position')

        camera.header['CAMPOSX'] = (camera_params['CAMPOSX'], '[kpc] X position of camera        ')
        camera.header['CAMPOSY'] = (camera_params['CAMPOSY'], '[kpc] Y position of camera        ')
        camera.header['CAMPOSZ'] = (camera_params['CAMPOSZ'], '[kpc] Z position of camera        ')
        camera.header['CAMDIRX'] = (camera_params['CAMDIRX'], '[kpc] X comp of camera viewing dir')
        camera.header['CAMDIRY'] = (camera_params['CAMDIRY'], '[kpc] Y comp of camera viewing dir')
        camera.header['CAMDIRZ'] = (camera_params['CAMDIRZ'], '[kpc] Z comp of camera viewing dir')
        camera.header['CAMUPX']  = (camera_params['CAMUPX'] , '[kpc] X comp of camera Y axis dir ')
        camera.header['CAMUPY']  = (camera_params['CAMUPY'] , '[kpc] Y comp of camera Y axis dir ')
        camera.header['CAMUPZ']  = (camera_params['CAMUPZ'] , '[kpc] Z comp of camera Y axis dir ')
        camera.header['linear_fov'] = (camera_params['linear_fov'], '[kpc] Linear field of view in y-dir at origin') 



    if True:
        camposx = upx*dist_to_cam
        camposy = upy*dist_to_cam
        camposz = upz*dist_to_cam
        camdirx = -camposx/dist_to_cam
        camdiry = -camposy/dist_to_cam
        camdirz = -camposz/dist_to_cam
        camupx = 1
        camupy = 0
        camupz = 0
        fov_kpc = 50
        camobj_1 = camera(camposx,camposy,camposz,camdirx,camdiry,camdirz,camupx,camupy,camupz,fov_kpc)
        global pixscale_kpc, pixscale_arcsec_z
        global z, surf_bright_lim, age_limit, sig_kernel, cube_vel_per_bin, rebin_cut, camera_n

        rebin_cut = 1.0
        cube_vel_per_bin= 1.


        redshift = ds.current_redshift
        pixscale_kpc        = camobj_1.fov*(abs(max_x-min_x))/(2*im_size) #kpc/pixel
        pixscale_arcsec_z   = pixscale_kpc*cosmo.arcsec_per_kpc_proper(redshift).value  #this is the arcsec/pixel scale at redshift 2


    if True:
        print 'running for gas...'
        vel_los_gas, vel_map_gas, disp_map_gas = make_maps(gas_mass, gas_x_cen.value[()],  gas_y_cen.value[()], gas_z_cen.value[()], 
                                                           gas_vx_cen, gas_vy_cen, gas_vz_cen, 
                                                           im_size, max_x, min_x, camobj_1, gas_temp)
        print 'running for young...'

        vel_los_yng, vel_map_yng, disp_map_yng = make_maps(star_mass, star_x_cen.value[()],  star_y_cen.value[()], star_z_cen.value[()], 
                                                           star_vx_cen, star_vy_cen, star_vz_cen, im_size, max_x, min_x, camobj_1, 
                                                           star_age = star_age_all,  star_age_min = 0., star_age_max = 1.e8)
        print 'running for old...'

        vel_los_old, vel_map_old, disp_map_old = make_maps(star_mass, star_x_cen.value[()],  star_y_cen.value[()], star_z_cen.value[()], 
                                                           star_vx_cen, star_vy_cen, star_vz_cen, im_size, max_x, min_x, camobj_1, 
                                                           star_age = star_age_all, star_age_min = 1.e9, star_age_max = 1.e11)

    return vel_map_gas, disp_map_gas, vel_map_yng, disp_map_yng, vel_map_old, disp_map_old

    kmaps = [[vel_map_gas, 'gas'], [vel_map_yng, 'young'], [vel_map_old, 'old']]
    write_fits('../fits/%s.fits'%gal_name, gal_name, kmaps)


if __name__ == "__main__":

    args = parse()
    import yt

    if args['snap_files'] is not None: snaps = [args['snap_files']]
    else: snaps = np.asarray(glob.glob("*.d"))
        
    print "Generating Sunrise Input for: ", snaps

    abssnap = os.path.abspath(snaps[0])
    assert os.path.lexists(abssnap)

    dirname = os.path.dirname(abssnap)
    simname = os.path.basename(dirname) #assumes directory name for simulation name



    rad_min = 0
    rad_max = 50.
    gals = ['VELA%.2i'%(i+1) for i in np.arange(35)]


    Parallel(n_jobs = -1, backend = 'threading')(delayed(make_simple_kmaps)(scale, simname) for scale in scales)















