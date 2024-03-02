import pandas as pd
import numpy as np

from astropy.io import fits
from astropy import wcs
import astropy.units as u
import sys
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser(description="Process some integers.")

# Adding the arguments
parser.add_argument('-method', type=int, help='The method to use. 0 for making a density map from scratch, 1 for using an X-ray FITS file grid.')
parser.add_argument('-data_file', type=str, help='The data file in CSV format containing the RA and Dec of the galaxies as "ra" and "dec" respectively.')
parser.add_argument('-output_file', type=str, help='The output file in FITS format to save the density map.')
parser.add_argument('-step', type=float, help='The step size for the grid in degrees (only for method 0).')
parser.add_argument('-margin', type=float, help='The margin to add to the grid in percentage (only for method 0 and default 10 percentage)', default=10.0)
parser.add_argument('-xray_fits', type=str, help='The X-ray FITS file to use for creating the density map', default=None)
parser.add_argument('-wavelet_exp_map', type=bool, help='Whether to create exposure map for the wavelet filtering or not', default=False)
args = parser.parse_args()

# Accessing the arguments
method = args.method
data_file = args.data_file
output_file = args.output_file
step = args.step
margin = args.margin
xray_fits = args.xray_fits
wavelet_exp_map = args.wavelet_exp_map

if method==None:
    print('ERROR: Please provide the method to use. 0 for making a density map from scratch, 1 for using an X-ray FITS file grid.')
    sys.exit()

if data_file==None or output_file==None:
    print('ERROR: Please provide the data file and the output file name.')
    sys.exit()

df = pd.read_csv(data_file)
if 'ra' not in df.columns:
    print('ERROR: Please provide the RA and Dec columns in the data file as "ra" and "dec" respectively.')
    sys.exit()


if method==0:
    if step==None:
        print('ERROR: Please provide the step size.')
        sys.exit()
    if margin==None:
        print('Using default margin of 10%')
        margin = 10
    # Creating the grid
    ra_min, ra_max = df['ra'].min(), df['ra'].max() 
    dec_min, dec_max = df['dec'].min(), df['dec'].max()

    ra_ext, dec_ext = ra_max-ra_min, dec_max-dec_min

    grid_ra = np.arange(ra_min-ra_ext/margin, ra_max+ra_ext/margin, step) #Here we add a 10% margin to the grid 
    grid_dec = np.arange(dec_min-dec_ext/margin, dec_max+dec_ext/margin, step)

    ra_size, dec_size = grid_ra.shape[0], grid_dec.shape[0]
    
    w = wcs.WCS(naxis=2)
    # First we add the projection type
    w.wcs.ctype = ["RA---SIN", "DEC--SIN"]
    # Then the units of the coordinates
    w.wcs.cunit = ["deg", "deg"]
    # Now we set the center/reference pixel of the grid
    w.wcs.crpix = [ra_size/2, dec_size/2] 
    # This sets the value of the reference pixel in the same units as the coordinates
    w.wcs.crval = [grid_ra[int(ra_size/2)], grid_dec[int(dec_size/2)]]
    # This is the pixel scale in degrees. We take the difference between two adjacent pixels in the grid 
    w.wcs.cdelt = np.array([(grid_ra[int(ra_size/2)]-grid_ra[int(ra_size/2)+1]), (grid_dec[int(dec_size/2)+1]-grid_dec[int(dec_size/2)])])
    
    # To add the size of the grid to the header
    new_header = w.to_header().copy() # copy the header from the wcs object
    new_header['NAXIS']=2
    new_header['NAXIS1']=ra_size
    new_header['NAXIS2']=dec_size
    new_wcs = wcs.WCS(new_header) # create a new wcs object with the new header
    # The new wcs has the same information as the old one, but with the NAXIS updated
    output_map_shape = (new_header['NAXIS1'], new_header['NAXIS2'])
    print('The output density map will have dimension of:', output_map_shape)
    
    print('Creating the density map...')
    density = np.zeros(output_map_shape)
    # Loop over the galaxies and add them to the grid
    for i in tqdm(range(len(df))):
        # Convert the world coordinates (ra, dec) to pixel coordinates
        app_ra, app_dec = new_wcs.world_to_pixel_values(df['ra'][i], df['dec'][i])
        # Round the pixel coordinates to the nearest integer
        pix_ra, pix_dec = int(np.rint(app_ra)), int(np.rint(app_dec))
        
        # Increment the density value at the corresponding pixel in the density map
        density[pix_ra,pix_dec] += 1
    # Save the density map as a FITS file
    hdu = fits.PrimaryHDU(density.T, new_header)
    hdu.writeto(output_file, overwrite=True)
    print('Density map saved as', output_file)
    print('Done')

if method==1:
    if xray_fits==None:
        print('ERROR: Please provide the X-ray FITS file.')
        sys.exit()
    # Open the xray FITS file
    ref_image = fits.open(xray_fits)[0]
    ref_wcs = wcs.WCS(ref_image.header)
    
    output_map_shape = (ref_image.header['NAXIS1'], ref_image.header['NAXIS2'])
    print('The output density map will have dimension of:', output_map_shape)
    
    def get_box_coord(pix_ra, pix_dec, image_wcs, ext=1):
        """
        This function calculates the coordinates of a box around a pixel in an image.
        
        Parameters:
        pix_ra (float): Pixel coordinate in the x-axis.
        pix_dec (float): Pixel coordinate in the y-axis.
        image_wcs (astropy.wcs.WCS): WCS of the image.
        ext (int): Extent of the box in terms of pixels.
        
        Returns:
        tuple: A tuple containing the coordinates of the box (ra_prev, dec_prev, ra_next, dec_next).
        """
        ra_next = image_wcs.pixel_to_world(pix_ra+(ext/2), pix_dec).ra.value
        dec_next = image_wcs.pixel_to_world(pix_ra, pix_dec+(ext/2)).dec.value
        ra_prev = image_wcs.pixel_to_world(pix_ra-(ext/2), pix_dec).ra.value
        dec_prev = image_wcs.pixel_to_world(pix_ra, pix_dec-(ext/2)).dec.value
        
        # Sometimes the coordinates are in decrasing order, so we need to swap them
        if ra_next < ra_prev:
            ra_next, ra_prev = ra_prev, ra_next
        if dec_next < dec_prev:
            dec_next, dec_prev = dec_prev, dec_next
        return ra_prev, dec_prev, ra_next, dec_next
    print('Creating the density map...')
    density = np.zeros(output_map_shape)
    for i in tqdm(range(len(df))):
        app_ra, app_dec = ref_wcs.world_to_pixel_values(df['ra'][i], df['dec'][i])
        pix_ra, pix_dec=int(np.rint(app_ra)), int(np.rint(app_dec))
        ra_min, dec_min, ra_max, dec_max = get_box_coord(pix_ra, pix_dec, ref_wcs)
        mask = (df['ra']>=ra_min) & (df['ra']<=ra_max) & (df['dec']>=dec_min) & (df['dec']<=dec_max)
        density[pix_ra,pix_dec] = mask.sum()
        # If the mask is empty due to edge effects, we set the density to +1
        if mask.sum()==0:
            density[pix_ra,pix_dec] = +1
    
    # Save the density map as a FITS file
    hdu = fits.PrimaryHDU(density.T, ref_image.header)
    hdu.writeto(output_file, overwrite=True)
    print('Density map saved as', output_file)
    
    if wavelet_exp_map:
        exp = np.ones(output_map_shape)
        hdu_exp = fits.PrimaryHDU(exp.T, ref_image.header)
        hdu_exp.writeto(output_file.replace('.fits', '_exp.fits'), overwrite=True)
        print('Exposure map saved as', output_file.replace('.fits', '_exp.fits'))
    
    print('Done')
    