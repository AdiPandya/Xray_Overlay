
#Install reproject package if not already installed
#pip install reproject
import sys
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from astropy.utils.data import get_pkg_data_filename
from reproject import reproject_interp

xray_image = sys.argv[1]
optical_data = sys.argv[2]
reproj_image = sys.argv[3]
reproj_fits = sys.argv[4]

# Open the xray FITS file
hdu_x = fits.open(get_pkg_data_filename(xray_image))[0]

# Extract the data and the WCS
data_x = hdu_x.data
header_x = hdu_x.header
data_x = np.nan_to_num(data_x) # Replace NaN values with 0

# Open the optical FITS file / wcs file
hdu_O = fits.open(get_pkg_data_filename(optical_data))[0]

# Extract the data and the WCS
data_O = hdu_O.data
header_O = hdu_O.header

if data_O == None:
    new_header = header_O.copy()
    new_header['NAXIS'] = 2 # Since it is a 2D image
    new_header['NAXIS1'], new_header['NAXIS2'] = new_header['IMAGEW'], new_header['IMAGEH'] # Set the dimensions same as the image height and width
else:
    new_header = header_O

print('Starting reprojection...')
array, footprint = reproject_interp(hdu_x, new_header)
print('Reprojection complete.')

ax1 = plt.subplot(1,2,1, projection=WCS(new_header))
ax1.imshow(array, origin='lower')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination')
ax1.set_title('Reprojected X-ray Image')

ax2 = plt.subplot(1,2,2, projection=WCS(new_header))
ax2.imshow(footprint, origin='lower')
ax2.coords['ra'].set_axislabel('Right Ascension')
ax2.coords['dec'].set_axislabel('Declination')
ax2.coords['dec'].set_axislabel_position('r')
ax2.coords['dec'].set_ticklabel_position('r')
ax2.set_title('X-ray image footprint')

print('Image saved as', reproj_image)

plt.savefig(reproj_image, dpi=300, bbox_inches='tight')

fits.writeto(reproj_fits, array, new_header, overwrite=True)
print('FITS file saved as', reproj_fits)
print('Done')
