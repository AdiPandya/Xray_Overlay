# Xray Overlay Guide

Guide for making X-ray overlay on the Optical Images. Read the instruction below to get started.

To overlay the X-ray images on the Optical images, basic requirements are:
1. X-ray data fits file with correct WCS (World Coordinate System) information (usually available with the data)
2. Optical data as fits file or an image.
3. Software: FITS Liberator (Download at: https://noirlab.edu/public/products/fitsliberator/)
4. Software: GIMP (Download at: https://www.gimp.org/downloads/)

The steps to overlay are as follows:

1. Astrometric Calibration of Optical Data:
    - If the optical data is in fits format, it should have the correct WCS information. If the optical data is in image format, one needs to use the [astrometry.net](https://nova.astrometry.net/upload) webpage to get the WCS information. Upload the image and once the astrometric calibration is done, download the fits file with the WCS information or just the WCS information as a fits file. For the purpose of this work, only the WCS information will be used.
2. Align the X-ray data with the Optical data:
    - Use the script `Align.py` to align the X-ray data with the Optical data. The script uses the WCS information of the Optical data to align the X-ray data. The script can takes 4 arguments:
        - The path to the X-ray data fits file.
        - The path to the Optical data fits file.
        - The path to reprojected X-ray image.
        - The path to the output reprojected fits file.
    - A tutorial on how to use the script is available in the `Align_tutorial.ipynb` notebook.
3. Using FITS Liberator to change the scale of the xray reprojected image:
    - Open the reprojected X-ray fits file in FITS Liberator.
    - Change the scaling of the image to power-law or arcsinh scaling to make the X-ray sources more visible.
    - Save the image as a .tiff file. 
    - More details on how to use FITS Liberator are available in the `Image_overlay_tutorial.ipynb` notebook.
4. Overlay the X-ray image on the Optical image using GIMP:
    - Open the Optical image in GIMP.
    - Open the X-ray image as a new layer.
    - Follow the steps in the `Image_overlay_tutorial.ipynb` to make the overlay.
    - Save the overlayed image as a .tiff/.png/.jpg file.

### Acknowledgement:
- The `Image_overlay_tutorial.ipynb` is created with the help of the guide available at: https://chandra.harvard.edu/photo/openFITS/
- The x-ray data used in the tutorial is from the Chandra X-ray Observatory. The data is created to NASA/CXC/Wesleyan Univ./R.Kilgard, et al
- The optical data used in the tutorial is from the Hubble Space Telescope. The data is created to NASA, ESA, S. Beckwith (STScI), and The Hubble Heritage Team (STScI/AURA)