# Density Maps

The notebook [Density_Maps_tutorial.ipynb](Density_Maps_tutorial.ipynb) contains the code to create the density maps from their celestial coordinates and save them as fits files with correct WCS information. The density maps can be used to create contours on the X-ray images.

Differnt methods are mentioned in the notebook to create the density maps. A scipt version of the notebook is available as [Density_Maps.py](Density_Maps.py).

The script takes several arguments to create the density maps. The arguments are: \
-method : The method to use. 0 for making a density map from scratch, 1 for using an X-ray FITS file grid. \
-data_file : The data file in CSV format containing the RA and Dec of the galaxies as "ra" and "dec" respectively. \
-output_file : The output file in FITS format to save the density map.

-step : The step size for the grid in degrees (only for method 0). \
-xray_fits : The X-ray FITS file to use for creating the density map (only for method 1)

Optional arguments: \
-margin : The margin to add to the grid in percentage (only for method 0 and default 10 percentage) \
-wavelet_exp_map : Whether to create exposure map for the wavelet filtering or not (only for method 1 and default False) \
-extent : The extent of the box in terms of pixels (only for method 1 and default 1 pixel)

The script for method 0 can be run as:
```bash
python Density_Maps.py -method 0 -data_file data.csv -output_file density_map.fits -step 0.1
```

For method 1, run the script as follows:
```bash
python Density_Maps.py -method 1 -data_file data.csv -output_file density_map.fits -xray_fits xray.fits
```

To see the help for the script, run:
```bash
python Density_Maps.py -h
```

There are methods to create the density maps. In method 0, the density map is created from scratch i.e. a new grid is created and wcs information is added to the fits file. Here you can set the step size of the grid in degrees and additional margin to add to the grid.

In method 1, the density map is created using the grid of an X-ray FITS file. The X-ray FITS file might have very fine grid and the density map created from it will majorly contain only 0s. This image can be passed through a wavelet filtering scipt to get a smooth density map. The wavelet filtering script also requires the exposure map of the X-ray FITS file. The exposure map can be created by passing the argument `-wavelet_exp_map True` in the script. (By default, the exposure map is not created.)

Using the method 1, extent of the box in terms of pixels is 1. This means that for each pixel, a box of 1 pixel is created and the density is calculated in that box. If larger extent is required, it can be set using the optional argument `-extent`.

In method 1, if the CSV file contains the coordinates that are not in the range of the X-ray FITS file, the density map will contain only 0s. If this happens, one can see how many galaxies are outside and inside the range of the X-ray FITS file once the script is run. 

## Acknowledgement
- The X-ray FITS file in the Data_Files directory is created from the XMM-Newton pointing observation of the galaxy cluster RXCJ0631.3-5610.
- The code for method 1 is written by referencing the work of Jakob Dietl and Caroline Mannes
