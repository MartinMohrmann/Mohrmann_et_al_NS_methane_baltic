# Mohrmann_et_al_NS_methane_baltic
Source codes for the paper manuscript "Nord Stream methane leaks spread across 14% of Baltic waters"

Steps to reproduce the work:
1. Download ocean model input files from https://data.marine.copernicus.eu/product/BALTICSEA_ANALYSISFORECAST_PHY_003_006/description
2. Download atmospheric forcing input files from https://cds.climate.copernicus.eu/cdsapp#!/software/app-era5-explorer?tab=overview
3. Adjust the Nordsstream_model_run.py to use one of the configuration files. The configuration files are inputparameters194.py (default run), and inputparameters[170,180.190] with different oxidation rates respectively.
4. run python Nordstream_model_run.py to generate output files in trajectory format.
5. run python create_histogram.py to generate a gridded output format from the trajectories generated in the previous step
6. Now the plots presented in the Ipython-notebooks can be reproduced.
