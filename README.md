# Range_size_winners_losers
This repository contains data and code necessary to reproduce results shown in the manuscript "Rich get richer: widespread species increase while narrow-ranged species decrease in biodiversity time series".

**Contacts:** 
  Wubing Xu – wubing.xu@idiv.de or wbingxu@gmail.com

## Data
This metacommunity time series analyzed were selected from four databases: BioTIME (19), RivFishTIME (20), InsectChange (21), and a previously unpublished database (Metacommunity Resurvey). The BioTIME data can be accessed on Zenodo (https://doi.org/10.5281/zenodo.2602708) or through the BioTIME website (http://biotime.standrews.ac.uk/); the RivFishTIME data can be accessed through the iDiv Biodiversity Portal: https://doi.org/10.25829/idiv.1873-10-4000; the InsectChange data can be accessed at http://onlinelibrary.wiley.com/doi/10.1002/ecy.3354/suppinfo; the ‘Metacommunity Resurvey’ data can be accessed through the iDiv Biodiversity Portal: https://doi.org/10.25829/idiv.3503-jevu6s (will be activated soon).
In this repository, the selected metacommunity time series can be accessed:

'''
Combined_assemblages_20220208.RDATA

Combined_assemblages_same_locations_20220208.RDATA
'''

The first RDATA file contains filtered dataset based on grid-based approach, and the second contains only sites in same geographic coordinates sampled across years (see manuscript for more details).

## R Analysis Files
These scripts were used to prepare the data for analysis, calculate occupancy, fit models and produce figures.

To run first few R scripts (named with prefix1-4), raw assemblage data are required to download from links described above. The output is the selected metacommunity time series and the estimates of species’ geographic range sizes (the file 'Species_rangesize.RDATA' in the folder data).

Please note that some of the code in this repository was written to run on a HPC cluster.
