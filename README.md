# Range_size_winners_losers
This repository contains data and code necessary to reproduce results shown in the paper:  
> Wu-Bing Xu, Shane A. Blowes, Viviana Brambilla, Cher F. Y. Chow, Ada Fontrodona-Eslava, Inês S. Martins, Daniel McGlinn, Faye Moyes, Alban Sagouis, Hideyasu Shimadzu, Roel van Klink, Anne E. Magurran, Nicholas J. Gotelli, Brian J. McGill, Maria Dornelas, Jonathan M. Chase. Regional occupancy increases for widespread species but decreases for narrowly distributed species in metacommunity time series. Nature Communications. (accepted)

**Contact:** wubing.xu@idiv.de or wbingxu@gmail.com (Wubing Xu)

## Data
This metacommunity time series analyzed were selected from four databases: BioTIME, RivFishTIME, InsectChange, and a previously unpublished database (Metacommunity Resurveys). The BioTIME data can be accessed on Zenodo (https://doi.org/10.5281/zenodo.2602708) or through the BioTIME website (http://biotime.standrews.ac.uk/); the RivFishTIME data can be accessed through the iDiv Biodiversity Portal (https://doi.org/10.25829/idiv.1873-10-4000); the InsectChange data can be accessed on KNB (https://doi.org/10.5063/F11V5C9V) or through the data paper (http://onlinelibrary.wiley.com/doi/10.1002/ecy.3354/suppinfo); the ‘Metacommunity Resurveys’ data was compiled using the R code available here (https://github.com/chase-lab/metacommunity_surveys/tree/version-2) and can be accessed through the iDiv Biodiversity Portal (https://doi.org/10.25829/idiv.3503-jevu6s).

In this repository, the selected metacommunity time series can be accessed:
```
  data/Combined_assemblages.RDATA
  data/Combined_assemblages_same_locations.RDATA
```

The first RDATA file contains filtered dataset based on grid-based approach, and the second contains only sites in same geographic coordinates sampled across years (see manuscript for more details).

## R Analysis Files
These scripts were used to prepare the data for analysis, calculate occupancy, fit models and produce figures.

To run first few R scripts (named with prefix1-4), raw assemblage data are required to download from links described above. The output is the selected metacommunity time series and the estimates of species’ geographic range sizes (the file `Species_rangesize.RDATA` in the folder data).

Please note that some of the code in this repository was written to run on a HPC cluster.
