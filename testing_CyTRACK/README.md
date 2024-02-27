# Testing CyTRACK
## Steps

1 - You can check the reproducible CyTRACK capsule in Code Ocean

OR

1 - After install CyTRACK, manually run the testing example on your local PC. 

* Clone the CyTRACK repository.

 ```
https://github.com/apalarcon/CyTRACK.git
  ```

* Go to CyTRACK/testing_CyTRACK directory
```
cd  CyTRACK/testing_CyTRACK
```
* Run the testing_CyTRACK.sh script

```
sh testing_CyTRACK.sh
```
### Important Notes
1 - The example run is for the dates of the Hurricane Irma in the North Atlantic basin from 28 August to 12 September 2017. Note that several tropical cyclones occurred at the same time during this period. Therefore CyTRACK captured all of then.

2 - The input data is from the ERA5 reanalysis. CyTRACK will automatically download the required ERA5 input data. Therefore, be sure you have installed and correctly configured the python CDS API (cdsapi) for data downloading (see <a href="https://cds.climate.copernicus.eu/api-how-to" target="blank"> How to use the CDS API - Climate Data Store - Copernicus </a>).

## Testing Results
1 - If CyTRACK runs successfully, the CyTRACK_OUTPUTS/CyTRACK_output directory should be created. In this directory, you should find the following file ```CyTRACK_AL_2017082500-2017091500_ERA5_TC.dat```, containing the information on the identified tropical cyclones.

2 - To plot the cyclones tracks, run the ```plotting_test_CyTRACK_outputs.py``` script.
```
python plotting_test_CyTRACK_outputs.py
```
3 - Your should obtain the following map

