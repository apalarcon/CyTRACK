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
1 - The example run is for the case of Hurricane Irma in the North Atlantic basin from 28 August to 12 September 2017.

2 - The input data is from the ERA5 reanalysis. CyTRACK will automatically download the required ERA5 input data. Therefore, be sure you have installed and correctly configured the python CDS API (cdsapi) for data downloading.
