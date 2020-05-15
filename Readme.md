<!-- dx-header -->
# expected_depth (DNAnexus Platform App)

## What does this app do?
Calculates the expected read depth per region per sample for a sequencing run.

## What are typical use cases for this app?
This app should be executed as part of a DNAnexus workflow for multiple samples.

## What data are required for this app to run?
This app requires all the region coverage files with their indexes and all the flagstat files for the whole run.
All region coverage files must contain the same set of regions

## What does this app output?
It outputs one file with the expected depth for each genomic region covered by the capture. The header contains the script name and all the samples used to calculate this. The rest of the file includes the genomic region (chrom, start, end), expected depth for that region and standard deviation. 

i.e. 

|                                                       | 
|-------------------------------------------------------| 
| 1       982953  983067  216.992222222   53.4980118248 | 
| 1       983156  983275  391.570555556   43.8037843078 | 
| 1       983392  983745  433.371111111   50.4381822929 | 
| 1       984247  984439  347.611111111   43.307407006  | 
| 1       984616  984831  432.396388889   55.1304850671 | 
| 1       984946  985175  286.7725        37.7168462835 | 
| 1       985283  985417  562.738611111   62.9361017065 | 
| 1       985613  985709  429.916111111   42.2017597655 | 
| 1       985807  985971  152.649166667   19.8170263554 | 
| 1       986106  986217  242.583055556   30.5412092397 | 
| 1       986633  986749  596.584166667   62.8941426015 | 
| 1       986833  987025  345.305555556   44.4336752246 | 
| 1       987108  987195  408.429722222   50.7999240993 | 

Output is normalised to represent a sample containing 100M usable reads

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

#### This app was made by EMEE GLH
