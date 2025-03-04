## R CMD check results

* Fixes to the notes:
1. "#VignetteBuilder" line was removed
2. Date field was updated
3. time taken for example >10 sec seems to be an internal wind-dev issue. I checked it manually with win-dev release and all examples execute in less than 10 seconds

This is the third update of the package. Several changes have been made to two major functions and new updates have been added. They're as follows:

1. parallelization enabled with parallel package
2. dupValidate function revised
3. per site Fis added to deviant detection
4. vstPermutation function added
5. maf modified to remove multi-allelic sites
5. FIT correction added



On check:  
There were no ERRORs, WARNINGs, or NOTES

## Downstream dependencies
There are no downstream dependencies included in the package so far.
  
