Varvar
--

[![GitHub version](https://badge.fury.io/gh/bretonics%2Fvarvar.svg)](http://badge.fury.io/gh/bretonics%2Fvarvar)
[![Github Issues](http://githubbadges.herokuapp.com/bretonics/varvar/issues.svg)](https://github.com/bretonics/varvar/issues)
[![Pending Pull-Requests](http://githubbadges.herokuapp.com/bretonics/varvar/pulls.svg)](https://github.com/bretonics/varvar/pulls)
[![GitHub license](https://img.shields.io/badge/License-MIT-orange.svg)](https://bretonics.mit-license.org/)
![](https://reposs.herokuapp.com/?path=bretonics/varvar&color=lightgrey)

Useful utilities for various variant call manipulations from NGS assembly/alignment results.

## Lasergene - DNASTAR Software
Manipulate SNP table from Lasergene results

### snpExtract
Filter out SNP calls accordingly. Set **SNP %** limits and **depth** boundaries to filter SNPS calls you wish to call as true and change with `variantChanger`.

`snpExtract -snp [0.00 100.00] -depth [0 Inf]`

### variantChanger
Use the output from `snpExtract` as input to make SNP changes in original sequence.

`variantChanger -seq <sequenceFile> -snp <filteredSNPsFile> -out <outFile`


## License [![GitHub license](https://img.shields.io/badge/License-MIT-orange.svg)](https://bretonics.mit-license.org/)
