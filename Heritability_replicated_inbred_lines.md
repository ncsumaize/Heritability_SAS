
EXAMPLE SAS PROGRAM #1 FOR ESTIMATING HERITABILITY AND ITS STANDARD ERROR
Example SAS codes correspond to J.B. Holland, W.E. Nyquist, and C.T. Cervantes-Martinez, "Estimating and Interpreting Heritability in Plant Breeding: An Update," submitted to(2003) Plant Breeding Reviews.  

Heritability based on variance among random S0:2 lines families evaluated in an unbalanced, multi-environment trial.  Options are given for estimating heritability from one environment or multiple environments, and from randomized complete block designs or incomplete block (lattice) designs within each environment.

Data are taken from Hoi et al. (1999)   Crop Sci. 39:1055-1059.

Briefly: yield was measured on 95 random S0:2 lines in five environments in Iowa.

To perform the analysis you first need to download the data files:

[Original data on a plot basis for all lines (including checks) is available in an Excel spreadsheet file](OatYield.xls). 

I have developed four different example programs, corresponding to estimating heritability from different experimental designs:

1. [SAS code for estimating heritability from lines evaluated in an RCB design in a single environment](OneEnvironRCBDHeritability.sas).  [Output from this analysis](MultiEnvironLatticeHeritability.sas).

2. [SAS code for estimating heritability from lines evaluated in an incomplete block (lattice) design in a single environment](OneEnvironLatticeHeritability.sas).   [Output from this analysis](OneEnvironLatticeHeritOut.lst).

3. [SAS code for estimating heritability from lines evaluated in RCB designs in multiple environments](MultiEnvironRCBDHeritability.sas).  [Output from this analysis](MultiEnvironRCBDHeritOut.lst).

4. [SAS code for estimating heritability from lines evaluated in incomplete block (lattice) designs in multiple environments](MultiEnvironLatticeHeritability.sas).   [Output from this analysis](MultiEnvironLatticeHeritOut.lst).
