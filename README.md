# Meta-Metabolome Ecology: A conceptual synthesis between ecosystem metabolomics and meta-community ecology
Unlike my other projects, this specific repository serves both to provide supplemental information for an
associated manuscript, but also will continue to serve as the location where all past and future meta-metabolomic,
ecological analyses are hosted. Specifically, those files that are not enclosed in the "Manuscript Specific Files" 
folder are broadly applicable to any FTICR report generated using Formularity (<a href="https://pubs.acs.org/doi/10.1021/acs.analchem.7b03318">Tolić et al., 2017</a>).

In order to perform analyses similar to those detailed in the main manuscript an order does need to be followed:
1) FTMS_Analysis.R - this script is based upon the ftmsRanalysis R package and will help process and organize the data
2) Transformation_Analysis_Dataset.R - this script will calculate all possible mass differences within a dataset and assign putative biochemical transformations based upon the provided database. This is necessary to generate 2 of the 3 dendrograms.
3) Generate_MCD.R - this script will create the molecular characteristics dendrogram, which will cluster metabolites based upon molecular properties (e.g., elemental composition, aromaticity index, etc.)
4) Generate_TD.R - this script will use information from step 2 to generate the transformation-based dendrogram, which estimates relationships based upon distances in a transformation network
5) Generate_TWCD.R - this script will use information from step 2 as well as moelcular information to generate the transformation-weighted characteristics dendrogram, which is a bit of a mixture between the MCD and TD.

From here you can perform other typcial (and atypical) ecological analyses, such as the included bNTI and Raup-Crick implementations.
