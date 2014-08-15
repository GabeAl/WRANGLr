WRANGLr
=======

Within-Rank Analysis of Novel Groups via Likeness in R
0.1 pre alpha


INTRODUCTION.
Use this package to produce a Newick-formatted tree of the relationships between OTUs (features) according to their mutual correlation with a binary outcome variable. A branch containing closely grouped OTUs may represent a single effective feature whose members, in concert, better correlate to the outcome variable. 

In practice, this tends to group related OTUs together by virtue of heightened inherent anti-correlation with each other's abundances but strong combined correlation with the outcome variable. This may have some biological utility in assigning relevance to groups of species found in frequent co-occurrance. It may also serve as a form of feature selection for machine-learning algorithms, as preliminary tests have indicated a small but consistent increase in the performance of Randomforest over using naive or non-existant feature selection.

USAGE.
Inputs are tab-delimited files in the format of QIIME's text files (without comments). For a full description of this format, consult QIIME references. 
Otherwise, "map" is a tab-delimited text file containing metadata elements in columns and named individual samples in rows. "Taxa" is a tab-delimited text file containing OTU/feature counts in columns and named individual samples in rows. They do not have to match in order or number of elements, but it is assumed that the mapping file contains a subset of the samples included in taxa (with matching names). 

Alternatively, the resulting parsed data structures can be fed into wranglR directly. Simply specify a NULL for the mapPath, a matrix for "taxaPath", and a name-order matching binary factor for "GroupCol".

wranglR = function(mapPath, taxaPath, GroupCol = "Group", taxThres = 0, randSubset = 0, dat=NA,
A=T, hardCut = F, shift = 1);

- mapPath is either a string path to a file containing the map, or NULL. 
- taxaPath is either a string path to a taxa file or a parsed, sorted matrix of OTUs (taxa) vs counts. 
- taxThres is the minimum average number of taxa that must be present across samples to be included. Specifying a cutoff of 0.00001 (fractional dataset) may significantly speed up computation time without significantly altering results.
- randSubset takes a random subset (fraction) of the data to operate on. It may speed up computation and is useful for comparing the outcome of the model on different subsets of input data. It may eventually have application for cross-validating a machine learning model (not yet integrated).
- dat takes a pre-computed setup object obtained from calling "setup" function with the same parameters as wranglR (without the dat parameter). It significantly speeds up computation of subsequent runs over the same dataset.
- A is a "method" parameter that alters the calculation type for correlation analysis. The default naive approach uses the variance on sums, while the alternative approach uses a granular variance calculation prior to assigning distance scores. A=T is more intuitive but produces more spurious results on test and real datasets, so it is recommended to set this to A=F, which produces results in-line with the non-tree based variants of wranglR.
- hardCut assigns all pairs of taxa for which distance is greater than 1 to 1. Effectively, this cap removes the algorithm's ability to distinguish between pairs that are "very distant" and "extremely distant." In practice this can result in cleaner trees without outlying distances. This is closer to an "additive tree" model where distances are capped between 0 and 1 (or 0 and shift). 
- shift specifies the amount to which highly distant pairs are offset. Values from 1 (recommended) to 5 are possible. Larger shift offsets result in a more tip-compacted tree preferentially excluding pairs known to reduce one another's correlation with the outcome variable when combined.

cutTree(XT, depth, useDist = F)
Not yet implemented. Takes input tree and cuts according to specified depth, returning a list of features and their associated P-value and Bonferroni corrected p-value. Its utility overlaps significantly with wranglR_clusters but provides less fine-grained control since it has to operate within an existing tree. It is also expected to overlap in functionality with existing, mature tree analysis software. 

USAGE EXAMPLES.
One example of the usage (single-run) is as follows: 
  wrangled = wranglr("map.txt","taxa.txt",isDiseased, A=TRUE);

For multiple runs on the same dataset, or to test different parameters (besides map, taxa, and outcome variable), considering saving the setup data object to reduce time complexity:
  dat = setup("map.txt","taxa.txt",isDiseased, A=FALSE);
  wrangled = wranglr("map.txt","taxa.txt",isDiseased, A=FALSE, dat=dat);
  write(wrangled$tree,"output.tre"); # for visualization or cutting/flow analysis in dedicated programs

The newick formatted string is stored in wranglR's output $tree. wranglR also provides other debugging objects as output. 
