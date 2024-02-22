For simplicity, the scripts in this folder are organized according to the
corresponding major pipeline:

	* scripts_getinputs/ - collects the sripts to process the raw RBP
	  CisBP files downloaded form http://cisbp-rna.ccbr.utoronto.ca, fetch
	  the RNA sequences around splicesites that will be scanned, and
	  derive the clusters of representative motifs.
	* scripts_affimx/ - contains all the steps to run AffiMx to scan the
	  previously processed RNA sequences with the set of representative
	  motifs obtained from the motif clustering analysis.
	* scripts_models/ - has the scripts needed to fit the models of
	  differential splicing regulation for all representative RBP motifs
	  accross cancer types.

The `input/` and `output/` folders are organized according to the
corresponding program or script, which could come from either of the three
scripts folders. 
