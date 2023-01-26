# wcw
The replication code for the "When Ties are Possible: Weak Condorcet Winners and Arrovian Rationality" by Diss, Dougherty, and Heckelman.

This is the code for the simulations used in "When Ties are Possible: Weak Condorcet Winners and Arrovian Rationality" by Diss, Dougherty, and Heckelman.  The primary R file wcw.R references an input file (in the same directory) and a single Rcpp file iterate2.cpp in the subdirectory <\c_functions>.  Output is written to the subdirectory <\output> which should be created on your computer.  To replicate the data, download the file and directory structure and run wcw.R.

Several functions are designed for three alternatives cases and may not work with more alternatives.  Further notes, are written in the header of wcw.R and in various comments.
