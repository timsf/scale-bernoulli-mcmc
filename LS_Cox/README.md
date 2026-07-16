# DCBF-LSCP
Ox code for DCBF MCMC for level-set Cox process
Manuscript: Scalable Bernoulli Factory MCMC for Intractable Marginalised Posteriors
Authors: Timothée Stumpf-Fétizon and Flávio B. Gonçalves


The algorithms are coded in Ox

Please, download Ox Console (including OxEdit) from the link [https://www.doornik.com/download.html](https://www.doornik.com/download/oxmetrics9/Ox_Console/).
Ox is free for academic use.

The main files are: mcmc_data*.ox, for *=1, 2, 3 and 4. These are the files to be run and contain all the required specifications to be set.

The first row of the data file should contain its number of rows and columns (Ox format -- see included data files).
The first column contains the x coordinates and the second one the y coordinates.

The covariance function is set in file cov_func.ox - default is the power exponential with exponent 1.95.
