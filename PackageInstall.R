## Please make sure that your install of R and/or RStudio are up to date. (at least 3.2.2, preferably 3.3.*)

# Install and Load package from GitHub
# Run this file only once!


install.packages(c('rrcov','doparallel','TInPosition','rpca','reshape','ggplot2','colorspace','robustbase', 'SIBER', 'psych'),dependencies = T)


	## install from either of the options below; running both is not necessary.

## INSTALL OPTION (1)
	## If the following block of code does not work (i.e., a failure to install certain packages) then use the line of code below this block.
	install.packages("devtools")
	require(devtools)
	install_github("derekbeaton/outlieRs", subdir="Package")



## INSTALL OPTION (2)
	install.packages('http://github.com/derekbeaton/outlieRs/raw/master/outlieRs_0.1.0.9003.tar.gz',repos=NULL)
