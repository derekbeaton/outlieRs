# Install and Load package from GitHub
# Run this file only once!

install.packages(c('rrcov','doparallel','TInPosition','rpca','reshape','ggplot2','colorspace','robustbase'),dependencies = T)

install.packages("devtools")
require(devtools)
install_github("derekbeaton/outlieRs", subdir="Package")
