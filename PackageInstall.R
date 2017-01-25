# Install and Load package from GitHub

install.packages("devtools")
require(devtools)

install.packages(c('rrcov','doparallel','TInPosition','rpca','reshape','ggplot2','colorspace'),dependencies = T)

install_github("derekbeaton/outlieRs", subdir="Package")