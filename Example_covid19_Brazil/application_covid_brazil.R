# 

# File with functions to perform SGWT
source("./../sgwt_functions.R")
# File with functions to perform SAGWT
source("./../sagw_functions.R")

require(RColorBrewer)
require(tidyverse)
require(plot3D)


###############################################################
# Covid data - October 2020
###############################################################


# Graph of Brazilian cities
citiesData <- readRDS('./covid_brazil_Oct2020.rds')
adjacency_matrix <- readRDS('./adjacencyMatrix_citiesBrazil.rds')

# checking for nodes without signals or isolated nodes. 
# We consider only data from October 6th, 2020 
node_degrees <- apply(adjacency_matrix, 2, sum)
id_valid_nodes <- which(!is.na(citiesData@data$cases1000people) & node_degrees!=0)
rm(node_degrees)

# number of cities used in our study
dim(citiesData@data[id_valid_nodes,])
# number of cities NOT used in our study (there are two NA's)
dim(citiesData@data[-id_valid_nodes,])

# sample size
n_nodes = length(id_valid_nodes)

# matrix with the indexes of the nonzero entries of mA
pairs_connected_nodes = which(adjacency_matrix[id_valid_nodes, id_valid_nodes] != 0, 
                              arr.ind = TRUE)
# number of edges in the graph
dim(pairs_connected_nodes)[1]/2

# plot of the signals in the graph nodes
horizontal_coord <- sp::coordinates(citiesData)[id_valid_nodes,1]
vertical_coord <- sp::coordinates(citiesData)[id_valid_nodes,2]

numberCases_1000people <- as.numeric(citiesData@data$cases1000people[id_valid_nodes])
hist(numberCases_1000people, breaks = 20)
hist(sqrt(numberCases_1000people), breaks = 20)

# plot of the number of cases
plot3D::scatter2D(horizontal_coord, vertical_coord, 
                  colvar=as.numeric(sqrt(numberCases_1000people)),
                  pch=16,cex=.5, 
                  col=rev(brewer.pal(9,name = "Spectral")), 
                  axes=FALSE, ann = FALSE, frame.plot=FALSE)
title(main='Square root of number of cases per 1000 people',cex.main=1)
#dev.off()


###############################################################
# Spectral graph wavelet analysis
###############################################################

# computing the graph Laplacian
graph_laplacian = sgwt_laplacian(mA =  adjacency_matrix[id_valid_nodes,id_valid_nodes])
# computing the largest eigenvalue. We need to convert mL the class "dgCMatrix"
# in order to use the function eigs
max_eigenvalue = eigs(as(graph_laplacian,"dgCMatrix"), 1, which = "LM")$values

# warping function generated from an approximation of the 
# cumulative spectral density function of the graph 
# Laplacian eigenvalues of the Minnesota graph
spectral_warp_fn = sagw_spectrum_cdf_approx(mL =  graph_laplacian, 
                                            lmax = max_eigenvalue,
                                            Q = 24)

# number of scales (including approximation) to be used in the SGWT
M = 10
R = 3

# Plot of the empirical distribution of the graph spectrum
seqence_graph_spectrum = seq(0,max_eigenvalue,.1)

par(mar=c(5,5,3,3))
plot(seqence_graph_spectrum, 
     spectral_warp_fn(seqence_graph_spectrum), 
     type='l', ylab="Empirical distribution function",
     xlab=expression(lambda), cex.lab=1.5)
title(sub="(a)", adj=1, line=3, font=2,cex.sub=1.5)

# Plot of spectrum adapted filters used to analyze the Brazilian graph
seqence_graph_spectrum = seq(0,max_eigenvalue,.1)
par(mar=c(5,5,3,3))
values_scaling_filter <- sagw_scaling_filter_approx(x = seqence_graph_spectrum, 
                                                    mL = graph_laplacian, 
                                                    P_lamb = spectral_warp_fn, 
                                                    lmax = max_eigenvalue,
                                                    M = M, R = R)
plot(seqence_graph_spectrum, values_scaling_filter,
     type='l', ylim=c(0,1.25), xlab="x", ylab=expression(g[m](x)), cex.lab=1.5)
for(m in 2:M){
  values_wavelet_filter <- sagw_wavelet_filter_approx(x = seqence_graph_spectrum,
                                                      m = m,
                                                      mL = graph_laplacian,
                                                      P_lamb = spectral_warp_fn,
                                                      lmax = max_eigenvalue,
                                                      M = M, R = R)
  lines(seqence_graph_spectrum, values_wavelet_filter, type='l')
}
title(sub="(b)", adj=1, line=3, font=2,cex.sub=1.5)


# Graph wavelet coefficients
wavelet_coefficients <- compute_coefficients_sagwt(numberCases_1000people,
                                                   adjacency_matrix[id_valid_nodes,id_valid_nodes], 
                                                   M=M, R=R, Q=24)
# Plot of the graph wavelet coefficients
plot3D::scatter2D(horizontal_coord, vertical_coord, 
                  colvar=as.numeric(wavelet_coefficients[[1]]),
                  pch=16,cex=.5, 
                  col=rev(brewer.pal(9,name = "Spectral")), 
                  axes=FALSE, ann = FALSE, frame.plot=FALSE)

plot3D::scatter2D(horizontal_coord, vertical_coord,
                  colvar=as.numeric(wavelet_coefficients[[10]]),
                  pch=16,cex=.5, 
                  col=rev(brewer.pal(9,name = "Spectral")), 
                  axes=FALSE, ann = FALSE, frame.plot=FALSE)


###############################################################
# Preparing plots
###############################################################

# Plot of number of cases
png(file="covid_cases_brazil.png", width=450, height=450)

par(mar=c(3,2,5,2))
plot3D::scatter2D(horizontal_coord, vertical_coord, 
                  colvar=as.numeric(sqrt(numberCases_1000people)),
                  pch=16,cex=.5, 
                  col=rev(brewer.pal(9,name = "Spectral")), 
                  axes=FALSE, ann = FALSE, frame.plot=FALSE)
title(main='Square root of number of\n cases per 1000 people',cex.main=2)
title(sub="(a)", adj=1, line=1, font=2,cex.sub=2)

dev.off()

# Plot of graph wavelet coefficients
png(file="sagwt_covid_cases_brazil.png", width=450, height=450)

par(mar=c(3,2,5,2))
plot3D::scatter2D(horizontal_coord, vertical_coord, 
                  colvar=as.numeric(wavelet_coefficients[[1]]),
                  pch=16,cex=.5, 
                  col=rev(brewer.pal(9,name = "Spectral")), 
                  axes=FALSE, ann = FALSE, frame.plot=FALSE)
title(main='Wavelet coefficients at lowest scale',cex.main=2)
title(sub="(b)", adj=1, line=1, font=2,cex.sub=2)

dev.off()
