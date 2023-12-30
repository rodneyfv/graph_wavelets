# 

# File with functions to perform SGWT
source("./../sgwt_functions.R")
# File with functions to perform SAGWT
source("./../sagw_functions.R")

# Function to create the adjacency matrix of line graph with n_nodes nodes
create_adjacency_matrix_linegragh <- function(n_node){
  adj_matrix <- matrix(0,n_nodes,n_nodes)
  adj_matrix[1,2] = 1
  adj_matrix[n_nodes, n_nodes-1] = 1
  for(i in 2:(n_nodes-1)){ 
    adj_matrix[i,i-1] = 1
    adj_matrix[i,i+1] = 1
  }
  return(adj_matrix)
}

###############################################################
###############################################################
# Comparing graph kernels of Hammond et al. (2011) and Shuman et al. (2015)              
###############################################################
###############################################################

# creating the adjacency matrix of a path graph
n_nodes <- 64  # number of nodes

adjacencyMatrix_LineGraph <- create_adjacency_matrix_linegragh(n_nodes)

degreeMatrix_LineGraph <- diag(apply(adjacencyMatrix_LineGraph,1,sum))
# Laplacian
laplacianMatrix_LineGraph <- degreeMatrix_LineGraph - adjacencyMatrix_LineGraph
# eigen decomposition of the Laplacian, usind a non-increasing
# order for the eigenvalues
eigenvalues_laplacian <- eigen(laplacianMatrix_LineGraph)$values[n_nodes:1]
eigenvectors_laplacian <- eigen(laplacianMatrix_LineGraph)$vectors[,n_nodes:1]

# largest eigenvalue of the graph Laplacian
max_eigenvalue <- max(eigenvalues_laplacian)
# Sequence of equally spaced values in the graph spectrum. 
# We will make plots evaluating the wavelet kernels at these values.
seqence_graph_spectrum <- seq(0, max_eigenvalue, length=1000)

# number of scales (including approximation) to be used in the SGWT
numberScales <- 5

################

# Filters obtained using the wavelet filter of SGWT proposed
# by Hammond et al (2011)
scale_sgwt_filter <- sgwt_setscales(max_eigenvalue, Nscales=(numberScales-1))
seqence_graph_spectrum <- seq(0, max_eigenvalue, length=1000)

# function to save the plot
png(file="sgwt_filter.png", width=600, height=450)

par(mar=c(5,5,2,2))

values_scaling_filter <- sgwt_filter_kernel(x = seqence_graph_spectrum, jj = 0, 
                                            vt = scale_sgwt_filter, 
                                            lmax = max_eigenvalue)
plot(seqence_graph_spectrum, values_scaling_filter, 
     ylim=c(0,1.7), type='l', lwd = 2,
     ylab=expression(g[m](lambda)), xlab=expression(lambda),
     cex.lab = 2, cex.sub = 2)
for(m in 1:(numberScales-1)){
  values_wavelet_filter <- sgwt_filter_kernel(x = seqence_graph_spectrum, jj = m, 
                                              vt = scale_sgwt_filter, 
                                              lmax = max_eigenvalue)
  lines(seqence_graph_spectrum,values_wavelet_filter,
        ylim=c(0,1.5), type='l', lwd = 2)
}
points(eigenvalues_laplacian, rep(0,n_nodes), pch='x', col='red')
grid(nx = NULL, ny = NULL, lty = 1, col = "gray", lwd = 1)

title(sub="(a)", adj=1, line=3, font=2,cex.sub=2)

dev.off()

################

# Filters obtained using the wavelet filter of SAGWT proposed
# by Shuman et al (2015)
R <- 3

# function to save the plot
png(file="sagwt_filter.png", width=600, height=450)

par(mar=c(5,5,2,2))

values_scaling_filter <- specadapt_scaling_filter(x=seqence_graph_spectrum, 
                                                  vLambda = eigenvalues_laplacian, 
                                                  M = numberScales,
                                                  R = R,
                                                  method="monocub")
plot(seqence_graph_spectrum, values_scaling_filter, 
     type='l', ylim=c(0,1.25), lwd = 2,
     ylab=expression(g[m](lambda)), xlab=expression(lambda),
     cex.lab = 2, cex.sub = 2)

for(m in 2:numberScales){
  values_wavelet_filter <- specadapt_wavelet_filter(x=seqence_graph_spectrum, 
                                                    m=m, 
                                                    vLambda = eigenvalues_laplacian, 
                                                    M = numberScales,
                                                    R =R,
                                                    method="monocub")
  
  lines(seqence_graph_spectrum, values_wavelet_filter, type='l', lwd = 2)
}

points(eigenvalues_laplacian, rep(0,n_nodes), pch='x', lwd=2, col='red')  # eigenvalues
grid(nx = NULL, ny = NULL, lty = 1, col = "gray", lwd = 1)
title(sub="(b)", adj=1, line=3, font=2, cex.sub=2)

dev.off()


