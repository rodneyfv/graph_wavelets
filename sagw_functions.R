# Functions to perform the Spectrum-adapted graph wavelet 
# transform (SAGWT) proposed by Shuman et al (2015).
# It requires the implemented SGWT functions, used to apply the
# method of Hammond et al (2011).

# The original Matlab codes prepared by Shuman et al (2015) are available 
# at the link below:
# http://faculty.olin.edu/dshuman/publications.html

##############################################################################

# Function to compute spectral adapted graph wavelet coeffients using the 
# method proposed by Shuman et al (2015) and the Chebyshev approximation
# proposed by Hammond et al (2011).
# x:  vector of values obtained for each vertices of the graph
# mA: the graph adjacency matrix
# M: number of filters
# R: some integer in {2,...,M}
# Q: number of points used to approximate the spectral cdf
compute_coefficients_sagwt <- function(vSignal, mA, M=5, R=3, Q=24){
  # sample size
  n = dim(mA)[1]
  # computing the graph Laplacian
  mL = sgwt_laplacian(mA)
  # computing the largest eigenvalue. We need to convert mL the class "dgCMatrix"
  # in order to use the function eigs
  lmax = eigs(as(mL,"dgCMatrix"), 1, which = "LM")$values
  
  # warping function generated from an approximation of the 
  # cumulative spectral density function of the graph 
  # Laplacian eigenvalues
  spectral_warp_fn = sagw_spectrum_cdf_approx(mL,lmax,Q)
  
  m0=40 # chebyshev polynomial order, for approximation
  # computing a list with the Chebyshev coefficients
  mc = list(0)
  # wavelet coefficients
  for(m in 1:(M-1)){
    mc[[m]] = sgwt_cheby_coeff(
      function(x) sagw_wavelet_filter_approx(x,m+1,mL,
                                             spectral_warp_fn,lmax,M,R),
      m0=m0,m1=1000,lmax)
  }
  # scaling coefficients
  mc[[M]] = sgwt_cheby_coeff(
    function(x) sagw_scaling_filter_approx(x,mL,
                                           spectral_warp_fn,lmax,M,R),
    m0=m0,m1=1000,lmax)
  
  # Computing the wavelet graph transform of the signal
  vW = sgwt_cheby_op(vSignal - mean(vSignal),mL,mc,lmax)

  return(vW = vW)
}


# Hann kernel, defined in the Corollary 1 of Shuman et al (2015),
# with K=1 and a_0=a_1=1/2
# Inputs
# y: vector of values to evaluate the function at
h_Hann = function(y){
  m = length(y)
  r1 = which(y>=0 & y<1)
  hh = rep(0,m)
  hh[r1] = .5 + .5*cos(2*pi*(y[r1] - .5))
  return(hh)
}

# Function used in the frame bounds obtained with the shifted
# Hann kernels. It is constant on the interval [0,(M+1-R)/R]
# Inputs
# y: vector of values to evaluate the function at
# M: number of filters
# R: some integer in {2,...,M}
H_Hann = function(y,M,R){
  m = length(y)
  r1 = which(y>=0 & (y<=(M+1-R)/R))
  hh = rep(0,m)
  for(i in 1:length(r1)){
    hh[r1[i]] = sum(h_Hann(y[r1[i]] - c((1-R):(M-R))/R)^2)
  }  
  return(hh)
}

# kernel function given in Eq. (9) of Shuman et al (2015), obtained
# from the Hann kernel. It is the same as used in the toobox provided
# by the authors, which is slightly different from the expression
# presented in the paper.
# Inputs
# x: vector of values to evaluate the function at
# m: index of the translated kernel
# gama: upper bound of the values in x
# M: number of filters
# R: some integer in {2,...,M}
gUhatm = function(x,m,gama,M,R){
  n = length(x)
  dilation = R*gama/(M+1-R)
  r1 = which(x>=0 & x<=gama)
  gg = rep(0,n)
  gg[r1] = h_Hann(x[r1]/dilation - (m-R)/R)
  return(gg)
}

# Frame bound function given in Eq. (10) of Shuman et al (2015).
# Inputs
# x: vector of values to evaluate the function at
# gama: upper bound of the values in x
# M: number of filters
# R: some integer in {2,...,M}
Glamb = function(x,gama,M,R){
  n = length(x)
  r1 = which(x>=0 & x<=gama)
  gg = rep(0,n)
  for(m in 1:M) gg[r1] = gg[r1] + gUhatm(x[r1],m,gama,M,R)^2
  return(gg)
}

# wavelet kernel with log wrap given in Eq. (12) of Shuman et al (2015), 
# using Hann kernel.
# Inputs
# x: vector of values to evaluate the function at
# m: index of the translated kernel. Must have value in {2,...,M}
# lmax: upper bound of the values in x
# M: number of filters
# R: some integer in {2,...,M}
ghat_log_wavelet = function(x,m,lmax,M,R){
  n = length(x)
  gama = log(lmax) + 1e-10
  r1 = which(x>0 & x<lmax)
  gg = rep(0,n)
  tmp = log(x[r1]) + 1e-10
  gg[r1] = gUhatm(tmp,m-1,gama,M,R)
  return(gg)
}

# scaling kernel with log wrap given in Eq. (12) of Shuman et al (2015), 
# using Hann kernel.
# Inputs
# x: vector of values to evaluate the function at
# lmax: upper bound of the values in x
# M: number of filters
# R: some integer in {2,...,M}
ghat_log_scaling = function(x,lmax,M,R){
  n = length(x)
  gama = log(lmax) + .001
  r1 = which(x>0 & x<lmax)
  gg = rep(0,n)
  gg[r1] = 3*R/8
  for(m in 2:M) gg[r1] = gg[r1] - ghat_log_wavelet(x[r1],m,lmax,M,R)^2
  r1 = which(abs(gg)<1e-15)
  gg[r1] = 0
  return(sqrt(gg))
}

# Spectral adapted wavelet filter using the empirical spectral distribution
# of the graph Laplacian eigenvalues. See section 5.A of Shuman et al (2015)
# Inputs
# x: vector of values to evaluate the function at
# m: index of the translated kernel. Must have value in {2,...,M}
# vLambda: vector of eigenvalues of the graph Laplacian
# M: number of filters
# R: some integer in {2,...,M}
# method: interpolation method to apply at points of the empirical
# spectral density. It can be 'linear' or 'monocub' for linear or
# monotonic cubic interpolation, respectivelly.
specadapt_wavelet_filter = function(x,m,vLambda,M,R,method){
  n = length(x)  # number of points to evaluate the filter at
  N = length(vLambda)  # number of vertices in the graph
  # interpolating the cumulative spectral density function
  if(method=="linear"){
    P_lamb = approxfun(vLambda,c(0:(N-1))/N,method="linear")
  }else{
    if(method=="monocub"){
      P_lamb = splinefun(vLambda,c(0:(N-1))/N,method="monoH.FC")
    }else{
      cat("Interpolation method invalid.")
      return(NULL)
    }
  }
  # using the interpolated spectral cdf as warping function to the
  # wavelet filter. See Eq. (12) of Shuman et al (2015)
  lmax = max(vLambda)
  gama = P_lamb(lmax)
  r1 = which(x>0 & x<lmax)
  gg = rep(0,n)
  tmp = P_lamb(x[r1])
  gg[r1] = gUhatm(tmp,m-1,gama,M,R)
  return(gg)
}

# Spectral adapted wavelet filter using the empirical spectral distribution
# of the graph Laplacian eigenvalues. See section 5.A of Shuman et al (2015)
# Inputs
# x: vector of values to evaluate the function at
# vLambda: vector of eigenvalues of the graph Laplacian
# M: number of filters
# R: some integer in {2,...,M}
# method: interpolation method to apply at points of the empirical
# spectral density. It can be 'linear' or 'monocub' for linear or
# monotonic cubic interpolation, respectivelly.
specadapt_scaling_filter = function(x,vLambda,M,R,method){
  n = length(x)  # number of points to evaluate the filter at
  N = length(vLambda)  # number of vertices in the graph
  # interpolating the cumulative spectral density function
  if(method=="linear"){
    P_lamb = approxfun(vLambda,c(0:(N-1))/N,method="linear")
  }else{
    if(method=="monocub"){
      P_lamb = splinefun(vLambda,c(0:(N-1))/N,method="monoH.FC")
    }else{
      cat("Interpolation method invalid.")
      return(NULL)
    }
  }
  # using the interpolated spectral cdf as warping function to the
  # scaling filter. See Eq. (13) of Shuman et al (2015)
  lmax = max(vLambda)
  gama = P_lamb(lmax)
  r1 = which(x>0 & x<lmax)
  gg = rep(0,n)
  gg[r1] = 3*R/8
  for(m in 2:M){
    gg[r1] = gg[r1] - specadapt_wavelet_filter(x[r1],m,vLambda,M,R,method)^2
  }
  # setting low absolute values to zero to avoid problem with the square root
  r1 = which(abs(gg)<1e-15)
  gg[r1] = 0
  return(sqrt(gg))
}

# Approximation of the cumulative spectral density function
# with a spectrum slicing method. See section 5.B
# of Shuman et al (2015).
# Inputs
# mL: graph Laplacian matrix
# lmax: upper bound of the values in x (largest eigenvalue)
# Q: number of points used to approximate the spectral cdf
# Output
# a function that takes one number and returns the approximated
# cdf evaluated at it
sagw_spectrum_cdf_approx = function(mL,lmax,Q=24){
  N = dim(mL)[1]  # number of vertices in the graph
  # sequence of points used in the sprectrum slicing method
  vq = seq(0,Q,1)*lmax/Q
  # vector to store the number of negative eigenvalues in the
  # diagonal matrix of the LDL' decomposition of mL-(q*lmax/Q)*mI
  vCounts = rep(0,Q+1)
  vCounts[Q+1] = N-1
  # sparse identity matrix
  mIdentity = Matrix(0,N,N, sparse=TRUE)
  diag(mIdentity) = 1
  for(i in 2:Q){
    mL0 = mL - vq[i]*mIdentity
    # using the Cholesky function of the package Matrix. That's
    # why we need to transform mL0 into a sparse matrix
    mL_chol = Cholesky(as(mL0,"dgCMatrix"))
    # The objetc .@x contains nonzero values of the Cholesky
    # decomposition and .@p contain the initial (zero-based) index
    # of elements of each colum. These initial indexes correspond
    # to the values in .@x in the diagonal matrix of the LDL'
    # decomposition, and that's all we need.
    vDelta = mL_chol@x[mL_chol@p[1:N]+1]
    # using the method described on Section 5.B of Shuman et al (2015)
    vCounts[i] = sum(vDelta<0)
  }
  # interpolating the points obtained from the spectrum slicing method
  P_lamb = splinefun(vq,vCounts/(N-1),method="monoH.FC")
  return(P_lamb)
}

# Spectral adapted wavelet filter using an approximation of the cumulative
# spectral density function with a spectrum slicing method. See section 5.B
# of Shuman et al (2015).
# Inputs
# x: vector of values to evaluate the function at
# m: index of the translated kernel. Must have value in {2,...,M}
# mL: graph Laplacian matrix
# P_lamb: approximated spectral cdf
# lmax: upper bound of the values in x (largest eigenvalue)
# M: number of filters
# R: some integer in {2,...,M}
sagw_wavelet_filter_approx = function(x,m,mL,P_lamb,lmax,M,R){
  n = length(x)  # number of points to evaluate the filter at
  # using the approximated spectral cdf as warping function to the
  # wavelet filter. See Eq. (12) of Shuman et al (2015)
  r1 = which(x>0 & x<lmax)
  gg = rep(0,n)
  tmp = P_lamb(x[r1])
  gg[r1] = gUhatm(tmp,m-1,1,M,R)
  return(gg)
}

# Spectral adapted scaling filter using an approximation of the cumulative
# spectral density function with a spectrum slicing method. See section 5.B
# of Shuman et al (2015).
# Inputs
# x: vector of values to evaluate the function at
# mL: graph Laplacian matrix
# P_lamb: approximated spectral cdf
# lmax: upper bound of the values in x (largest eigenvalue)
# M: number of filters
# R: some integer in {2,...,M}
sagw_scaling_filter_approx = function(x,mL,P_lamb,lmax,M,R){
  n = length(x)  # number of points to evaluate the filter at
  # using the approximated spectral cdf as warping function to the
  # scaling filter. See Eq. (13) of Shuman et al (2015)
  gama = P_lamb(lmax)
  r1 = which(x>0 & x<lmax)
  gg = rep(0,n)
  gg[r1] = 3*R/8
  for(m in 2:M){
    gg[r1] = gg[r1] - (sagw_wavelet_filter_approx(x[r1],m,mL,P_lamb,lmax,M,R)^2)
  }
  r1 = which(abs(gg)<1e-10)
  gg[r1] = 0
  return(sqrt(gg))
}

