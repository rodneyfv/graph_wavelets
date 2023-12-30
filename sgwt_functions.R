# Functions to perform the Spectral graph wavelet transform (SGWT) proposed by
# Hammond et al (2011)

# The following functions are R implementations based on the Matlab toolbox
# written by Hammond et al (2011), which is available at the link below:
# https://wiki.epfl.ch/sgwt

##############################################################################

# package used to read an image file
if(!require(OpenImageR)){
  install.packages("OpenImageR")
  library(OpenImageR)
}else{
  require(OpenImageR)
}

# package used to create a sparse matrix
if(!require(Matrix)){
  install.packages("Matrix")
  library(Matrix)
}else{
  require(Matrix)
}

# package with the eigs function, used to compute the largest eigenvalue of the
# Laplacian matrix
if(!require(rARPACK)){
  install.packages("rARPACK")
  library(rARPACK)
}else{
  require(rARPACK)
}

# package with the dfsane function, which solves nonlinear system of equantions.
# It is used to solve the problem (W^*W)W^*f=c, to find the pseudo-inverse of
# the SGWT
if(!require(BB)){
  install.packages("BB")
  library(BB)
}else{
  require(BB)
}

##############################################################################

# Spectral graph wavelet transform computed using the full eigendecomposition
# of the Laplacian matrix
# Inputs
# gft: graph Fourier transform of the signal
# mU: matrix with the eigenvectors of the Laplacian
# vLambda: vector with the eigenvalues of the Laplacian
# vt: scale factors used in the kernel
# Output
# mW: matrix whose rows correspond to each scale of the kernel. Columns
# correspond to vertices
sgwt0 = function(gft,mU,vLambda,vt){
  M = length(vt) + 1  # total number of scales in the SGWT
  N = length(vLambda)  # number of vertices
  mW = matrix(0,M,N)  # matrix for the SGWT coefficients
  for(m in 1:M){
    # for this part, check section 4.1 of Hammond et al (2011)
    # kernel or scaling function applied to eigenvalues
    vgm = sgwt_filter_kernel(vLambda,(m-1),vt,max(vLambda))
    mW[m,] = t(mU%*%(vgm*gft))  # SGWT coefficient for scale m
  }
  return(mW)
}

# Adjoint of the spectral graph wavelet transform computed using the full 
# eigendecomposition of the Laplacian matrix
# Inputs
# mW: matrix with SGWT coefficients for each scale and vertices
# mU: matrix with the eigenvectors of the Laplacian
# vLambda: vector with the eigenvalues of the Laplacian
# vt: scale factors used in the kernel
# Output
# vector with the coefficients corresponding to the adjoint of mW
sgwt0_adjoint = function(mW,mU,vLambda,vt){
  M = length(vt) + 1  # total number of scales in the SGWT
  N = length(vLambda)  # number of vertices
  mW2 = matrix(0,M,N) # matrix for coefficients of the SGWT's adjoint
  for(m in 1:M){
    # for this part, check section 6.1 of Hammond et al (2011)
    # kernel or scaling function applied to eigenvalues
    vgm = sgwt_filter_kernel(vLambda,(m-1),vt,max(vLambda))
    # we re-apply the corresponding wavelet operator to each subband of mW
    mW2[m,] = t(mU%*%(vgm*(t(mU)%*%mW[m,])))
  }
  # the adjoint is given by the sum of subbands in mW2 in all scales
  return(apply(mW2,2,sum))
}

# (Generalized) inverse of the spectral graph wavelet transform, computed using 
# the full eigendecomposition of the Laplacian matrix
# Inputs
# mW: matrix with SGWT coefficients for each scale and vertices
# mU: matrix with the eigenvectors of the Laplacian
# vLambda: vector with the eigenvalues of the Laplacian
# vt: scale factors used in the kernel
# Output
# vector with the inverse transform corresponding to mW's coefficients
sgwt0_inverse = function(mW,mU,vLambda,vt){
  M = length(vt) + 1  # total number of scales in the SGWT
  N = length(vLambda)  # number of vertices
  vadj = sgwt0_adjoint(mW,mU,vLambda,vt)  # adjoint of mW
  # for this part, check section 7 of Hammond et al (2011)
  # function to compute (nW^*mW)x for a given signal vector x, and return its
  # difference with the vectors vadj above
  wstarw = function(x){
    gft_x = t(mU)%*%x  # GFT of x
    mWx = sgwt0(gft_x,mU,vLambda,vt)  # SGWT of x
    vadjx = sgwt0_adjoint(mWx,mU,vLambda,vt)  # adjoint of the SGWT of x
    return(vadjx - vadj)
  }
  # solving the pseudo-inverse proble (mW^*mW)x=vadj
  vr = dfsane(par=vadj/2,fn=wstarw)
  
  return(vr$par)
}

# Blocks test function of Donoho and Johnstone (1994)
fBlocks = function(t){
  vt = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
  vh = c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
  n = length(t)
  y = rep(NA,n)
  for(i in 1:n) y[i] = sum(vh*(1+sign(t[i] - vt))/2)
  return(y)
}


# Bumps test function of Donoho and Johnstone (1994)
fBumps = function(t){
  vt = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
  vh = c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
  vw = c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)
  n = length(t)
  y = rep(NA,n)
  for(i in 1:n) y[i] = sum(vh*((1+abs((t[i] - vt)/vw))^(-4)))
  return(y)
}

# Heavisine test function of Donoho and Johnstone (1994)
fHeavysine = function(t){
  4*sin(4*pi*t) - sign(t - 0.3) - sign(0.72 - t)
}

# Doppler test function of Donoho and Johnstone (1994)
fDoppler = function(t){
  sqrt(t*(1-t))*sin(2*pi*1.05/(t + 0.05))
}


# function to compute the adjacency matrix of a 2d mesh. For each element (pixel),
# the neighbors above, below, on the left and on the right receive weight one,
# and all the others positions receive zero.
# Inputs
# dimmat: dimension of the 2d mesh
# Output
# mA: An adjacency matrix. A sparse matrix is returned.
sgwt_mesh_rectangle = function(dimmat){
  # case when the dimension of the square matrix is given as one number
  if(length(dimmat)==1) dimmat = dimmat*c(1,1)
  # total number of pairs the matrix has
  N = prod(dimmat)
  # all possible combinations of pairs, e.g., if dim=(n,n), then it gives a
  # matrix with rows (1,1),...,(1,n),(2,1),...,(2,n),...,(n,1),...,(n,n)
  alli = rep(1:dimmat[1],times=dimmat[1])
  allj = rep(1:dimmat[1],each=dimmat[1])
  # (ci(k),cj(k)) has neighbor (ni(k),nj(k)). 
  # Imagine the matrix [ci,cj,ni,nj]. The first half of the rows
  # are the neighbors in the column at the right, and the second half are
  # the neighbors at the row below
  ci = c(alli,alli)
  cj = c(allj,allj)
  ni = c(alli, alli+1)
  nj = c(allj+1, allj)
  # since we basically added one to the columns and rows indexes, we will
  # get indexes greater than the actual number of rows and columns. Here
  # we get rid of these numbers
  valid = (ni>=1 & ni<=dimmat[1] & nj>=1 & nj<=dimmat[2])
  ni = ni[valid]
  nj = nj[valid]
  ci = ci[valid]
  cj = cj[valid]
  # index of a vectorized version of the matrix of indexes, going row by
  # row, e.g., (1,1)=1,(1,2)=2,...,(1,n)=n,(2,1)=n+1,...
  cind=dimmat[1]*(cj-1)+ci
  nind=dimmat[1]*(nj-1)+ni
  # The dimension of A will the the double of
  # the number of elements in ni, because, e.g., (1,1) is nighbor of (1,2),
  # but (2,1) is also neighbor of (1,1). This makes A be a symmetric
  # matrix. Each row corresponds to a element (pixel), and the columns will
  # have a one on its neighbors and zero otherwise.
  mA = Matrix(0,N,N, sparse=TRUE)
  mA[cbind(c(cind,nind),c(nind,cind))] = 1
  
  return(mA)
}

# Laplacian of an adjacency matrix mA
# Inputs
# mA: an adjacency matrix
# case: when "raw", the Laplacian is computed as (D-A), and when "normalized", 
# it is computed as D^(-1/2)*(D-A)*D^(-1/2)
# Output
# mL: the Laplacian of mA
sgwt_laplacian = function(mA,case="raw"){
  N = dim(mA)[1]  # sample size
  vD = colSums(mA)  # degree of each vertex
  diagw = diag(mA)  # diagonal of mA
  
  # getting indexes of the values different of 0 in mA
  mKij = which(mA !=0, arr.ind = T)
  # checking if there is a loop (nonzero diagonal element in mA)
  ndind = which(mKij[,1]!=mKij[,2])
  # indexes of rows and columns and weights corresponding to these indexes
  ni = mKij[ndind,1]
  nj = mKij[ndind,2]
  w = mA[cbind(ni,nj)]
  di = c(1:N)  # diagonal indices
  
  mL = Matrix(0,N,N, sparse=TRUE)
  if(case=="raw"){
    # raw laplacian (D-A)
    mL[cbind(c(ni,di),c(nj,di))] = c(-w,vD - diagw)
  }else{
    if(case=="normalized"){
      # normalized laplacian D^(-1/2)*(D-A)*D^(-1/2)
      # diagonal entries
      dL=(1-diagw/vD) # will produce NA for vD==0 locations
      dL[vD==0] = 0  # which will be fixed here
      # nondiagonal entries
      ndL = -w/c( sqrt(vD[ni]*vD[nj]) );
      mL[cbind(c(ni,di),c(nj,di))] = c(ndL,dL)
    }else{
      cat("choose raw or normalized")
      return(NULL)
    }
  }
  return(mL)
}

# Wavelet function kernel.
# Inputs
# x: a vector of points to evaluate the kernel g
# alpha, beta, t1, t2: See Eq. (65) of Hammond et al (2011)
# Outputs
# y: a vector of g evaluated at the points of x
filter_g = function(x,alpha=2,beta=2,t1=1,t2=2){
  # variables to store the index of values that lie 
  # on the interval (t1,t2) or not
  r1 = which(x>=0 & x<t1)
  r2 = which(x>=t1 & x<t2)
  r3 = which(x>=t2)
  # getting the coefficients of the cubic spline fit in (t1,t2)
  M = matrix(c(1, t1, t1^2, t1^3,
               1, t2, t2^2, t2^3,
               0, 1, 2*t1, 3*t1^2,
               0, 1, 2*t2, 3*t2^2),4,4,byrow=TRUE)
  v=c(1, 1, (t1^(-alpha))*alpha*(t1^(alpha-1)), -beta*(t2^(-beta-1))*(t2^beta))
  a=solve(M,v)
  # computing the wavelet function according to where x_j lies
  y=x
  y[r1] = (x[r1]/t1)^alpha
  y[r2] = a[1]+a[2]*x[r2]+a[3]*(x[r2]^2)+a[4]*(x[r2]^3)
  y[r3] = (t2/x[r3])^beta
  return(y)
}

# scaling function kernel.
# Inputs
# x: a vector of points to evaluate the kernel g
# lmax: largest eigenvalue of the graph Laplacian
# Nscales: number of scales to evaluate the wavelet kernel
# lpfactor: See Section 8.1 of Hammond et al (2011)
# Outputs
# a vector of h evaluated at the points of x
filter_h = function(x,lmax,Nscales,lpfactor=20){
  gamma = 1.3849 # maximum value of filter_g
  lmin = lmax/lpfactor
  return(gamma*exp(-(x/(0.6*lmin))^4))
}

# Wavelet and scaling function considering the scale values used
# Inputs
# x: a vector of points to evaluate the kernel g(x*t_jj) or h(x)
# jj: wavelet scale to evaluate with the wavelet kernel g
# vt: vector with the scale values used
# lmax: largest eigenvalue of the graph Laplacian
# Nscales: number of scales to evaluate the wavelet kernel
# lpfactor: See Section 8.1 of Hammond et al (2011)
# Remaining parameters can be consulted on Section 8.1 of Hammond et al (2011)
# Outputs
# a vector of g(x*t_jj) or h(x) evaluated at the points of x
sgwt_filter_kernel = function(x,jj,vt,lmax,lpfactor=20,
                              alpha=2,beta=2,t1=1,t2=2){
  Nscales = length(vt)
  if(jj==0){
    # scaling kernel
    return(filter_h(x,lmax,Nscales,lpfactor))
  }else{
    if(jj<=Nscales){
      # wavelet kernel at scale jj
      return(filter_g(vt[jj]*x,alpha,beta,t1,t2))
    }else{
      cat("invalid scale index j")
      return(NULL)
    }
  }
}

# Function to compute the Chebyshev coefficients. 
# Inputs
# g: kernel function to approximate with Chebyshev polynomials
# m0: maximum order Chebyshev coefficient to compute
# m1: grid order used to compute quadrature (default is m+1)
# lmax: largest eigenvalue of the graph Laplacian
# Outputs
# vc: a vector of coefficients of a truncated Chebyshev expansion
sgwt_cheby_coeff = function(g,m0=50,m1=NULL,lmax){
  if(is.null(m1)) m1 = m0 + 1
  a_lmax = lmax/2
  # Computing the Chebyshev coefficients through numerical integration. Check
  # Eq. (58) of Hammond et al (2011)
  vc = rep(NA,m0 + 1)
  for(j in 1:(m0+1)){
    vc[j] = sum( g(a_lmax*cos( (pi*(c(1:m1)-0.5))/m1) + a_lmax) * 
                   cos(pi*(j-1)*(c(1:m1)-.5)/m1) )*2/m1
  }
  return(vc)
}

# vector version of the delta function
# Inputs
# N: dimension of the vector
# j: position where we have one, having zero elsewhere
# Output
# f: a Nx1 vector with one in the j-th entry and zero elsewhere
sgwt_delta = function(N,j){
  f = rep(0,N)
  f[j] = 1
  return(f)
}

# Chebyshev polynomial of Laplacian applied to vector
# For more details, check Section 6 of Hammond et al (2011). 
# Inputs
# f: vector of values obtained for each vertices of the graph
# mL: the graph Laplacian
# mc: list with Chebyshev coefficients for scaling and wavelet kernels
# lmax: largest eigenvalue of the graph Laplacian
# Outputs
# mr: a list with the SGWT coefficients for each kernel (scaling and wavelets) 
# considered
sgwt_cheby_op = function(f,mL,mc,lmax){
  a1 = lmax/2
  # total number of kernels, scaling and wavelets
  Nscales = length(mc)
  
  # number of coefficients computed for each kernel
  vM = rep(0,Nscales)
  for(j in 1:Nscales) vM[j] = length(mc[[j]])
  # maximum number of Chebyshev coefficients.
  maxM = max(vM)
  
  mr = list(0) # list for the SGWT coefficients
  # initial operators of the recurrence relation.
  Twf_old = f # j=0
  Twf_cur = (mL%*%f - a1*f)/a1 # j=1
  # First two terms of the series that computes the wavelet coefficients. Check
  # Eq. (57) of Hammond et al (2011)
  for(j in 1:Nscales){
    mr[[j]] = .5*mc[[j]][1]*Twf_old + mc[[j]][2]*Twf_cur
  }
  # Remaining terms of the series that computes the wavelet coefficients
  for(k in 2:maxM){
    # recurrence relation of the Chebyshev polynomials. Eq. (58)
    Twf_new = (2/a1)*(mL%*%Twf_cur-a1*Twf_cur) - Twf_old
    for(j in 1:Nscales){
      if(k+1<=vM[j]) mr[[j]] = mr[[j]] + mc[[j]][k+1]*Twf_new
    }
    # updating the current and new Chebyshev polynomials
    Twf_old = Twf_cur
    Twf_cur = Twf_new
  }
  return(mr)
}

# Compute a set of wavelet scales adapted to spectrum bounds
# Inputs
# lmax: largest eigenvalue of the graph Laplacian
# Nscales: number of scales to evaluate the wavelet kernel
# lpfactor: See Section 8.1 of Hammond et al (2011)
# Outputs
# vt: a vector of logarithmically equispaced values
sgwt_setscales = function(lmax,Nscales,lpfactor=20,t1=1,t2=2){
  lmin = lmax/lpfactor
  smin = t1/lmax
  smax = t2/lmin
  vt = exp(seq(log(smax),log(smin),length = Nscales))
  return(vt)
}


# Compute adjoint of SGWT. See eq. (59) of Hammond et al (2011)
# Inputs
# vy: vector of SGWT coefficients
# mL: the graph Laplacian
# mc: list with Chebyshev coefficients for scaling and wavelet kernels
# lmax: largest eigenvalue of the graph Laplacian
# Outputs
# v_adj: a vector with the adjoint coefficients of the SGWT
sgwt_adjoint = function(vy,mL,mc,lmax){
  N = dim(mL)[1]
  Nc = length(mc)
  v_adj = rep(0,N)
  for(j in 1:Nc){
    tmp = sgwt_cheby_op(vy[[j]],mL,list(mc[[j]]),lmax)
    v_adj = v_adj + tmp[[1]]
  }
  return(v_adj)
}

# Chebyshev coefficients for square of polynomial of the kernel corresponding to
# a vector of Chebyshev coefficients vc
# See Eq. (62) of Hammond et al (2011)
# Inputs
# vc: vector with Chebyshev coefficients for a kernel (scaling or wavelets)
# Outputs
# vdp: vector of Chebyshev coefficients for square of polynomial
sgwt_cheby_square = function(vc){
  # chebyshev polynomial order used in the truncation
  M = length(vc) - 1
  vcp = vc
  vcp[1] = 0.5*vc[1]
  # adjust cp so that p(x) = sum cp(1+k) T_k(x)
  # for all k>=0 (rather than with special case for k=0)
  vdp = rep(0,2*M + 1)
  for(m in 0:(2*M)){
    if(m==0){
      vdp[1+m] = vdp[1+m] + .5*vcp[1]^2
      for(i in 0:M) vdp[1+m] = vdp[1+m] + .5*vcp[i+1]^2
    }else{
      if(m<=M){
        for(i in 0:m) vdp[1+m] = vdp[1+m] + .5*vcp[i+1]*vcp[m-i+1]
        for(i in 0:(M-m)) vdp[1+m] = vdp[1+m] + .5*vcp[i+1]*vcp[m+i+1]
        for(i in m:M) vdp[1+m] = vdp[1+m] + .5*vcp[i+1]*vcp[i-m+1]
      }else{
        for(i in (m-M):M) vdp[1+m] = vdp[1+m] + .5*vcp[i+1]*vcp[m-i+1]
      }
    }
  }
  vdp[1] = 2*vdp[1]
  return(vdp)
}

# Compute inverse SGWT via conjugate gradients
# Inputs
# vy: list with SGWT coefficients
# mL: the graph Laplacian
# mc: list with Chebyshev coefficients for scaling and wavelet kernels
# lmax: largest eigenvalue of the graph Laplacian
# Outputs
# vr: a vector containing the (pseudo) inverse SGWT of vy
sgwt_inverse = function(vy,mL,mc,lmax){
  N = dim(mL)[1]
  # computing coefficients of adj = W^*y
  vadj = sgwt_adjoint(vy,mL,mc,lmax)
  # computing coefficients of (W^*W)
  # total number of kernels, scaling and wavelets
  Nscales = length(mc)
  # number of coefficients computed for each kernel
  vM = rep(0,Nscales)
  for(j in 1:Nscales) vM[j] = length(mc[[j]])
  maxM = max(vM)
  # computing P(x) = p(x)^2, see Eq. (63) of Hammond et al (2011)
  vd = rep(0,2*maxM-1)
  for(j in 1:Nscales){
    cpad = rep(0,maxM)
    cpad[c(1:vM[j])] = mc[[j]]
    # denoting the order of the Chebyshev polynomial by mj, one can note that
    # sgwt_cheby_square will return a vector of length
    # 2*mj+1 = 2*(mj+1)-1 = 2*vM[j]-1
    vd = vd + sgwt_cheby_square(cpad)
  }
  ## matrix corresponding to the operator (W^*W)
  #wstarw = sgwt_cheby_op(.sparseDiagonal(N),mL,list(vd),lmax)
  ## using conjugate gradients to solve (W^*W)f=adj for f
  #vr = solve(wstarw[[1]],vadj)
  wstarw = function(x) as.vector(sgwt_cheby_op(x,mL,list(vd),lmax)[[1]]) - as.vector(vadj)
  vr = dfsane(par=vadj/2,fn=wstarw)
  
  return(vr$par)
}
