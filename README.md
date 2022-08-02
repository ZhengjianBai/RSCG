# GSCG-HM
Geometric Spectral Conjugate Gradient Method for Monotone Vector Fields on Hadamard Manifolds

This package includes Matlab codes for geometric spectral conjugate gradient algorthms.
They are efficient algorithms for approximating a zero of a monotone tangent vector field
on a constant curvature Hadamard manifold. 
The test cases and scripts for running the experiments in paper
"Geometric Spectral Conjugate Gradient Method for Monotone Vector Fields on Hadamard Manifold" 
by Zhi Zhao, Xiaoqing Jin, Zengjian Bai, and Tengteng Yao, are also included.


1. Main algorithms

RSCGPerry1.m -- geometric spectral conjugate gradient method using the Perry formula for Example 4.1

RSCGPerry2.m -- geometric spectral conjugate gradient method using the Perry formula for Example 4.2

2. Test data and codes

TRANS.m -- the vector transport on H^2 

PHI.m-- backtraking line search

LORENTZ.m -- the Riemannian metric on H^2 

INPRODUCT.m -- the Riemannian metric on R^n_++  

VV1.m -- the value of the underlying vector field of Example 4.1

VV2.m -- the value of the underlying vector field of Example 4.2

EXPINV1-- the inverse of the exponential mapping on H^2 

EXPP2-- the exponential mapping on R^n_++ 

EXPINV2-- the inverse of the exponential mapping on R^n_++ 
 
DIST1-- the Riemannian distance function on H^2 

DIST2-- the Riemannian distance function on R^n_++  

 

 