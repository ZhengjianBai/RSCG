function [ XX ] = TRANS1(X, XK)
%  This is the algorithm used to calculate
%  the vector transport On the Riemannian Manifold;
%  Vector Transport Defined by Projection.


%% Translation on the compact Stiefel manifold ;
%  XX = XK - X*sym(X^T*XK) ;

         aa = sqrt(LORENTZ(XK,XK)) ;
  
         XX = XK + LORENTZ(X,XK)*X ;
         
%          XX = (aa/sqrt(LORENTZ(XX,XX)))*XX ;
         
    

end

