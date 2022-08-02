function [ Z, FZ, lfro, rfro ] = PHI(X, XK, ak)
% This is the algorithm used to calculate
%    F(X) = (-X_1X3, -X_2X_3, 1-X_3^2) ;


         cc =  sqrt(LORENTZ(XK,XK)) ;

             XK = XK/cc ;
             
             ak = cc*ak ;
        
        VX = sinh(ak)*X + cosh(ak)*XK ;
        
        Z = cosh(ak)*X + sinh(ak)*XK ;    

               FZ = VV1(Z)  ;


          lfro = -LORENTZ(FZ, VX) ;
  
                 rfro = ak ;
        

end

