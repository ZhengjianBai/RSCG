function [ SK ] = EXPINV(NX,X)
% This is the algorithm used to calculate


                cc = LORENTZ(NX,X) ;
         
         SK = (acosh(-cc)/sqrt(cc^2-1))*(X + cc*NX) ;

end

