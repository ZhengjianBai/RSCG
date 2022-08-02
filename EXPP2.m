function [ Z ] = EXPP2(X, XK, ak)
%UNTITLED2 

             VK = ak*XK ;
  
          Z = X.*exp(VK./X) ;

end

