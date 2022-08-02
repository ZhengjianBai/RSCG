function [ FX ] = VV2(X)
% This is the algorithm used to calculate
%  F(X) ;


%                A = prod(X) ;
%            
%              B = A/(A + 1) ;
%            
%                C = X.*B ;
%            
%            FX =-0.05*C + X.*log(X)-X ;
           
           
                 A = prod(X) ;
           
%                B = A/(A + 1) ;
%            
%                 C = X.*B ;
           
           FX =-0.01*X + 0.01*X/(A+1) + X.*log(X)-X ;

end

