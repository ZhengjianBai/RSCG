function [ ff ] = INPRODUCT(VK, WK, X)
% This is the algorithm used to calculate
% the Lorentz metric ;



        ff =  sum((VK.*WK)./(X.^2)) ;



end