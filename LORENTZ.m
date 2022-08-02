function [ ff ] = LORENTZ(VK, WK)
% This is the algorithm used to calculate
% the Lorentz metric ;



     ff =  VK(1:2)'*WK(1:2) - VK(3)*WK(3) ;



end