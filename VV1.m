function [ FX ] = VV1(X)
% This is the algorithm used to calculate
%    F(X) = (-X_1X3, -X_2X_3, 1-X_3^2) ;



     FX = -[ -X(1)*X(3), -X(2)*X(3), 1-X(3)^2 ]' ;


end

