function [ dd ] = DIST2(XS, X)

         
          dd = sqrt(sum((log(XS./X)).^2)) ;


end