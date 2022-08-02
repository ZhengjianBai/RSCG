% function  RSCGPerry1

% This is a Riemannian Spectral Conjugate Gradient Method used to
% find the zeros of (X_1X_3, X_2X_3, X_3^2-1) on H^2.
% Xtrue: the true zero of the monotone tangent vector V   
% X: the computed zero of the monotone tangent vector V   

%% Parameter Settings ;
  
            ITmax = 1e3 ;     Tolerence = 1e-6 ;
           
              sigma = 1e-5 ;    dd = 1e-8 ;
               
                       rho = 0.5 ;

                               
%% Initial Guesses; 
                 
                    X = 0.5*rand(2,1) ;
                
               X = [ X ; sqrt(X'*X+1) ] ;

                    Xtrue = [ 0 0 1]' ; 
            
            
%% Begining;

                    t0 =clock ;           
                    
          %  F(X) = (-X_1X3, -X_2X_3, 1-X_3^2) ;    
               
                    FX = VV1(X) ;
                   
          %    The objective function  ;               
           
               beta0  = sqrt(LORENTZ(FX, FX)) ;

                     beta = beta0 ;
                    
                rdist0 = DIST1(Xtrue, X) ;
                    
                    
%% Initial Nonlinear Conjugate Direction;    
                    
                      XK = -FX ;
          
                       k = 0 ;   
                   
                     f_eval = 1 ;  
                               
                     
% Stopping Criterion ;
                                        
              while   beta > Tolerence && k < ITmax       
                  
                  
                        ak = 1 ;
   
                        
%% Intial Guess of The Step-size; 


                [ Z, FZ, lfro, rfro ] = PHI(X, XK, ak) ;        
             
                      f_eval = f_eval + 1 ;                                                      

                         l = 0 ;
                                       
                              while  lfro < sigma*rfro 
                                                               
                                         ak = rho*ak ;
                                          
                                 [ Z, FZ, lfro, rfro ] = PHI(X, XK, ak) ;
                         
                                       f_eval = f_eval + 1 ;
                                       
                                       l = l + 1 ;
                                                                    
                              end
                              
                              
%% Projection Step ;                         
                            
                       beta2 = LORENTZ(FZ,FZ)^(1/2) ; 
 
             if   beta2  > Tolerence

                                                        
                    FZ = FZ/sqrt(LORENTZ(FZ,FZ)) ;
            
                      aa = LORENTZ(FZ,FZ) ;
                      
                      bb = LORENTZ(X,FZ) ;
                       
                      PX = aa*X - bb*FZ ;
                      
                      if PX(3) > 0 ;
                          
                          mm = 1/sqrt(aa*(1+bb^2))  ;                 
                          
                      else
                          
                         mm = -1/sqrt(aa*(1+bb^2)) ;

                      end                    
                      
                          NX = mm*PX ;
        
                          NFX = VV1(NX) ;   
                         
                         
%%  Spectral Perry Strategy;


             %  Vector Transport ;   
       
                      TFX = TRANS(NX, FX) ; 

             %  The New Value ;   
                        
                   nbeta = LORENTZ(NFX,NFX)^(1/2) ;  
                                          
                        SK = -EXPINV1(NX,X) ;
                        
                   YK = NFX - TFX + nbeta*SK ;   
                   
                   
                      aa = LORENTZ(SK,YK) ;
                       
                    theta = LORENTZ(SK,SK)/aa  ;           
 
                bk = LORENTZ(theta*YK-SK,NFX)/aa ;
                
        

              %  The Spectral Conjugate Direction --- The Search Direction ;
              
                          XK = TRANS(NX, XK) ; 
              
                        XK = -theta*NFX + bk*XK ;   

                    if  LORENTZ(XK,NFX) > -dd*max(LORENTZ(XK,XK),nbeta^2)
   
                        XK = -theta*NFX  ;
        
                    end
                    
                    
                          beta = nbeta ;
                          
         fprintf('RSCG: Norm of (X_1X_3, X_2X_3, X_3^2-1) %d \n',nbeta)   

                             X = NX ;
                         
                            FX = NFX ;
                           
             else
                 
                          beta = beta2 ;
                          
         fprintf('RSCG: Norm of (X_1X_3, X_2X_3, X_3^2-1) %d \n',beta)  
                    
                          
                              X = Z ;
                          
                             FX = FZ ;
                 
                 
             end
                           
                           
                           k = k + 1 ;
                                       
                         
              end              
              
              
              
                       rdist = DIST1(Xtrue, X)   ;
                
                            iterk = k ;    
  
                     TimeCost = etime(clock,t0) ;

                            
          fprintf('\n');
          fprintf('RSCG: Number of Iterations %d \n', iterk)
          fprintf('RSCG: Number of function estimations %d \n',f_eval)  
          fprintf('RSCG: Initial norm of (X_1X_3, X_2X_3, X_3^2-1) ==========%d \n',beta0)
          fprintf('RSCG: Final norm of (X_1X_3, X_2X_3, X_3^2-1) ==========%d \n',beta)
          fprintf('RSCG: Initial Riemannian distance of X_0 and X_* ==========%d \n',rdist0)
          fprintf('RSCG: Final Riemannian distance of X_k and X_* ==========%d \n',rdist)
          fprintf('RSCG: Computing time used ========== %d \n', TimeCost)                
        
         Xtrue % the true zero 
         X % the computed zero 
         



