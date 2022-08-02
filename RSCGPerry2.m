% function  RSCGPerry2

% This is a Riemannian SpectralConjugate Gradient Method used to
% find the zeros of F(X):
% X: the computed zero of the monotone tangent vector V   
% Xtrue=X: We see the computed zero as the true zero of the monotone tangent vector V   

%% Parameter Settings ;
  
            ITmax = 1e3 ;     Tolerence = 1e-6 ;
           
              sigma = 1e-5 ;    dd = 1e-8 ;
               
                       rho = 0.5 ;
                       
                        N = 200 ;

                               
%% Initial Guesses; 
                 
                    X = 0.5*rand(N,1) ;
                    
                       X0 = X ;
            
%% Begining;

                    t0 =clock ;           
                    
          %  F(X) ;    
               
                    FX = VV2(X) ;
                   
          %    The objective function  ;               
           
               beta0  = sqrt(INPRODUCT(FX, FX, X)) ;

                     beta = beta0 ;
                    
                    
%% Initial Nonlinear Conjugate Direction;    
                    
                      XK = -FX ;
          
                       k = 0 ;   
                   
                     f_eval = 1 ;  
                               
                     
% Stopping Criterion ;
                                        
            while   beta > Tolerence && k < ITmax       
                  
                  
                        ak = 1 ;
                         
%% Intial Guess of The Step-size; 


                   rfro = INPRODUCT(XK, XK, X) ;
                   
                       Z = EXPP2(X, XK, ak) ;
                      
                           FZ = VV2(Z) ;
                         
                    lfro  = INPRODUCT(FZ, XK, Z) ;
             
                      f_eval = f_eval + 1 ;                                                  

                         l = 0 ;
                                       
                              while  -lfro < sigma*ak*rfro ; 
                                                               
                                         ak = rho*ak ;
                                          
                                       Z = EXPP2(X, XK, ak) ;
                      
                                           FZ = VV2(Z) ;
                          
                                  lfro  = INPRODUCT(FZ, XK, Z) ;
                         
                                        f_eval = f_eval + 1 ;
                                       
                                       l = l + 1 ;
                                                                    
                              end
                              
                              
%% Projection Step ;                         
                            
                       beta2 = INPRODUCT(FZ,FZ,Z)^(1/2) ;
 
             if   beta2  > Tolerence

                                                        
                               NX = X ;
                 
                          IND1 = find(FZ<0) ;
                                                
                      IND2 = find(X(IND1) < Z(IND1)) ;
                      
                           NX(IND2) = Z(IND2) ;
                           
                           IND3 = find(FZ>0) ;
                           
                       IND4 = find(X(IND3) > Z(IND3)) ;
                       
                           NX(IND4) = Z(IND4) ;
                         
                              NFX = VV2(NX) ;  
                         
%%  Spectral Perry Strategy;
 
             %  The New Value ;   
                        
                  nbeta = sqrt(INPRODUCT(NFX, NFX, NX)) ;  
                                          
                        SK = -EXPINV2(NX,X) ;
                        
                   YK = NFX - FX + nbeta*SK ;   
                   
                   
                        aa = INPRODUCT(SK,YK, NX) ;
                       
                     theta = INPRODUCT(SK,SK,NX)/aa  ;   

                  bk =  INPRODUCT(theta*YK-SK,NFX,NX)/aa ;
                
        

              %  The Spectral Conjugate Direction --- The Search Direction ;
 
              
                        XK = -theta*NFX + bk*XK ;   

                    if  INPRODUCT(XK,NFX,NX) > -dd*max(INPRODUCT(XK,XK,NX),nbeta^2)
   
                           XK = -theta*NFX  ;
        
                    end
            
                          beta = nbeta ;
                          
         fprintf('RSCG: Norm of F(X) %d \n',nbeta)   

                             X = NX ;
                         
                            FX = NFX ;
                           
             else
                 
                          beta = beta2 ;
                          
         fprintf('RSCG: Norm of F(X) %d \n',beta)  

                              X = Z ;
                          
                             FX = FZ ;
                 
                 
             end
                           
                           
                           k = k + 1 ;
                                       
                         
            end              
              
                            iterk = k ;    
  
                     TimeCost = etime(clock,t0) ;
                     Xtrue = X;
                       rdist0 = DIST2(X0,Xtrue) ;
                       
                       rdist = DIST2(X,Xtrue) ;
 
        
                            
          fprintf('\n');
          fprintf('RSCG: Number of Iterations %d \n', iterk)
          fprintf('RSCG: Number of function estimations %d \n',f_eval)  
          fprintf('RSCG: Initial norm of F(X) ==========%d \n',beta0)
          fprintf('RSCG: Final norm of F(X) ==========%d \n',beta)
          fprintf('RSCG: Initial Riemannian distance of X_0 and X_* ==========%d \n',rdist0)
          fprintf('RSCG: Final Riemannian distance of X_k and X_* ==========%d \n',rdist)
          fprintf('RSCG: Computing time used ========== %d \n', TimeCost)    
        
         Xtrue; % the true zero 
         X; % the computed zero 
                          
         

% end

