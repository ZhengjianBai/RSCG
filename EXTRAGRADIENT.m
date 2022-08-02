function  EXTRAGRADIENT

% This is a Riemannian extragradient type method used to
% find the zeros of F(X).

%% Parameter Settings ;
  
            ITmax = 1e3 ;     Tolerence = 1e-7 ;
           
              sigma = 1e-5 ;      rho = 0.5 ;
              
                         N = 800 ;
                     
                               
%% Initial Guesses; 
                 
                    X = 0.5*rand(N,1) ;

%% Begining;

                    t0 =clock ;           
                    
          %  F(X) ;    
               
                     FX = FF(X) ;
                   
          %    The objective function  ;               
           
               beta0  = sqrt(INPRODUCT(FX, FX, X)) ;

                        beta = beta0 ;    

%% Extragradient Method;

                   
                       XK = -FX ;

                         k = 0 ;   
                   
                       f_eval = 1 ;            
                     
% Stopping Criterion ;
                                        
             while   beta > Tolerence && k < ITmax                                             
                  
                        ak = 1 ;
   
%% Intial Guess of The Step-size; 

                   rfro = INPRODUCT(XK, XK, X) ;
                   
                       Z = EXPP(X, XK, ak) ;
                      
                           FZ = FF(Z) ;
                         
                    lfro  = INPRODUCT(FZ, XK, Z) ;
             
                      f_eval = f_eval + 1 ;                                                      

                         l = 0 ;
                                       
                              while  -lfro < sigma*ak*rfro ; 
                                                               
                                         ak = rho*ak ;
                                          
                                       Z = EXPP(X, XK, ak) ;
                      
                                           FZ = FF(Z) ;
                          
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
                         
                              FX = FF(NX) ;  
           
                               XK = -FX ;

                                X = NX ;
                          
          %    The objective function  ;               
           
                 beta  =  sqrt(INPRODUCT(FX, FX, X)) ;
                 
         fprintf('ETRGRAD: Norm of F(X) %d \n',beta)   
                                                       
              else
                 
                           beta = beta2 ;
                          
         fprintf('ETRGRAD: Norm of F(X) %d \n',beta)  
              
                               X = Z ;
         
              end

                            k = k + 1 ;               
 
             end              
     
              
                            iterk = k ;    
  
                     TimeCost = etime(clock,t0) ;

                            
          fprintf('\n');
          fprintf('ETRGRAD: Number of Iterations %d \n', iterk)
          fprintf('ETRGRAD: Number of function estimations %d \n',f_eval)  
          fprintf('ETRGRAD: Initial norm of F(X) ==========%d \n',beta0)
          fprintf('ETRGRAD: Final norm of F(X) ==========%d \n',beta)
          fprintf('ETRGRAD: Computing time used ========== %d \n', TimeCost)                
                          
         

end

