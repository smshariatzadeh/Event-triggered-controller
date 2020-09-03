classdef LTIObserver < handle
    % EventGenerator , a simple class for observer for discerete time system
    % 
    % https://github.com/smshariatzadeh/Event-triggered-controller   
    %
    % Copyright 2019-2025 smshariatzadeh.
    % smshariatzadeh@yahoo.com
    %
    properties (SetAccess = private)
        % common notation for linear observer 
        
        Obsys   % LTI system
        L;      % observer Gain
        Al;     %
        Ld;     % observer gain for disturbance estimation
        Xhat;   % 
        Yhat;   % observer output
        Dhat ;  % estimation of disturbance
        Uold;
        Yold;
    end

    methods (Access = public)

        function obj = LTIObserver(sys, L )
            % LTIObserver Construct an instance of this class
 
            %check dimension of L and sys
            [na,ma]=size(sys.A);
            [nb,mb]=size(sys.B);
            [nc,mc]=size(sys.C);            
            
            if na~=ma %Verify A be square
                error('A must be square')
            else
                [n, n] =size (sys.A) ;
            end

            [nl,ml]=size(L);
            if (n ~= nl)  
               %Check that A and L are the same size
               error('A and L must be the same size');
            end
            if (mb ~= ml)
               %Check that A and L are the same size
               error('B and L must be the same column number');
            end

            %Data Checks observability for system setup
            % if n > rank(obsv(sys.A,sys.B))
            %     %Check for observability
            %     error('System Not observable')
            %     return
            % end

            %check observer stability
            if any(abs(eig(sys.A-L*sys.C))>=1)
                error('ERROR: A-L*C is not Schur, A is Unstable (eig out of unit circle)')
            end            
            obj.Obsys = sys;
            
            obj.L = L;
            obj.Ld = 0.01;  %disturbance estimator gain for siso system
            
            obj.Al = (sys.A + L*sys.C );
            
            %initialize Xhat
            obj.Xhat = zeros(na, 1 ); % for generate event in the first time step
            obj.Uold = zeros(mb, 1 ); % for generate event in the first time step
            obj.Yold = zeros(nc, 1 ); % for generate event in the first time step
            
        end
        
        function []= InitObserver(obj, X0hat )
            %function for initialize Xhat
            
            %check dimension of Xhat and sys
            na=size(obj.Obsys.A,1) ;
            if size( X0hat ,1 ) == na
               %initialize Xhat
               obj.Xhat = X0hat ; % for generate event in the first time step
% %             obj.Uold = zeros(mb, 1 ); % for generate event in the first time step
% %             obj.Yold = zeros(nc, 1 ); % for generate event in the first time step
            else
               error('ERROR: Xhat size is not consistante with the system state size')
            end    
        end
        
        function [Xhatkp1]  = EstimateCurrentState(obj, Ykp1 , Uk)
            % Read current output & previous U and calculate current X
            %
            % X[k+1]= Ob (Y[k+1] , BU[k] )
            % X[k+1]= AX[k]+BU[t]- L{Y[k+1]- C(AX[k]+BU[k]) }
            %
            sy = size(Ykp1);
            su = size(Uk);
            if sy(1) ~= obj.Obsys.ny
               error ('New output is not consistante with the system state space model') 
            end    
            if su(1) ~= obj.Obsys.nu
               error ('New input is not consistante with the system state space model') 
            end    
             
            % calculation 
            dh_i=0;
            xh_ipls =  obj.Obsys.A*obj.Xhat + obj.Obsys.B*Uk + obj.L*( Ykp1 - obj.Obsys.C*((obj.Obsys.A*obj.Xhat + obj.Obsys.B*Uk ) - dh_i )) ;

            dh_ipls = dh_i + (obj.Ld*(Ykp1 - obj.Obsys.C*obj.Xhat - dh_i));


            % old data is output for this step
            obj.Xhat = xh_ipls;
            Xhatkp1 = xh_ipls; %output
            obj.Dhat=dh_ipls; % current estimated disturbance (Output argument)
             
        end
        
        function [Xhatkp1]  = EstimateNextState(obj, Yk , Uk)
            % Get current output &  U and calculate next X
            %
            % X[k+1]= Ob (Y[k] , BU[k] )
            % X[k+1]= AX[k]+BU[t]- L{Y[k]- CX[k] }
            %
            sy = size(Ykp1);
            su = size(Uk);
            if sy(1) ~= obj.Obsys.ny
               error ('New output is not consistante with the system state space model') 
            end    
            if su(1) ~= obj.Obsys.nu
               error ('New input is not consistante with the system state space model') 
            end    
             
            % calculation 
            dh_i = 0;
            obj.Yold = Yk;
            obj.Uold = Uk;
            xh_ipls =  (obj.Obsys.A*obj.Xhat + obj.Obsys.B*obj.Uold + obj.L*(obj.Yold - obj.Obsys.C*obj.Xhat - dh_i )) ;
            Xhatkp1 = xh_ipls;
            %dh_ipls = dh_i + (Ld*(yi - Cd1*xh_i - dh_i));

            % old data is output for this step
            obj.Xhat = xh_ipls;
            %obj.Dhat=dh_ipls; % current estimated disturbance (Output argument)

            %xh_i=xh_ipls; 
            % saved data for new 
             
        end
        
        function Z_approx = ApproxmRPIset(obj,  W, n_order, alpha)
            % compute approximate minimal robust positively invariant set
            % for observer
            % which takes the form of Z = alpha*(W + Al*W + Al^2*W + ... Ak^n_ordr*W).
            % where + denotes Minkowski addition.
            %
            %input:
            %   Al: closed loop system matrix 
            %   W : system noise (as Polyhedron ) (from obj)
            %   n : degree of approximate (n >=1)
            %       choose bigger n for best approximation 
            %   alpha :a magnifier constant  (alpha>=1)
            %
            %output:
            %  Zapprox: approximation of minimal robust positively invarient set as
            %           polyhedron
            %
            % We could obtain approximate minimal RPI (Z) by computing an infinite geometric series,
            % which is not practicall to get.
            % we approximate this by trancating the polynomial.
            
            if (n_order <1)
               error ('degree of approximate must be (n >=1)')
            end
            if (alpha<1)
               error ('magnifier constant must be (alpha>=1)') 
            end
            Z_approx = W;
            
            for n = 1:n_order
                Z_approx = Z_approx + obj.Al^n*W;
            end
            Z_approx = Z_approx*alpha;
            % which takes the form of Z = alpha*(W + Al*W + Al^2*W + ... Al^n_ordr*W).
            % where + denotes Minkowski addition.
        end
        
        
    end
end

