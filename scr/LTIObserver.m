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
        
        Obsys   % Lti system
        L;      % Gain
        Ld;     % observer gain or disturbance estimation
        Xhat;   % 
        Yhat;   %
        Dhat ;   % disturbance estimation 
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
            
            %initialize Xhat
            obj.Xhat = zeros(na, 1 ); % for generate event in the first time step
            obj.Uold = zeros(mb, 1 ); % for generate event in the first time step
            obj.Yold = zeros(nc, 1 ); % for generate event in the first time step
            
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
            dh_i=0;
            xh_ipls =  (obj.Obsys.A*obj.Xhat + obj.Obsys.B*obj.Uk + obj.L*(obj.Yk - obj.Obsys.C*obj.Xhat - dh_i )) ;
            Xhatkp1 = xh_ipls;
            %dh_ipls = dh_i + (Ld*(yi - Cd1*xh_i - dh_i));

            % old data is output for this step
            obj.Xhat = xh_ipls;
            %obj.Dhat=dh_ipls; % current estimated disturbance (Output argument)

            %xh_i=xh_ipls; 
            % saved data for new 
             
        end
        
    end
end

