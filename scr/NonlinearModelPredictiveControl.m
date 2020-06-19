classdef NonlinearModelPredictiveControl < handle
%
% A class for the closed loop solution of the Nonlinear Model Predictive
% Control(NMPC) problem for discrete time system or continious time system. 
%
% https://github.com/smshariatzadeh/Event-triggered-controller   
%
% Copyright 2019-2025 smshariatzadeh.
% smshariatzadeh@yahoo.com
%
% class arguments:
%  
%   N:              Length of optimization horizon
%   T:              Sampling interval
%   u0:             Initial guess of open loop control
%
%   tol_opt:         Tolerance of the optimization method
%   opt_option: = 0: Active-set method used for optimization (default)
%                   = 1: Interior-point method used for optimization
%                   = 2: Trust-region reflective method used for
%                       optimization
%   type:               Type of dynamic, 'either difference equation' or
%                       'differential equation' can be used
%   atol_ode_real:      Absolute tolerance of the ODE solver for the
%                       simulated process
%   rtol_ode_real:      Relative tolerance of the ODE solver for the
%                       simulated process
%   atol_ode_sim:       Absolute tolerance of the ODE solver for the
%                       simulated NMPC prediction
%   rtol_ode_sim:       Relative tolerance of the ODE solver for the
%                       simulated NMPC prediction
%
% Internal Functions:
%   solve:                        solves the optimal control problem of the
%                                 horizon N with sampling length T for the
%                                 given initial values t0 and x0 and the
%                                 initial guess u0 using the specified
%                                 algorithm.
%                                 The function returns the computed optimal
%                                 control u, the corresponding value of the
%                                 cost function V as well as possible exit
%                                 flags and additional output of the
%                                 optimization method.
%   costfunction:                 evaluates the cost function of the
%                                 optimal control problem over the horizon
%                                 N with sampling time T for the current
%                                 data of the optimization method t0, x0
%                                 and u.
%                                 The function return the computed cost
%                                 function value.
%   nonlinearconstraints:         computes the value of the restrictions
%                                 for all sampling instances provided the
%                                 data t0, x0 and u given by the
%                                 optimization method.
%                                 The function returns the value of the
%                                 restrictions for all sampling instances
%                                 separated for inequality restrictions c
%                                 and equality restrictions ceq.
%   computeNstepSolution:         computes the system response over the
%                                 horizon N with sampling time T for the
%                                 initial values t0 and x0 as well as the
%                                 control u.
%                                 The function returns the complete open
%                                 loop solution over the requested horizon.
%   dynamic:                      evaluates the dynamic of the system for
%                                 given initial values t0 and x0 over the
%                                 interval [t0, tf] using the control u.
%                                 The function returns the state vector x
%                                 at time instant tf as well as an output
%                                 of all intermediate evaluated time
%                                 instances.
%
% nmpc problem defined by the following functions:
%
%   runningcosts:         evaluates the running costs for state and control
%                         at one sampling instant.
%                         The function returns the running costs for one
%                         sampling instant.
%          Usage: [cost] = runningcosts(t, x, u)
%                 with time t, state x and control u
%
%   terminalcosts:        evaluates the terminal costs for state at the end
%                         of the open loop horizon.
%                         The function returns value of the terminal costs.
%          Usage: cost = terminalcosts(t, x)
%                 with time t and state x
%
%   constraints:          computes the value of the restrictions for a
%                         sampling instance provided the data t, x and u
%                         given by the optimization method.
%                         The function returns the value of the
%                         restrictions for a sampling instance separated
%                         for inequality restrictions c and equality
%                         restrictions ceq.
%          Usage: [c,ceq] = constraints(t, x, u)
%                 with time t, state x and control u
%
%   terminalconstraints:  computes the value of the terminal restrictions
%                         provided the data t, x and u given by the
%                         optimization method.
%                         The function returns the value of the
%                         terminal restriction for inequality restrictions
%                         c and equality restrictions ceq.
%          Usage: [c,ceq] = terminalconstraints(t, x)
%                 with time t and state x
%
%   linearconstraints:    sets the linear constraints of the discretized
%                         optimal control problem. This is particularly
%                         useful to set control and state bounds.
%                         The function returns the required matrices for
%                         the linear inequality and equality constraints A
%                         and Aeq, the corresponding right hand sides b and
%                         beq as well as the lower and upper bound of the
%                         control.
%          Usage: [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, u)
%                 with time t, state x and control u
%
%   system:               evaluates the difference equation describing the
%                         process given time t, state vector x and control
%                         u.
%                         The function returns the state vector x at the
%                         next time instant.
%          Usage: [y] = system(t, x, u, T)
%                         with time t, state x, control u and sampling interval T
% for a given number of NMPC iteration steps. For
% the open loop problem, the horizon is defined by the number of
% time instances N and the sampling time T. Moreover, the
% initial time tmeasure, the state measurement xmeasure and a guess of
% the optimal control u0 are required.
%

% note:  Note that the dynamic can also be the solution of
% a differential equation or a difference equation.
%
% 
%
    properties (SetAccess = private)
        system % non linear system
        runningcosts;
        terminalcosts;
        constraints;
        terminalconstraints;
        linearconstraints;
        w_min; w_max; % lower and upper bound of system noise
                      % each vector has the same dim as that of system
        tol_opt;
        opt_option;
        options;
        type;   %system type  :'difference equation' or 'differential equation'
        atol_ode_real;
        rtol_ode_real;
        atol_ode_sim;
        rtol_ode_sim;
        N;
        T;
        u0; % vector of u0(Nu:N), initial guess for warm start
        x;  % vector of future x(Nx:N)
    end
    
    methods (Access = public)
        
        function obj = NonlinearModelPredictiveControl(runningcosts, terminalcosts, ...
              constraints, terminalconstraints, ...
              linearconstraints, system, systemtype,...
              N, T, u0, ...
              varargin)
          
                obj.N=N;
                obj.T=T;
                obj.u0=u0;
                obj.system =  system;
                if ( strcmp( systemtype, 'difference equation') || ...
                        strcmp( systemtype, 'differential equation') )
                    obj.type = systemtype;
                else
                    error([' Wrong input for type of systemtype: use either ', ...
                        '"difference equation" or "differential equation".']);
                end
                
                obj.runningcosts = runningcosts;
                obj.terminalcosts = terminalcosts;
                obj.constraints = constraints;
                obj.terminalconstraints = terminalconstraints;
                obj.linearconstraints = linearconstraints;
                
                if (nargin>=10)
                    obj.tol_opt = varargin{1};
                else
                    obj.tol_opt = 1e-6;
                end;
                if (nargin>=11)
                    obj.opt_option = varargin{2};
                else
                    obj.opt_option = 0;
                end;
                if (nargin>=12)
                    obj.atol_ode_real = varargin{3};
                else
                    obj.atol_ode_real = 1e-8;
                end;
                if (nargin>=13)
                    obj.rtol_ode_real = varargin{4};
                else
                    obj.rtol_ode_real = 1e-8;
                end;
                if (nargin>=14)
                    obj.atol_ode_sim = varargin{5};
                else
                    obj.atol_ode_sim = atol_ode_real;
                end;
                if (nargin>=15)
                    obj.rtol_ode_sim = varargin{6};
                else
                    obj.rtol_ode_sim = rtol_ode_real;
                end;

                % Determine MATLAB Version and
                % specify and configure optimization method
                vs = version('-release');
                vyear = str2num(vs(1:4));
                if (vyear <= 2007)
                    %MATLAB version R2007 or earlier detected
                    if ( obj.opt_option == 0 )
                        obj.options = optimset('Display','off',...
                            'TolFun', obj.tol_opt,...
                            'MaxIter', 2000,...
                            'LargeScale', 'off',...
                            'RelLineSrchBnd', [],...
                            'RelLineSrchBndDuration', 1);
                    elseif ( obj.opt_option == 1 )
                        error('nmpc:WrongArgument', '%s\n%s', ...
                              'Interior point method not supported in MATLAB R2007', ...
                              'Please use opt_option = 0 or opt_option = 2');
                    elseif ( obj.opt_option == 2 )
                         obj.options = optimset('Display','off',...
                             'TolFun', obj.tol_opt,...
                             'MaxIter', 2000,...
                             'LargeScale', 'on',...
                             'Hessian', 'off',...
                             'MaxPCGIter', max(1,floor(size(u0,1)*size(u0,2)/2)),...
                             'PrecondBandWidth', 0,...
                             'TolPCG', 1e-1);
                    end
                else
                    % MATLAB version R2008 or newer detected
                    if ( obj.opt_option == 0 )
                        obj.options = optimset('Display','off',...
                            'TolFun', obj.tol_opt,...
                            'MaxIter', 10000,...
                            'Algorithm', 'active-set',...
                            'FinDiffType', 'forward',...
                            'RelLineSrchBnd', [],...
                            'RelLineSrchBndDuration', 1,...
                            'TolConSQP', 1e-6);
                    elseif ( obj.opt_option == 1 )
                        obj.options = optimset('Display','off',...
                            'TolFun', obj.tol_opt,...
                            'MaxIter', 2000,...
                            'Algorithm', 'interior-point',...
                            'AlwaysHonorConstraints', 'bounds',...
                            'FinDiffType', 'forward',...
                            'HessFcn', [],...
                            'Hessian', 'bfgs',...
                            'HessMult', [],...
                            'InitBarrierParam', 0.1,...
                            'InitTrustRegionRadius', sqrt(size(u0,1)*size(u0,2)),...
                            'MaxProjCGIter', 2*size(u0,1)*size(u0,2),...
                            'ObjectiveLimit', -1e20,...
                            'ScaleProblem', 'obj-and-constr',...
                            'SubproblemAlgorithm', 'cg',...
                            'TolProjCG', 1e-2,...
                            'TolProjCGAbs', 1e-10);
                    %                       'UseParallel','always',...
                    elseif ( obj.opt_option == 2 )
                        obj.options = optimset('Display','off',...
                            'TolFun', obj.tol_opt,...
                            'MaxIter', 2000,...
                            'Algorithm', 'trust-region-reflective',...
                            'Hessian', 'off',...
                            'MaxPCGIter', max(1,floor(size(u0,1)*size(u0,2)/2)),...
                            'PrecondBandWidth', 0,...
                            'TolPCG', 1e-1);
                    end
                end

        end
        
        function x = computeNstepSolution(obj,t0, x0 , u)
            % compute Open loop Solution and save it in x
            x(1,:) = x0;
            for k=1:obj.N
                [xx, t_intermediate, x_intermediate] = obj.dynamic( t0 + k* obj.T, x(k,:), u(:,k));
                x(k+1,:)=xx;
            end
        end

        function [x, t_intermediate, x_intermediate] = dynamic(obj, t0, x0 , u)
            %Find system response 
            if ( strcmp(obj.type, 'difference equation') )
                x = obj.system(t0, x0, u, obj.T);
                x_intermediate = [x0; x];
                t_intermediate = [t0, t0+obj.T];
                 
            elseif ( strcmp(obj.type, 'differential equation') )
                obj.options = odeset('AbsTol', obj.atol_ode, 'RelTol', obj.rtol_ode);
                [t_intermediate,x_intermediate] = ode45(obj.system, ...
                    [t0, t0+obj.T], x0, obj.options, u);
                x = x_intermediate(size(x_intermediate,1),:);
            end
        end

        
        function cost = costfunction(obj, t0, x0, u )
            cost = 0;
            x = zeros(obj.N+1, length(x0));
            x = obj.computeNstepSolution( t0, x0 , u);
            for k=1:obj.N
                cost = cost + obj.runningcosts(t0+k*obj.T, x(k,:), u(:,k));
            end
            cost = cost + obj.terminalcosts(t0+(obj.N+1)*obj.T, x(obj.N+1,:));
        end

        function [c,ceq] = nonlinearconstraints(obj, t0, x0, u )
            x = zeros(obj.N+1, length(x0));
            x = obj.computeNstepSolution( t0, x0, u );
            c = [];
            ceq = [];
            for k=1:obj.N
                [cnew, ceqnew] = obj. constraints(t0+k*obj.T,x(k,:),u(:,k));
                c = [c cnew];
                ceq = [ceq ceqnew];
            end
            [cnew, ceqnew] = obj.terminalconstraints(t0+(obj.N+1)*obj.T,x(obj.N+1,:));
            c = [c cnew];
            ceq = [ceq ceqnew];
        end


        function [u, V, exitflag, output] = solve (obj,  t0, x0)
            
                obj.x = zeros(obj.N+1, length(x0));                
                % Step (1) of the NMPC algorithm:
                %   Obtain new initial value
                
                obj.x =obj.computeNstepSolution( t0, x0 , obj.u0);

                % Step (2) of the NMPC algorithm:                
                % Set control and linear bounds
                A = [];
                b = [];
                Aeq = [];
                beq = [];
                lb = [];
                ub = [];
                for k=1:obj.N
                    [Anew, bnew, Aeqnew, beqnew, lbnew, ubnew] = ...
                           obj.linearconstraints(t0+k*obj.T,obj.x(k,:),obj.u0(:,k));
                    A = blkdiag(A,Anew);
                    b = [b, bnew];
                    Aeq = blkdiag(Aeq,Aeqnew);
                    beq = [beq, beqnew];
                    lb = [lb, lbnew];
                    ub = [ub, ubnew];
                end

                % Step (3) of the NMPC algorithm:                
                % Solve optimization problem by Matlab solver
                % use u0 for initializing u
                
                [u, V, exitflag, output] = fmincon(...
                    @(u) obj.costfunction( t0, x0, u ),...
                    obj.u0,...
                    A, b, Aeq, beq, lb, ub,...
                    @(u) obj.nonlinearconstraints( t0, x0, u ),...
                    obj.options);
                
                % Step (4) of the NMPC algorithm:                
                %   Prepare initial guess to warm restart for next step
                obj.u0 = obj.shiftHorizon(u); %initial value of next step
        
        end

        

        function u0 = shiftHorizon(obj, u)
        %   shift horizon data:           applies the shift method to the open loop
        %                                 control (u) in order to ease the restart.
        %                                 The function returns a new initial guess
        %                                 u0 of the control for warm start.
        u0 = [u(:,2:end) u(:,end)];
        obj.u0 = u0;
        end


        function [] = simulate(obj, Tsimu, x_init)
%             x = x_init;
%             [x_nominal_seq] = obj.optcon.solve(x, obj.u0);
%             x_seq_real = [x];
%             u_seq_real = [];
%             
%             for i=1:Tsimu
%                 [x_nominal_seq, u_nominal_seq] = obj.optcon.solve(x);
%                 u = u_nominal_seq(:, 1);
%                 w = rand(2, 1).*(obj.w_max - obj.w_min)+obj.w_min;
%                 x = dynamic(x, u, w);
%                 x_seq_real = [x_seq_real, x];
%                 u_seq_real = [u_seq_real, u];
% 
%             end
        end
    end
    
end
