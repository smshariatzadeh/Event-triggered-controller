function example_et_nmpc_inverted_pend

% example of event triggered of nonlinear model predictive control for
% inverted pendulum

    addpath('./src');
    clear ;
    close all;

    global R           % set point
    global D           % Disturbance
    global landa       % ratio in cost function

    systemNx = 4;       % number of system state
    systemNu = 1;       % number of system input
    systemNy = 4;       % number of system output
    
    t0            = 0.0;
    x0            = [0.0 0.00 0.10  .10];
    
    xmeasure = x0;
    tmeasure = t0;
    
    Tsimu=5;
    t=0;
    dt=0.1;             % sample time   
    nmpcLastIterations = round(Tsimu/dt);    

    %% Nmpc controller design
    N             = 11;  % prediction horizon
    u0            = 0.2*ones(systemNu,N);    
    tol_opt       = 1e-8;
    opt_option    = 0; % Active-set method used for optimization (default)
    iprint        = 5; % Print closed loop data and errors and warnings of
                       % the optimization method as well as graphical
                       % output of closed loop state trajectories
    landa=.01;   % ratio in cost function
    nR=1;        % set point kind
    R=SetPoint(nmpcLastIterations,nR);
    
    nD=0;
    D=Disturbance(nmpcLastIterations,nD);
    
    systemtype          = 'difference equation'; %for discerete time system
    atol_ode_real = 1e-12;
    rtol_ode_real = 1e-12;
    atol_ode_sim  = 1e-4;
    rtol_ode_sim  = 1e-4;

    % make NMPC controller object
    nmpc = NonlinearModelPredictiveControl(@runningcosts, @terminalcosts, @constraints, ...
          @terminalconstraints, @linearconstraints, @system, systemtype,...
          N, dt, u0, ...
          tol_opt, opt_option, ...
          atol_ode_real, rtol_ode_real, atol_ode_sim, rtol_ode_sim);
      
    %% make event generator 
    Sigma = 0.05;
    et = EventGenerator( Sigma , systemNx);
    lastevent = 0;
    %% simulation  
    warning off all


    %% Start of the NMPC iteration
    for i = 1:nmpcLastIterations
        
         %init simulation loop       
         msg="";
         
         %% event generator part
         %  get new system state , recognize event  and save event message in event_array for  using in the controller
         [Event, x_error, snormX, normX]  = et.Check_event( xmeasure');
         
         %save result for plot curve
         x_error_array(i) = x_error;
         snormX_array(i) = snormX;
         normX_array(i) = (norm(xmeasure));
        
         event_array(i) = Event;  % message to controller         

         %controller calculation
         [t0, x0] = measureInitialValue ( tmeasure, xmeasure );

         if event_array(i) == 1  || i==1
             %event triggered
             event_time_array(i)= (i- lastevent)*dt ; %save event time for calculation of ineter event time
             lastevent = i;             
             msg = sprintf(' event ');
             
             % event occured so generate new u (run NMPC algorithm)
             % Solve the optimal control problem
             t_Start = tic;
             [u_new, V_current, exitflag, output] = nmpc.solve( t0, x0 );
             t_Elapsed = toc( t_Start );

             uold = u_new(:,1); %save u for next step             
             u = u_new(:,1);    %save u for apply to process
         else
             % event not triggered so use old u
             u = uold;              
        end
         
        %save result for plot curve
        u_array(:, i) = u;
        x_array(:,i) = xmeasure';
        t_array(i) = tmeasure;
        %y_array(i)=y; %save output
        %r_array(i)=0; %save setpoint
         
        %   Print result of this iteration
        if ( iprint >= 1 )
             Graphics.printSolution(i, iprint, exitflag, t_Elapsed, u_array, t_array, x_array, msg)
        end
        
        
        % Apply control to process and find cloed loop response
        if ( strcmp(systemtype, 'difference equation') )
            x1 = system(t0, x0, u, dt);
        elseif ( strcmp(systemtype, 'differential equation') )
            obj.options = odeset('AbsTol', atol_ode, 'RelTol', rtol_ode);
            [t_intermediate, x_intermediate] = ode45(system, ...
                [t0, dt], x0, options,  u);
            x1= x_intermediate;
        end
        xmeasure = x1;
        tmeasure = t0 + dt;
        
    end

     figure
     subplot(2, 1,1)
     plot(t_array,x_array, '-o'),title('x'),legend('x1','x2','x3','x4')
     subplot(2, 1,2)
     plot(t_array,u_array, '-ob'),title('u')

     
     
     
     %% simulation result
   
    figure(50);
    %for tracking
    %subplot(3,1,1)
    %plot(t_array,r_array,t_array,y_array);
    %xlabel('time(s)');ylabel('Output Value');
    %legend('y_{setpoint}','y')
    %grid on
    subplot(2,1,1)
    stairs(t_array,u_array);
    xlabel('time(s)');ylabel('Input Value');
    legend('u')
    grid on
    subplot(2,1,2)
    stem(t_array,event_array,'r','fill');
    title('Event-triggered')
    xlabel('time(s)');

    figure(60)
    clf
    subplot(3,1,1)
    stem(t_array,x_error_array,'fill'),hold on
    plot(t_array, snormX_array, '-', 'Marker' ,'.')
    plot(t_array, normX_array , '--', 'Marker' ,'.' ),hold off
    legend('x_error','sigma*normX','normX')
    xlabel('time(s)');    

    subplot(3,1,2)    
    stem(t_array,event_time_array)
    legend('Inter-event sampling times')
    xlabel('time(s)');
    grid on

    subplot(3,1,3)    
    stem(t_array,event_array,'r','fill')
    legend('event')
    xlabel('time(s)');
    grid on

     



end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the system dynamic and NMPC cost functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = runningcosts(t, x, u)
    global mpciter     % current time in mpc algorithm
    global R           % set point
    global D           % Disturbance
    global landa       % ratio in cost function

    cost = norm(x-[0, 0, 0, 0],2)^2+norm(u,2)^2;
  %  cost = norm(x-[0 0 R(mpciter+1)-D(mpciter+1)  0],2)^2+landa*norm(u,2)^2;
end

function cost = terminalcosts(t, x)
    cost = 0.0;
end

function [c,ceq] = constraints(t, x, u)
    c   = [];
    ceq = [];
end

function [c,ceq] = terminalconstraints(t, x)
    c   = [];
    ceq = [];
end

function [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, u)

    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    
    lb  = 0;
    ub  =  0.2;
    %ub  =  [];
    %lb  = [];
end

function yk = system(t, Xk, uk, Ts)
    %% Discrete-time nonlinear dynamic model of a pendulum on a cart at time k
    %
    % 1 inputs: (uk)
    % 4 outputs: (yk)
    %   same as states (i.e. all the states are measureable)
    %
    % xk1 is the states at time k+1.
    % Repeat application of Euler method sampled at Ts/M.
    xk= Xk';
    M = 5;
    delta = Ts/M;
    xk1 = xk;
    u = uk(1);
    for ct=1:M
        xk1 = xk1 + delta*pendulumCT1(xk1,u);
    end
    yk = xk1';
    % Note that we choose the Euler method (first oder Runge-Kutta method)
    % because it is more efficient for plant with non-stiff ODEs.  You can
    % choose other ODE solvers such as ode23, ode45 for better accuracy or
    % ode15s and ode23s for stiff ODEs.  Those solvers are available from
    % MATLAB.
end


 function [dxdt, y, A, B, C, D] = pendulumCT1(x, u)
    %System Parameters
    mp=0.127;
    Lr=0.216;
    Lp=0.337;
    Jr=0.00100;
    Br=0.00240;
    Jp=0.00120;
    Bp=0.00240;
    g=9.81;
    %------------------------------------------
    Rm=2.60;
    km=0.00767;
    kt=0.00767;
    Jm=4.60E-7;
    Dm=0.00;
    Jg=5.28E-5;
    JL=5.00E-7;
    DL=0.00150;
    kg=70;
    %--------------------------------------------------------------------
    % constant terms are defined
    a = mp*Lr^2+Jr;
    b = 0.25*mp*Lp^2;
    cs = 0.5*mp*Lp*Lr;
    % d = 0.5*mp*Lp^2;
    d = Jp+0.25*mp*Lp^2;
    f = 0.5*mp*Lp*g;
    % syms theta_dot alpha_dot alpha u

    % a = 0.0069;
    % b = 0.0036;
    % cs = 0.0046;
    % d = 0.0048;
    % f = 0.2099;

    %-------------------------------------------------------
    %% Obtain x, u and y
    % x
    % x=[0
    %     0
    %     0
    %     0];
    %   u=0;
    theta = x(1);
    alpha = x(2);
    theta_dot = x(3);
    alpha_dot = x(4);
    %  y = x;
    %% Compute dxdt
    % dxdt = x;
    % dxdt(1) = theta_dot;
    % dxdt(2) = alpha_dot;
    % % theta_dot
    % dxdt(3) = (d*u+f*cs*cos(alpha)*sin(alpha)-Br*theta_dot*d-2*b*d*sin(alpha)*cos(alpha)*theta_dot*alpha_dot-d*cs*sin(alpha)*alpha_dot^2-cs*cos(alpha)*Bp*alpha_dot+b*cs*cos(alpha)*cos(alpha)*sin(alpha)*theta_dot^2)/(d*(a+b-cos(alpha)^2)-cs^2*cos(alpha)*cos(alpha));
    % dxdt(4) = (cs*cos(alpha)*u-Br*cs*cos(alpha)*theta_dot-2*b*cs*cos(alpha)*cos(alpha)*sin(alpha)*theta_dot*alpha_dot-cs^2*cos(alpha)*sin(alpha)*alpha_dot^2-((a+b)-b*cos(alpha)^2)*Bp*alpha_dot+((a+b)-b*cos(alpha)^2)*b*cos(alpha)*sin(alpha)*theta_dot^2+((a+b)-b*cos(alpha)^2)*f*sin(alpha))/(d*(a+b-b*cos(alpha)^2)-cs^2*cos(alpha)*cos(alpha));
    %------------------------------------------------------------------------------------------------------------------
    %% Obtain A/B/C/D from Jacobian
    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ddx3dx2 = (d*(alpha_dot^2*cs*cos(alpha) + 2*alpha_dot*b*theta_dot*cos(alpha)^2 - 2*alpha_dot*b*theta_dot*sin(alpha)^2) + cs*sin(alpha)*(b*cos(alpha)*sin(alpha)*theta_dot^2 - Bp*alpha_dot + f*sin(alpha)) - cs*cos(alpha)*(f*cos(alpha) + b*theta_dot^2*cos(alpha)^2 - b*theta_dot^2*sin(alpha)^2))/(cs^2*cos(alpha)^2 - d*(a + b - b*cos(alpha)^2)) + ((d*(cs*sin(alpha)*alpha_dot^2 + 2*b*theta_dot*cos(alpha)*sin(alpha)*alpha_dot - u + Br*theta_dot) - cs*cos(alpha)*(b*cos(alpha)*sin(alpha)*theta_dot^2 - Bp*alpha_dot + f*sin(alpha)))*(2*cos(alpha)*sin(alpha)*cs^2 + 2*b*d*cos(alpha)*sin(alpha)))/(cs^2*cos(alpha)^2 - d*(- b*cos(alpha)^2 + a + b))^2;
    ddx3dx3 = (2*b*cs*theta_dot*cos(alpha)^2*sin(alpha))/(a*d + b*d - cs^2*cos(alpha)^2 - b*d*cos(alpha)^2) - (d*(Br + 2*alpha_dot*b*cos(alpha)*sin(alpha)))/(a*d + b*d - cs^2*cos(alpha)^2 - b*d*cos(alpha)^2);
    ddx3dx4 = (d*(2*alpha_dot*cs*sin(alpha) + 2*b*theta_dot*cos(alpha)*sin(alpha)) + Bp*cs*cos(alpha))/(cs^2*cos(alpha)^2 - d*(a + b - b*cos(alpha)^2));
    ddx4dx2 =((2*cos(alpha)*sin(alpha)*cs^2 + 2*b*d*cos(alpha)*sin(alpha))*(Bp*alpha_dot*(a + b - b*cos(alpha)^2) - f*sin(alpha)*(a + b - b*cos(alpha)^2) - cs*u*cos(alpha) + Br*cs*theta_dot*cos(alpha) + alpha_dot^2*cs^2*cos(alpha)*sin(alpha) - b*theta_dot^2*cos(alpha)*sin(alpha)*(a + b - b*cos(alpha)^2) + 2*alpha_dot*b*cs*theta_dot*cos(alpha)^2*sin(alpha)))/(cs^2*cos(alpha)^2 - d*(- b*cos(alpha)^2 + a + b))^2 - (alpha_dot^2*cs^2*sin(alpha)^2 - cs*u*sin(alpha) - alpha_dot^2*cs^2*cos(alpha)^2 + f*cos(alpha)*(a + b - b*cos(alpha)^2) + 2*b*f*cos(alpha)*sin(alpha)^2 + 2*b^2*theta_dot^2*cos(alpha)^2*sin(alpha)^2 + Br*cs*theta_dot*sin(alpha) + b*theta_dot^2*cos(alpha)^2*(a + b - b*cos(alpha)^2) - b*theta_dot^2*sin(alpha)^2*(a + b - b*cos(alpha)^2) - 2*alpha_dot*b*cs*theta_dot*cos(alpha)^3 - 2*Bp*alpha_dot*b*cos(alpha)*sin(alpha) + 4*alpha_dot*b*cs*theta_dot*cos(alpha)*sin(alpha)^2)/(cs^2*cos(alpha)^2 - d*(a + b - b*cos(alpha)^2));
    ddx4dx3 =(Br*cs*cos(alpha) + 2*alpha_dot*b*cs*cos(alpha)^2*sin(alpha) - 2*b*theta_dot*cos(alpha)*sin(alpha)*(a + b - b*cos(alpha)^2))/(cs^2*cos(alpha)^2 - d*(a + b - b*cos(alpha)^2));
    ddx4dx4 =(Bp*(a + b - b*cos(alpha)^2) + 2*alpha_dot*cs^2*cos(alpha)*sin(alpha) + 2*b*cs*theta_dot*cos(alpha)^2*sin(alpha))/(cs^2*cos(alpha)^2 - d*(a + b - b*cos(alpha)^2));
     %---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % used by A
    ddx3du = d/(a*d + b*d - cs^2*cos(alpha)^2 - b*d*cos(alpha)^2);
    ddx4du = -(cs*cos(alpha))/(cs^2*cos(alpha)^2 - d*(a + b - b*cos(alpha)^2));

    % % LTI
    A = [0 0 1 0;
         0 0 0 1;
         0 ddx3dx2 ddx3dx3 ddx3dx4;
         0 ddx4dx2 ddx4dx3 ddx4dx4];
    B = [0;0;ddx3du;ddx4du];
    % C = eye(4);
    % D = zeros(4,1);
    %Add actuator dynamics
    A(3,3) = A(3,3) - kg^2*kt*km/Rm*B(3);
    A(4,3) = A(4,3) - kg^2*kt*km/Rm*B(4);
     B=kg*kt*B/Rm;
    %-----------------------------------------------------------------------
     dxdt=A*x+B*u;
    %  s=ss(A,B,C,D)
    %  [V,D] = eig(A');
    %  pole(s)


 end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  other input of system (SetPoint and  Disturbance)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function     R=SetPoint(nmpcLastIterations,nR)

    R=zeros(nmpcLastIterations,1);
    if (nR == 0)  %% No set point
    end
    if (nR == 1)  %% for step 
        for i=1:nmpcLastIterations
            if  i>=5         
               R(i)=0.2;
            else
               R(i)=0;
            end    
        end
    end;
    if (nR == 2) %%% for pulse
        for i=1:nmpcLastIterations
            if mod(i,50)<=25  
               R(i)=0.2;
            else
               R(i)=-0.2;
            end    
        end
    end;
    if (nR == 3) %%% for pulse 
        for i=1:nmpcLastIterations
            if i>=25  
                if mod(i,50)<=25  
                   R(i)=0.2*i;
                else
                   R(i)=-0.2;
                end    
            end    
        end
    end
    if (nR == 4) %%% for pulse 
        for i=1:nmpcLastIterations
            if  i>=30          %% for step 
               R(i)=0.2*sin((i-30)/10);
            else
               R(i)=0;
            end    
        end
    end
end

function     D=Disturbance(nmpcLastIterations,nD)
% make Disturbance, Disturbance only affect on runningcost
% nD kind of Disturbance
    D=zeros(nmpcLastIterations,1);
    if (nD == 0)  %% No Disturbance
    end
    if (nD==1)
        for i=1:nmpcLastIterations
            if  i>=40          %% for step disturbance
               D(i)=1;
            else
               D(i)=0;
            end    
        end
    end
    if (nD==2)
        for i=1:nmpcLastIterations
            if mod(i,50)<=25  %% for pulse disturbance
               D(i)=1;
            else
               D(i)=0;
            end    
        end
    end
    
end


%   applyControl:                 applies the first control element of u to
%                                 the simulated process for one sampling
%                                 interval T.
%                                 The function returns closed loop state
%                                 vector xapplied at sampling instant
%                                 tapplied.
function [x1]= applyControl (t0, x0, u, dt)
% Apply control to process and find cloed loop response
        if ( strcmp(systemtype, 'difference equation') )
            x1 = system(t0, x0, u, dt);
        elseif ( strcmp(systemtype, 'differential equation') )
            obj.options = odeset('AbsTol', atol_ode, 'RelTol', rtol_ode);
            [t_intermediate, x_intermediate] = ode45(system, ...
                [t0, dt], x0, options,  u);
            x1= x_intermediate;
        end
end        


function [t0, x0] = measureInitialValue ( tmeasure, xmeasure )
%   measureInitialValue:          measures the new initial values for t0
%                                 and x0 by adopting values computed by
%                                 method applyControl.
%                                 The function returns new initial state
%                                 vector x0 at sampling instant t0.    
    t0 = tmeasure;
    x0 = xmeasure;
end




