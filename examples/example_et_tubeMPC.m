%% example_ET_tubeMPC.m
%
% Example of event triggered Robust model predictive control of constrained discrete-time state-space system with bounded disturbance.
%
% This example shows Robust model predictive control techniques over network controlled
% system for discerte-time linear system 
%
%% https://github.com/smshariatzadeh/Event-triggered-controller/blob/master/fig/event-trigger-TubeMPC-without-observer.png
%
% use Matlab R2020a, MPT3
% 
% Copyright 2020-2024 smshariatzadeh .

% make your own discrete linear system
A = [1 1; 0 1];
B = [0.5; 1];
f = [0  ; 0];
C = [1  1];
D = [0];
g = [0];
    
mysys = LTISystem('A', A, 'B', B, 'C', C, 'D', D, 'f', f, 'g', g);% discrete-time state-space model


% constraints on state Xc and input Uc
%Xc_vertex = [2, -2;  2 2; -10 2; -10 -2];
Xc_vertex = [2 2; 2, -1.5; 1,-2;  -10 2; -10 -2];
Uc_vertex = [1; -1];
Xc = Polyhedron(Xc_vertex);
Uc = Polyhedron(Uc_vertex);


% construct a convex set of system noise (2dim here)
%W_vertex = [0.015, 0.015; 0.015, -0.015; -0.015, -0.015; -0.015, 0.015];  % low noise
W_vertex = [0.15, 0.15; 0.15, -0.15; -0.15, -0.15; -0.15, 0.15];  % high noise  
W = Polyhedron(W_vertex);

% actual noise
% Please note that rectangle defined by w_min and w_max must be included in W, otherwise, of course, the robustness is not guranteed.
w_min = [-0.1; -0.10];
w_max = [0.1; 0.10];

% create a tube_mpc simulater
% if N_step is too small and initial state is far from origin, the path will never reach inside the robust MRPI set (X_MPI - Z) in time step N_step,
% then the problem becomes infeasible.

N_step = 10;
Q = diag([1, 1]);
R = 0.1;
TubeMPC = TubeModelPredictiveControl(mysys, Q, R, Xc, Uc, W, N_step);

Tsimu = 30;
dt= 1; 
iprint = 5; % show detail output
exitflag = 1; % no warning flag in tube MPC 
t_Elapsed = 0; % No calculation in Tube MPC

 
Z = TubeMPC.Z; %minimal robust invariant set

%% make event generator 
Sigma = 0.54;
et = EventGenerator( Sigma , mysys.nx);

%%network delay
hdelay=0;

%% get initial state
if mysys.nx==2
    %% get initial state interactive
    figure(2)
    plot(TubeMPC.Xc, 'color', 'm');
    hold on
    plot(TubeMPC.Xc_robust, 'color', 'r');
    plot(TubeMPC.Xmpi, 'color',  [0.3, 0.3, 0.3]);
    plot(TubeMPC.Xmpi_robust, 'color', [0.5, 0.5, 0.5]); % gray
    hold off
    text(0,0,'O')
    title('Select a point and then press ENTER')
    hold off
    [x, y] = getpts
    x_init = [x(end), y(end)]';
else
    %x_init = [-7; -2];  % fix initial state
end

%% find nominal trajectory of TubeMPC
[x_nominal_seq, u_nominal_seq] = TubeMPC.optcon.solve(x_init);

u = zeros(mysys.nu,1);
n=round(Tsimu/dt);
event_array =  zeros(1,n);  
event_time_array = zeros(1,n); %save event time for calculation of sample time
r_array =  zeros(1,n);  
y_array =  zeros(1,n);
t_array =  zeros(1,n);
x_array =  zeros(mysys.nx,n);  
u_array =  zeros(mysys.nu,n);  
xhat_array =  zeros(mysys.nx,n);  
yhat_array =  zeros(1,n);
x_error_array =  zeros(1,n);
snormX_array =  zeros(1,n);
normX_array =  zeros(1,n);
lastevent=0;

propagate = @(x, u, w) mysys.A*x+mysys.B*u + w;
y = mysys.C*x;

x = x_init;
t0 = 0;
u_array(:, 1) = 0;
x_array(:,1) = x;
t_array(1) = t0;


msg ="";  %message show in output text

%% simulation loop
for i=1:Tsimu
    
    t = t0+ dt; 
    

     %% event generator part
     %  get new system state , recognize event  and save event message in event_array for  using in the controller
     [Event, x_error, snormX, normX]  = et.Check_event( x );

     %save result
     event_array(i) = Event;  % message to controller         
     x_error_array(i) = x_error;
     snormX_array(i) = snormX;
     normX_array(i) = (norm(x));

     %% simulation of the network delay
     % (apply delay to send Xdnew to controller inorder to perform calculation)
     if (i-hdelay)<=0
        Xdnew = x_init;
     else    
        Xdnew = x_array(:,i-hdelay);
     end          

     %% simulation of the controller part 
     % At the moment of event occurrence, this part receives Xdnew and calculates u for plant use

     if event_array(i) == 1  || i==1 || i==2
         %event triggered
         event_time_array(i)= (i- lastevent)*dt ; %save event time for calculation of ineter event time
         lastevent = i;             
         msg= ' event ';

        % The event has occured so generate new u
        % The robust MPC guidances the path inside the robust MaxPI-set so that the path will reach the robust MPI-set exactly at N_step.
        % After that (meaning that t > N_step), the system will be stabilized around the origin by just using LQR.

        if i<=TubeMPC.N
            u = u_nominal_seq(:, i) + TubeMPC.optcon.K*(Xdnew-x_nominal_seq(:, i));
            msg = msg + "dual control";
        else 
            u = TubeMPC.optcon.K*Xdnew;
            msg = msg + "state feedback control";
        end
        uold = u; %save u in the memory for next step if necessary
     else
         % Event has not triggered so use old u
         u = uold; 
         msg= 'no event ';
     end         

     %save data for plot curve
     u_array(:, i)=u;

    
    %save result for plot curve
    t0 = t;
    u_array(:, i) = u;
    t_array(i) = t;
    %% show this step resutl
    
    % Print result of this iteration
    if ( iprint >= 1 )
         Graphics.printSolution(i, iprint, exitflag, t_Elapsed, u_array, t_array, x_array, msg)
    end

    w = rand(2, 1).*(w_max  - w_min) + w_min;
    
    %calculate new X (Next time step)
    x = mysys.A*x + mysys.B*u + mysys.f + w;
    x_array(:,i+1) = x;
    
    clf; % real time plot
    Graphics.show_convex(TubeMPC.Xc, 'm');
    Graphics.show_convex(TubeMPC.Xc_robust, 'r');
    Graphics.show_convex(TubeMPC.Xmpi, [0.2, 0.2, 0.2]*1.5);
    Graphics.show_convex(TubeMPC.Xmpi_robust, [0.5, 0.5, 0.5]); % gray
    for j=1:TubeMPC.N+1
        Graphics.show_convex(x_nominal_seq(:, j)+TubeMPC.Z, 'g', 'FaceAlpha', 0.3);
    end
    Graphics.show_trajectory(x_nominal_seq, 'gs-');
    if i<=TubeMPC.N
         Graphics.show_trajectory(x, 'k*-');
    else
         Graphics.show_trajectory(x, 'b*-');        
    end
    pause(0.2)
end
x_array(:,i+1) = []; %remove extra state (x)

% time slice plot after simulation
figure(20)
Graphics.show_convex_timeslice(TubeMPC.Xc, -0.04, 'm');
Graphics.show_convex_timeslice(TubeMPC.Xc_robust, -0.03, 'r');
Graphics.show_convex_timeslice(TubeMPC.Xmpi, -0.02, [0.2, 0.2, 0.2]*1.5);
Graphics.show_convex_timeslice(TubeMPC.Xmpi_robust, -0.01, [0.5, 0.5, 0.5]);
Graphics.show_convex_timeslice(x_nominal_seq(:, 1)+TubeMPC.Z, TubeMPC.optcon.N, 'g', 'FaceAlpha', .3);
Graphics.show_trajectory_timeslice(x_nominal_seq, 'gs-', 'LineWidth', 1.2);
Graphics.show_trajectory_timeslice(x_array(:, 1:TubeMPC.N+1), 'b*-', 'LineWidth', 1.2);
leg = legend('$X_c$', '$X_c\ominus Z$', '$X_f (= X_{MPI})$', '$X_f\ominus Z$', 'Tube', 'Nominal', 'Real');
set(leg, 'Interpreter', 'latex')

for i=2:TubeMPC.N+1 % show remaining tubes.
        Graphics.show_convex_timeslice(x_nominal_seq(:, i)+TubeMPC.Z, TubeMPC.N-i+1, 'g', 'FaceAlpha', .3);
end

xlabel('x1');
ylabel('x2');
zlabel('time (minus)');
xlim([TubeMPC.optcon.x_min(1), TubeMPC.optcon.x_max(1)]);
ylim([TubeMPC.optcon.x_min(2), TubeMPC.optcon.x_max(2)]);
zlim([-0.05, TubeMPC.optcon.N]);
set( gca, 'FontSize',12); 
view([10, 10])
grid on;



%% simulation result
figure(21);
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


figure(22);
subplot(2,1,1)
stairs(t_array,x_array(1,:),'r');
xlabel('time(s)');ylabel('x1');
grid on
title('x 1')

subplot(2,1,2)
stairs(t_array,x_array(2,:),'r');
title('x 2')
xlabel('time(s)');ylabel('x2');
grid on


figure(3)
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


disp ('Event Ratio:')
NumberofEvent=sum(event_array);
R=NumberofEvent/n

