%% example_ET_fo_MPC_discrete_LTI.m
%
% Example of event triggered of output feedback MPC 
%%
% This example shows ET-feedback output MPC control techniques over network controlled
% system for discerte-time linear system 
%
%% https://github.com/smshariatzadeh/Event-triggered-controller/blob/master/fig/event-trigger-MPC-with-observer.png
%
% use MPT3
% 
%   Copyright 2019-2025 smshariatzadeh .


clc
clear
addpath('../src/')

%% define system
%the usage is mostly same as tubeMPC
A = [1 1; 0 1];
B = [0.5; 1]; 
C = [ 1 1];
D= [0];
f=[0;0];
g=[0];
    
mysys = LTISystem('A', A, 'B', B, 'C', C, 'D', D, 'f', f, 'g', g);

%% make observer
%The command lqr can be adapted to calculate an optimal observer gain in a dual way: L = lqr(A',C',Q,R)'
Q = diag([1, 1]); 
R = 0.1;
L = dlqr(mysys.A', mysys.C', Q, R)';
ob = LTIObserver(mysys,L);

%% make MPC controller
Q = diag([1, 1]);
R = 0.1;
Xc_vertex = [4, -20; 4 15; -10 15; -10 -20];
Uc_vertex = [1; -2];
Xc = Polyhedron(Xc_vertex);
Uc = Polyhedron(Uc_vertex);

figure(1);
Graphics.show_convex(Xc, 'm');


%make mpc controller
mpc = ModelPredictiveControl(mysys, Q, R, Xc, Uc, 15);

%% make event generator 
Sigma = 0.4;
et = EventGenerator( Sigma , mysys.nx);

%% simulation
x_init = [-2.5; -2.55];
Xnew = zeros(size(x_init));
Xold = zeros(size(x_init));

%w_min = [0; -0.06];  % has error 
%w_max = [0; 0.05];

w_min = [0; 0];
w_max = [0; 0];

Tsimu=30;
t=0;
dt=1;
hdelay=0; %network delay , applicable to continous time system 

lastevent=0; %save the time of the last event for  calculation of inter event time
n=round(Tsimu/dt);
t_array = zeros(1,n);
u_array = zeros(mysys.nu,n);
uold = zeros(mysys.nu,1);
u = zeros(mysys.nu,1);
event_array =  zeros(1,n);  
event_time_array = zeros(1,n); %save event time for calculation of sample time
r_array =  zeros(1,n);  
y_array =  zeros(1,n);
x_array =  zeros(2,n);  
x_error_array =  zeros(1,n);
snormX_array =  zeros(1,n);
normX_array =  zeros(1,n);


x = x_init;
[x_nominal_seq] = mpc.optcon.solve(x);  % save nominal_seq
x_seq_real = [x];
x_array(: ,1)= x;
u_seq_real = [];
propagate = @(x, u, w) mysys.A*x+mysys.B*u + w;
Xnew = x;
y = mysys.C*x;

    for i=1:dt:Tsimu
         t = t + dt;  
         t_array(i)=t;
         if mod(i,1)==0 
             fprintf('\nRunning... Time: %f of %f',t , Tsimu);
         end    
         
         Xhat=  ob.EstimateCurrentState( y , u);
         Yhat = mysys.C*Xhat;
         xhat_array (:, i ) = Xhat; %save X for ploting result
         yhat_array (:, i ) = Yhat; %save Y for ploting result

         %% event generator part
         %  get new estimation of system state , recognize event  and save event message in event_array for sending to the controller
         [Event, x_error, snormX, normX]  = et.Check_event( Xhat);
         
         %save result
         event_array(i) = Event;  % message to controller         
         x_error_array(i) = x_error;
         snormX_array(i) = snormX;
         normX_array(i) = (norm(Xhat));
         
         %% controller calculation based on observer estimation
         if event_array(i) == 1  || i==1
             %event triggered
             event_time_array(i)= (i- lastevent)*dt ; %save event time for calculation of ineter event time
             lastevent = i;             
             fprintf(' controller gets event message') 
             
             % event occured so generate new u
             % At the moment of event occurrence, this part receives Xdnew and calculates u for plant use

             u = mpc.solve(Xhat);
             %save data for plot curve
             u_array(:, i)=u;
             
             %save data for plot curve             
             u_array(:, i)=u;
             uold = u; %save u for next step
         else
             % event not triggered so use old u
             u = uold; 
             
             %save data for plot curve
             u_array(:, i)=u;
         end         
              
         %% apply u to the system and find system response(  x+ )
                  
         w = rand(2, 1).*(mpc.w_max - mpc.w_min) + mpc.w_min;  % add small noise to system state
         Xnew = propagate(x, u, w); % calculate Xnew , immediately send data (Xnew) to observer without delay
         x = Xnew ; %save X for ploting result         
         Ynew = mysys.C*x;
         y=Ynew;
         x_array (:, i+1) = Xnew; %save X for ploting result
         y_array (:, i+1) = Ynew; %save X for ploting result
         
         %% plot mpc trajectory
         titl = sprintf('Running... Time: %f of %f',t , Tsimu);
         clf;
         Graphics.show_convex(mpc.Xc, 'm');
         Graphics.show_trajectory(x_nominal_seq, 'gs-');
         Graphics.show_trajectory(Xhat, 'b*-')
         Graphics.show_trajectory(x, 'bo-'), xlabel('x1'),ylabel('x2'), title(titl)
         legend('Xc','nominal sequence','Xhat(estimated state)' ,'X (real state)')
         pause(0.2)
         
    end
    
% remove extra state and output
x_array (:, i+1) = [];
y_array (:, i+1) = [];
    
%% simulation result 1
figure(2);
subplot(3,1,1)
stairs(t_array,x_array(1,:),'r');
xlabel('time(s)');ylabel('x1');
grid on
hold on
stairs(t_array,xhat_array(1,:),'g');
title('xhat1')
legend('x1','xhat1')
hold off

subplot(3,1,2)
stairs(t_array,x_array(2,:),'r');
hold on
title('xhat2')
xlabel('time(s)');ylabel('x2');
grid on
stairs(t_array,xhat_array(2,:),'g');
legend('x2','xhat2')
hold off

subplot(3,1,3)
stairs(t_array,u_array,'r');
hold on
title('u')
xlabel('time(s)');ylabel('u');
grid on
hold off



%% simulation results PART2
figure(3);
subplot(3,1,1)
stairs(t_array,x_array(1,:)- xhat_array(1,:) ,'r');
title('estimation error 1')
xlabel('time(s)');ylabel('estimation error x1');
grid on

subplot(3,1,2)
stairs(t_array,x_array(2,:) - xhat_array(2,:),'r');
title('estimation error 2')
xlabel('time(s)');ylabel('estimation error x2');

subplot(3,1,3)
stairs(t_array,y_array,'r');
title('y')
xlabel('time(s)');ylabel('y');
grid on

    
    
%% simulation results part3
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

figure(31)
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

fprintf('\n Event Ratio:') 
NumberofEvent=sum(event_array);
R=NumberofEvent/n  
