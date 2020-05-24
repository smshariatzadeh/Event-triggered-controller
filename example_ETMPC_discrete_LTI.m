%% example_ETMPC_discrete_LTI.m
%
% Example of event triggered MPC 
%
% This example shows ET-MMPC control techniques over network controlled
% system for discerte-time linear system 
%
%% https://github.com/smshariatzadeh/Event-triggered-controller/blob/master/event-trigger-MPC-without-observer.png
%
% use MPT3
% 
%   Copyright 2019-2025 smshariatzadeh .


clc
clear
addpath('../src/')

%the usage is mostly same as tubeMPC
A = [1 1; 0 1];
B = [0.5; 1]; 
C = [ 1 1];
D= [0];
f=[0;0];
g=[0];
    
mysys = LTISystem('A', A, 'B', B, 'C', C, 'D', D, 'f', f, 'g', g);

%% make MPC controller
Q = diag([1, 1]);
R = 0.1;
Xc_vertex = [2, -2; 2 2; -10 2; -10 -2];
Uc_vertex = [1; -1];
Xc = Polyhedron(Xc_vertex);
Uc = Polyhedron(Uc_vertex);


%make mpc controller
mpc = ModelPredictiveControl(mysys, Q, R, Xc, Uc, 13);

% make event generator 
Sigma = 0.4;
et = EventGenerator( Sigma , mysys.nx);
x_init = [-7; -2];
Xnew = zeros(size(x_init));
Xold = zeros(size(x_init));


% This generic MPC doesn't guarantee the robustness. So if you
% add some noise (determined by w_min and w_max), and running the
% simulation multiple times, the state may hit the constraint Xc sometimes,
% which results in error, However, tube-MPC is robust to noise
% Please do some experiments:

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

    for i=1:dt:Tsimu
         t = t + dt;  
         t_array(i)=t;
         x = Xnew;
         if mod(i,1)==0 
             fprintf('\nRunning... Time: %f of %f',t , Tsimu);
         end    
         
         %% event generator part 
         %  get new system state , recognize event  and save event message in event_array for  using in the controller
         [Event, x_error, snormX, normX]  = et.Check_event( x);
         
         %save result
         event_array(i) = Event;  % message to controller         
         x_error_array(i) = x_error;
         snormX_array(i) = snormX;
         normX_array(i) = (norm(Xnew));
         
         %% simulation of the network delay
         % (apply delay to send Xdnew to controller inorder to perform calculation)
         if (i-hdelay)<=0
            Xdnew = x_init;
         else    
            Xdnew = x_array(:,i-hdelay);
         end          

         %% simulation of the controller part 
         % At the moment of event occurrence, this part receives Xdnew and calculates u for plant use
         
         if event_array(i) == 1  || i==1
             %event triggered
             event_time_array(i)= (i- lastevent)*dt ; %save event time for calculation of ineter event time
             lastevent = i;             
             fprintf(' event ') 
             
             % event occured so generate new u
             [x_nominal_seq, u_nominal_seq] = mpc.optcon.solve(Xdnew);
             u = u_nominal_seq(:, 1);    
             %save data for plot curve
             u_array(:, i)=u;
             uold = u; %save u for next step
         else
             % event not triggered so use old u
             u = uold; 
             
             %save data for plot curve
             u_array(:, i)=u;
         end         
         
         
         %% apply u to the system and find system response ( new x )
         w = rand(2, 1).*(mpc.w_max - mpc.w_min) + mpc.w_min;  % add small noise to system state
         Xnew = propagate(x, u, w); % calculate Xnew , immediately send data (Xnew) to observer without delay
         x_array (:, i+1) = Xnew; %save X for ploting result
         
         
         %% plot mpc trajectory
         titl = sprintf('Running... Time: %f of %f',t , Tsimu);
         clf;
         show_convex(mpc.Xc, 'm');
         show_trajectory(x_nominal_seq, 'gs-');
         show_trajectory(Xnew, 'b*-'), title(titl)
         pause(0.2)
    end

%% simulation result
figure(2);
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

 disp ('\nEvent Ratio:')
 NumberofEvent=sum(event_array);
 R=NumberofEvent/n
 %R_array(iter)=R
 %error_array(iter)=sum(abs(y_array));
  
