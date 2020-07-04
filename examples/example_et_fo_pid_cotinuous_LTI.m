% example_et_pid_cotinuous_LTI.m
%
% PID(output feedback) event trigger for continuous time linear time-invariant system 
%
% simulation is run for different event trigger parameters i.e.: in sigma(event-triggered parameter)  
%
% event triggered rule is based on comparing norm Y with norm e
%
% sigma : event-triggered parameter in event rule 
%
% 
% System Configuration figure:
% https://github.com/smshariatzadeh/Event-triggered-controller/blob/master/fig/event-trigger-control-outputfeedback-PID.png
%         
% By S.M.Shariatzadeh
% Date :25 Mars 2020 

clear
%define linear time-invarient system (open loop unstable system)
 
 sys.A= [0.1  0;
          0   -0.1 ];
 sys.B= [.2;
         2];
 sys.C= [1  1];
 
 na=size(sys.A);
 nb=size(sys.B);
 nc=size(sys.C);

 %% calculate system data & Find Dimensions of A, B , Q, R
%
[na,ma]=size(sys.A);
[nb,mb]=size(sys.B);
[nc,mc]=size(sys.C);
 
%Nbar =  sys.C*(inv(sys.A - sys.B * K))*sys.B ; 

% designing PID gain 
 ki = 0.5;     % integral action
 ti = 1;
 kp = 2.5;
 kd = 0;   % PID action

%% network parameters
hdelay=3; %network delay

%% run simulation for different event trigger rule
sigma_array=linspace(0,0.4,10); % 0<sigma<1
[p,totalRun]=size(sigma_array);
Tsimu=5;
dt=0.005; % minimum inter-event time = very small sample time that system is linear for it

for iter = 1:totalRun
    Sigma = sigma_array(iter);
    fprintf('step %d of %d: simulation for sigma %d:\n', iter, totalRun ,Sigma);
    
    %initialize variable for each run
    x0=[1 2 3]'; %initial state 
    x0=[1 3]'; %initial state 
    y0 = sys.C*x0;    
    x=x0;
    t=0.0;
    s=0;
    y=y0;
    yb=0.0;
    r=0.0; % set point    
    u=zeros(mb,1);
    uold=zeros(mb,1);
    dist=0; % disturbance
    lastevent=0; %save the time of the last event
    n=round(Tsimu/dt);
    t_array = zeros(1,n);
    u_array = zeros(mb,n);
    event_array =  zeros(1,n);
    event_time_array = zeros(1,n); %save event time for calculation of sample time    
    r_array =  zeros(nc,n);  
    y_array =  zeros(nc,n);
    x_array =  zeros(na,n);  
    y_error_array=  zeros(1,n);
    normY_array=  zeros(1,n);
    snormY_array=  zeros(1,n);
    Ynew=zeros(nc,1);
    Ynew=y0;
    Yold=zeros(nc,1);

    % start of simulation loop 
    for i=1:n
         t = t + dt;       
         if mod(i,100)==0 
             fprintf('\nRunning... Time: %f',t) 
         end    
         %only output is available to controller
         

         %% event generator part (only recognize event and save event message in event_array for  using in the controller
         y_error = Ynew-Yold;  
         y_error_array(i) = abs(norm(y_error));
         snormY_array(i)= Sigma*(norm(Ynew));
         normY_array(i)= (norm(Ynew));
         if (Sigma*(norm(Ynew))<= abs(norm(y_error)))
             % event 
             event_array(i)=1;
             event_time_array(i)= (i- lastevent)*dt ; %save event time for calculation of sample time
             lastevent=i; 
             Yold=Ynew;                          
         else
             % event not triggered 
             event_array(i)=0;
         end

         %% simulation of the network delay
         if (i-hdelay)<=0
            Ydnew = y0;
         else    
            Ydnew = y_array(:,i-hdelay);
         end          
                 
         %% controller part
         if event_array(i)==1
             % there is an event , so get y and generate new u(PID)
             
             e = (r-Ydnew);
             s = s + (ki/ti)*(e) * dt;     % integral summation
             %calculate U 
             u= kp*(e)+ s +kd*(e)/dt; % PID action
            
             %save data for plot curve
             u_array(:, i)=u;
             uold = u; %save u for next step
         else
             % event not triggered so use old u
             u = uold; 
             
             %save data for plot curve
             u_array(:, i)=u;
         end
         
         
         %% find system response by applling u and calculate new x by integration method
         xdot = sys.A*x + sys.B*u;
         x = x + xdot*dt;
         y = sys.C*x;
         Ynew=y; %save output in sensor memory and event generator part and controller
         
         %save result
         x_array(:,i) = x;
         t_array(i) = t;
         y_array(i) = y;
         r_array(i) = r;
         
    end

    figure(1);
    subplot(3,1,1)
    plot(t_array,r_array,t_array,y_array);
    xlabel('time(s)');ylabel('Output Value');
    legend('y_{setpoint}','y')
    grid on
    subplot(3,1,2)
    plot(t_array,u_array);
    xlabel('time(s)');ylabel('Input Value');
    legend('u')
    grid on
    subplot(3,1,3)
    stem(t_array,event_array);
    title('Event-triggered')
    xlabel('time(s)');
    
    figure(2)
    subplot(3,1,1)
    plot(t_array,y_error_array,t_array, snormY_array, t_array, normY_array  );
    xlabel('time(s)');    
    legend('x_error','sigma*normX','normX')
    
    subplot(3,1,2)    
    stem(t_array,event_time_array)
    legend('Inter-event sampling times')
    xlabel('time(s)');
    grid on
    
    subplot(3,1,3)    
    stem(t_array,event_array)
    xlabel('time(s)');    
    legend('event')
    grid on
    
        

    NumberofEvent=sum(event_array);
    R=NumberofEvent/n;
    R_array(iter)=R;
    error_array(iter)=sum(abs(y_array));
    
    fprintf('\n\nR(Event Ratio)= %6.3f\n', R)
    disp('------------------------')
    disp('press any key to run another simulation')
    pause
end

figure(5)
subplot(3,1,1)
set(gca,'FontSize',10);
bar(sigma_array,R_array,'b');
grid on
xlabel('Sigma'); ylabel('R'); title('Event Ratio (R) vs Error Criterion (Sigma)');

subplot(3,1,2)
bar(sigma_array,error_array,'r');
grid on
xlabel('Sigma'); ylabel('abs(error)'); title(' Output Error vs Error Criterion (Sigma)');
print('sigma_vs_R_and_E','-dpng');

subplot(3,1,3)
bar(sigma_array,R_array*n,'g');
grid on
xlabel('Sigma'); ylabel(''); title({'Number of sampled data sent from plant to controller','vs Error Criterion (Sigma)'} );
print('sigma_vs_R_and_E3','-dpng');

disp('end of simulation')
