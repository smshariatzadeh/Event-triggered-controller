% example_etst_cotinuous_LTI.m
%
% static event trigger for continuous time linear time-invariant system
%
% simulation is run for different event trigger parameters i.e.: in C0
% (static error threshold in event-triggered rule)  
%
% event triggered rule is based on comparing norm e with C0( constant 
% function)
%
% C0 : event-triggered parameter in event rule 
%         
% By S.Majid.Shariatzadeh
% Date :12 Mars 2020 


%define linear time-invarient system
sys.A= [0 1 0;
    4 5 0;
    1 0 1 ];
sys.B= [0;
    1;
    2];
sys.C= [1 1 1];

na=size(sys.A);
nb=size(sys.B);
nc=size(sys.C);

% designing feedback gain based on Linear-Quadratic Regulator (LQR)  
Q=eye(3);
R= 50;
[K,S,e] = lqr(sys.A,sys.B,Q,R)

%run simulation for different event trigger rule
c0_array=linspace(0,0.4,10); % 0<C0<1
[p,totalRun]=size(c0_array);
Tsimu=5;
dt=0.005; %very small sample time that system is linear for it

for iter = 1:totalRun
    C0 = c0_array(iter);
    fprintf('step %d of %d: simulation for C0= %d:\n', iter, totalRun ,C0);
    
    %initialize variable for each run
    x0=[1 2 3]'; %initial state 
    x=x0;
    t=0.0;
    y=0.0;
    yb=0.0;
    r=0.0; % set point    
    u=zeros(nb(2),1);
    uold=zeros(nb(2),1);
    dist=0; % disturbance   
    lastevent=0;  %save the time of the last event
    n=round(Tsimu/dt);
    t_array = zeros(1,n);
    u_array = zeros(nb(2),n);
    event_array =  zeros(1,n);  
    event_time_array = zeros(1,n); %save event time for calculation of sample time
    r_array =  zeros(nc(1),n);  
    y_array =  zeros(nc(1),n);
    x_array =  zeros(na(1),n);  
    x_error_array=  zeros(1,n);
    normX_array=  zeros(1,n);
    Xnew=zeros(na(1),1);
    Xnew=x0;
    Xold=zeros(na(1),1);

    % start of simulation loop 
    for i=1:n
         t = t + dt;       
         if mod(i,100)==0 
             fprintf('\nRunning... Time: %f',t) 
         end    
         
         x_error = Xnew-Xold;  
         x_error_array(i) = abs(norm(x_error));
         normX_array(i)= (norm(Xnew));
         if (C0 <= abs(norm(x_error)))
             % event occured so generate new u
             u=r-K*Xnew;
             
             %save data for plot curve
             u_array(:, i)=u;
             event_array(i)=1;
             event_time_array(i)= (i- lastevent)*dt ; %save event time for calculation of sample time
             lastevent=i;
             uold = u; %save u for next step
         else
             % event not triggered so use old u
             u = uold; 
             
             %save data for plot curve
             u_array(:, i)=u;
             event_array(i)=0;
         end
         
         %find system responce and new x
         xdot = sys.A*x + sys.B*u;
         x = x + xdot*dt;
         Xnew=x;
         y = sys.C*x;


         if event_array(i)==1
            Xold=Xnew;             
         end    
         
         %save result
         x_array(:,i)=Xnew;
         t_array(i)=t;
         y_array(i)=y;
         r_array(i)=r;
         
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
    plot(t_array,x_error_array,  t_array, normX_array  ),hold on;
    plot (t_array,C0*(ones(size(x_error_array))),'--'),hold off
    legend('x_error' ,'normX', 'error threshold')
    xlabel('time(s)');
    
    subplot(3,1,2)    
    stem(t_array,event_array)
    xlabel('time(s)');
    legend('event')
    grid on
    
    subplot(3,1,3)    
    stem(t_array,event_time_array)
    legend('sample time')
    xlabel('time(s)');
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
subplot(2,1,1)
set(gca,'FontSize',10);
plot(c0_array,R_array,'LineWidth',2,'Color','b');
grid on
xlabel('Sigma'); ylabel('R'); title('Event Ratio (R) vs Error Criterion (C0)');
subplot(2,1,2)
plot(c0_array,error_array,'LineWidth',2,'Color','r');
grid on
xlabel('Sigma'); ylabel('abs(error)'); title(' Output Error vs Error Criterion (C0)');
print('sigma_vs_R_and_E','-dpng');
