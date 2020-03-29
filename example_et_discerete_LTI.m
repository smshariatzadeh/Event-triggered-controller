% example_et_discrete_LTI.m
%
% event trigger for discrete-time linear time-invariant system
%
% simulation is run for different event trigger parameters i.e.: in sigma(event-triggered parameter)  
%
% event triggered rule is based on comparing norm X with norm e
%
% sigma : event-triggered parameter in event rule 
%
% Ref: Event-triggered control for Discrete-Time Systems
% Alina Eqtami, et al (ACC2010)
%         
% By S.M.Shariatzadeh
% Date :20 Mars 2020 


%% define linear time-invarient system
sys.A= [0.1   1.2  ;
        0.007 1.05 ];
sys.B= [300  200  ;
        0.5  0.001];
sys.C= [1 1];

%optimal vontroller gain
Q=0.001*eye(2);
R= [ 0.01 0.01;
     0.0  0.01];

%% calculate system data & Find Dimensions of A, B , Q, R
%
[na,ma]=size(sys.A);
[nb,mb]=size(sys.B);
[nq,mq] =size (Q) ;
[nr,mr]=size(R);
%[nf ,mf] =size (F) ;
if na~=ma %Verify A is square
    error('A must besquare')
else
    [na, na] =size (sys.A) ;
end
%
%Data Checks for proper system setup
% if length(sys.A)>rank(ctrb(sys.A,sys.B))
%     %Check for controllability
%     error('System Not Controllable')
%     return
% end
if (na ~= nq) || (na ~= mq)
   %Check that A and Q are the same size
   error('A and Q must be the same size');
   return
end
% if ~(mf==1&nf==1)
%     if (nq ~= nf) || (mq ~= mf)
%         %Check that Q and F are the same size
%         error('Q and F must be the same size');
%         return
%     end
% end
if ~(mr==1&nr==1)
    if (mr ~= nr) || (mb ~= nr)
        error('R must be consistent with B');
        return
    end
end
mq = norm(Q,1);
% Check if Q is positive semi-definite and symmetric
if any(eig(Q) < -eps*mq) || (norm(Q'-Q,1)/mq > eps)
    disp('Warning: Q is not symmetric and positive ... semi-definite');
end
mr = norm(R,1);
% Check if R is positive definite and symmetric
if any(eig(R) <= -eps*mr) || (norm(R'-R,1)/mr > eps)
    disp('Warning: R is not symmetric and positive ... definite');
end


% designing feedback gain based on Linear-Quadratic Regulator (LQR)  
[K,S,e] = dlqr(sys.A,sys.B,Q,R)

%run simulation for different event trigger rule
sigma_array=linspace(0,0.8,20); % 0<sigma<1
[p,totalRun]=size(sigma_array);
Tsimu=10;
dt=0.5; % sample time 

for iter = 1:totalRun
    Sigma=sigma_array(iter);
    fprintf('step %d of %d: simulation for sigma %d:\n', iter, totalRun ,Sigma);
    
    %initialize variable for each run
    x0=[-0.2  0.5]'; %initial state 
    x=x0;
    t=0.0;
    y=0.0;
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
    r_array =  zeros(nc(1),n);  
    y_array =  zeros(nc(1),n);
    x_array =  zeros(na,n);  
    x_error_array=  zeros(1,n);
    normX_array=  zeros(1,n);
    snormX_array=  zeros(1,n);
    Xnew=zeros(na,1);
    Xold=zeros(na,1);
    Xnew=x0;
   
    % start of simulation loop 
    for i=1:n
         t = t + dt;       
         if mod(i,100)==0 
             fprintf('\nRunning... Time: %f',t) 
         end    
         
         x_error = Xnew-Xold;  
         x_error_array(i) = abs(norm(x_error));
         snormX_array(i)= Sigma*(norm(Xnew));
         normX_array(i)= (norm(Xnew));
         if (Sigma*(norm(Xnew))<= abs(norm(x_error)))
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
         xplus1 = sys.A*x + sys.B*u;
         x = xplus1;
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
    stem(t_array,event_array,'r','fill');
    title('Event-triggered')
    xlabel('time(s)');
    
    figure(2)
    clf
    subplot(3,1,1)
    stem(t_array,x_error_array,'fill'),hold on
    plot(t_array, snormX_array, '-', 'Marker' ,'.')
    plot(t_array, normX_array , '--', 'Marker' ,'.' ),hold off
    legend('x_error','sigma*normX','normX')
    xlabel('time(s)');    
    
    subplot(3,1,2)    
    stem(t_array,event_array,'r','fill')
    legend('event')
    xlabel('time(s)');
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
plot(sigma_array,R_array,'LineWidth',2,'Color','b');
grid on
xlabel('Sigma'); ylabel('R'); title('Event Ratio (R) vs Error Criterion (Sigma)');
subplot(2,1,2)
plot(sigma_array,error_array,'LineWidth',2,'Color','r');
grid on
xlabel('Sigma'); ylabel('abs(error)'); title(' Output Error vs Error Criterion (Sigma)');
print('sigma_vs_R_and_E','-dpng');
