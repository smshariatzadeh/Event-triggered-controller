% example_mRPI.m

% example of minimal Robust Positively Invariant(mRPI) Set for closed loop
% optimal feedback
%
%
% 1. Finding minimal Robust Positively Invariant Set for closed loop
% optimal control by different approximation
%
% 2. Simulation of a closed loop control system  
%
% 3. Show graphically that state is inside the mRPI
% 
% requires  mpt3 toolbox
% 
% Copyright 2019-2024 smshariatzadeh@yahoo.com .
%
close all

%% make your own discrete linear system (plant)
A = [.9 0.5; 0 0.7];
B = [0.5; 1]; 
C = [0.5  0.2]; 
D= [0];
f=[0;0];
g=[0];
    
myopnlopsys = LTISystem('A', A, 'B', B, 'C', C, 'D', D, 'f', f, 'g', g);% discrete-time state-space model

% constraints on state Xc and input Uc
Xc_vertex = [2, -2; 2 2; -10 2; -10 -2];
Uc_vertex = [1; -1];
Xc = Polyhedron(Xc_vertex);
Uc = Polyhedron(Uc_vertex);

% construct a convex set of system noise (2dim here)
W_vertex = [0.25, 0.25; 0.25, -0.25; -0.25, -0.25; -0.25, 0.25]; %low nois
%W_vertex = [0.15, 0.15; 0.15, -0.15; -0.15, -0.15; -0.15, 0.15];
W = Polyhedron(W_vertex);
% set boundary of system noise which corresponds to W.
w_max = [0.25; 0.25];
w_min = [-0.25; -0.25];

figure(1)
plot(W)
title("noise of system"); 

%% make optimal controller
Q = diag([1, 1]);
R = 0.1;
K = -dlqr(myopnlopsys.A, myopnlopsys.B, Q, R);

%% closed loop system

Ak= A+B*K;
closedsys = LTISystem('A', Ak, 'B', B, 'C', C, 'D', D, 'f', f, 'g', g);% discrete-time state-space model

Q = diag([1, 1]); 
R = 0.1;
[K_tmp, P] = dlqr(myopnlopsys.A, myopnlopsys.B, Q, R);
K = -K_tmp;
Ak = (myopnlopsys.A + myopnlopsys.B * K);

% compute minimal disturbance(robust) invariant set Z.
Z = approx_minRPIset(Ak, W, 3, 1.45);
figure(20)
title(" minimal disturbance invariant set for controller"); 
plot(Z)



%% make observer
%The command lqr can be adapted to calculate an optimal observer gain in a dual way: L = lqr(A',C',Q,R)'
Q = diag([1, 1]); 
R = 0.1;
L = dlqr(myopnlopsys.A', myopnlopsys.C', Q, R)';
ob = LTIObserver(myopnlopsys,L);

% compute observer invariant set (Zob)
Zob=ob.ApproxmRPIset(W,2,1.001);


%% find minimal Robust Positively Invariant (mRPI) Set for closed loop system 
% based on controller RPI and observer RPI


Zsys = Z + Zob;
figure(3)
hold on
Graphics.show_convex(Zsys, 'g', 'FaceAlpha', .4); % show Z
Graphics.show_convex(Z, 'r', 'FaceAlpha', .8); % show Z
Graphics.show_convex(Zob, 'b', 'FaceAlpha', .3); % show Z


hold off
legend(" Z system", "Z controller" ,"Z observer"  )

%approximation of Z controller

myclosedloopsys.A=Ak;
myclosedloopsys.f=[0;0];
myclosedloopsys.nx=size(Ak,1);
alpha=.2;
[Z_inner,Z_outer] = minrpiset(myclosedloopsys,W,alpha);
figure(2)
hold on
Graphics.show_convex(Z_outer,'m' , 'FaceAlpha', .4); % show Z
Graphics.show_convex(Z_inner,'r' , 'FaceAlpha', .7); % show Z
hold off
legend("Z outer approx.", "Z inner approx." )
title("approximation of minimal disturbance invariant set for controller"); 
%% simulation

% propagate particles many times following the discrete stochastic dynamics,
% and we will see that particles never go outside of distervance invariant set Z.
Nptcl = 100;
x = zeros(2, Nptcl); % particles 
y = myopnlopsys.C*x;
u=0;
t=0;
while true
     %% observer estimation
     Xhat=  ob.EstimateCurrentState( y , u);
     Yhat = myopnlopsys.C*Xhat;
    
    %% optimal control state feed back based on observer estimation
    u = K * Xhat;
    
    %% system response
    sys_noise = rand(2, Nptcl).*repmat(w_max - w_min, 1, Nptcl) + repmat(w_min, 1, Nptcl);
    x = myopnlopsys.A*x + myopnlopsys.B*u +sys_noise;
    y = myopnlopsys.C*x;
    %% show result
    titl = sprintf('Running... Time: %f (press any key to exit)',t );
    
    % Note that system state is placed inside the Z of controller.
    h=figure(4);
    clf;
    Graphics.show_convex(Zsys,'g' , 'FaceAlpha', .3); % show Z
    Graphics.show_convex(Z, 'r', 'FaceAlpha', .7); % show Z    
    scatter(x(1, :), x(2, :)); % show particles
    legend("Z system(including observer)", "Z controller ")
    title(titl); 
    drawnow
    isKeyPressed = ~isempty(get(h,'CurrentCharacter'));
    if isKeyPressed
       break
    end    
    pause(0.01);
    t=t+1;
    
end
