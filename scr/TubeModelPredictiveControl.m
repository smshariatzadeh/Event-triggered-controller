classdef TubeModelPredictiveControl
% a class to simulation of Tube Model Predictive Control for discsrete-time linear
% system

% use Matlab R2020a & MPT3

    properties (SetAccess = public)
        sys % MMPT3 LTIsystem
        optcon; % optimal contol solver object
        Xc;
        Uc;
        Xc_robust; % Xc-Z (Pontryagin diff.)
        Uc_robust; % Uc-K*Z (Pontryagin diff.)
        W; % disturbance
        Z; % disturbance invariant set (minimal robust positvely invariant set =mRPI)
        Xmpi;
        Xmpi_robust; % Xmpi-Z (Pontryagin diff. : Maximal robust positvely invariant set=MPRI)
        N;
        nc;
        
    end
    
    methods (Access = public)
        function obj = TubeModelPredictiveControl(sys, Q, R, Xc, Uc, W,  N,  varargin)
            %fill properteis
            obj.sys = sys;
            obj.W = W;
            obj.Xc = Xc;
            obj.Uc = Uc;
            obj.W = W;

            [K_tmp, P] = dlqr(sys.A, sys.B, Q, R);
            K = -K_tmp;
            Ak = (sys.A + sys.B * K);

            %----------approximation of d-inv set--------%
            Z = obj.ApproxmRPIset( Ak, 3, 1.01);
            
            %create robust X and U constraints, and construct solver using X and U
            Xc_robust = Xc - Z;
            Uc_robust = Uc - K*Z;

            optcon = OptimalControler(sys, Q, R, Xc_robust, Uc_robust, N);
            obj.optcon = optcon;
            
            optcon.reconstruct_ineq_constraint(Xc_robust, Uc_robust);

            %robustize Xmpi set and set it as a terminal constraint
            Xmpi=obj.compute_MPIset();
            Xmpi_robust = Xmpi - Z;
            optcon.add_terminal_constraint(Xmpi_robust);

            obj.Xc_robust = Xc_robust;
            obj.Z = Z;
            obj.Xmpi=Xmpi;
            obj.Xmpi_robust = Xmpi_robust;
            obj.N = N;
           
        end
                           
        function Z_approx = ApproxmRPIset(obj, Ak, n_order, alpha)
            %compute approximate minimal robust positively invarient set (Z)
            % which takes the form of Z = alpha*(W + Ak*W + Ak^2*W + ... Ak^n_ordr*W).
            % where + denotes Minkowski addition.
            %
            %input:
            %   Ak: closed loop system matrix 
            %   W : system noise (as Polyhedron ) (from obj)
            %   n : degree of approximate (n >=1)
            %       choose bigger n for best approximation 
            %   alpha :a magnifier constant  (alpha>=1)
            %
            %output:
            %  Zapprox: approximation of minimal robust positively invarient set as
            %           polyhedron
            %
            % We could obtain approximate minimal RPI (Z) by computing an infinite geometric series,
            % which is not practicall to get.
            % we approximate this by trancating the polynomial.
            
            if (n_order <1)
               error ('degree of approximate must be (n >=1)')
            end
            if (alpha<1)
               error ('magnifier constant must be (alpha>=1)') 
            end
            Z_approx = obj.W;
            
            for n = 1:n_order
                Z_approx = Z_approx + Ak^n*obj.W;
            end
            Z_approx = Z_approx*alpha;
            % which takes the form of Z = alpha*(W + Ak*W + Ak^2*W + ... Ak^n_ordr*W).
            % where + denotes Minkowski addition.
        end
        
        function Xmpi = compute_MPIset(obj)
            %Compute Maximal robust positvely invariant set 
            % MPIset is computed only once in the constructor;
            %
            [F, G, obj.nc] = convert_Poly2Mat(obj.Xc, obj.Uc);
            K=obj.optcon.K;
            Ak=obj.optcon.Ak;
            Fpi = @(i) (F+G*(K))*(Ak)^i;
            Xpi = @(i) Polyhedron(Fpi(i), ones(size(Fpi(i), 1), 1));
            obj.Xmpi = Xpi(0);
            i= 0;
            while(1) % 
                i = i + 1;
                Xmpi_tmp = and(obj.Xmpi, Xpi(i));
                if Xmpi_tmp == obj.Xmpi
                    break;
                else
                    Xmpi = Xmpi_tmp;
                end
            end
            
        end
        
        function [] = simulate(obj, Tsimu, x_init)
            obj.optcon.remove_initial_eq_constraint()
            Xinit = x_init+obj.Z;
            obj.optcon.add_initial_constraint(Xinit);

            [x_nominal_seq, u_nominal_seq] = obj.optcon.solve(x_init);
            
            x = x_init;
            x_real_seq = [x];
            u_real_seq = [];
            propagate = @(x, u, w) obj.sys.A*x+obj.sys.B*u + w;
            % simulation loop
            for i=1:Tsimu
                if i<=obj.N
                    u = u_nominal_seq(:, i) + obj.sys.K*(x-x_nominal_seq(:, i));
                else 
                    u = obj.sys.K*x;
                end
                w = rand(2, 1).*(obj.w_max - obj.w_min)+obj.w_min;
                x = propagate(x, u, w);
                x_real_seq = [x_real_seq, x];
                u_real_seq = [u_real_seq, u];

                clf; % real time plot
                Graphics.show_convex(obj.Xc, 'm');
                Graphics.show_convex(obj.Xc_robust, 'r');
                Graphics.show_convex(obj.optcon.Xmpi, [0.2, 0.2, 0.2]*1.5);
                Graphics.show_convex(obj.Xmpi_robust, [0.5, 0.5, 0.5]); % gray
                for j=1:obj.N+1
                    Graphics.show_convex(x_nominal_seq(:, j)+obj.Z, 'g', 'FaceAlpha', 0.3);
                end
                Graphics.show_trajectory(x_nominal_seq, 'gs-');
                Graphics.show_trajectory(x, 'b*-');
                pause(0.2)
            end

            % time slice plot after simulation
            figure(2)
            Graphics.show_convex_timeslice(obj.Xc, -0.04, 'm');
            Graphics.show_convex_timeslice(obj.Xc_robust, -0.03, 'r');
            Graphics.show_convex_timeslice(obj.optcon.Xmpi, -0.02, [0.2, 0.2, 0.2]*1.5);
            Graphics.show_convex_timeslice(obj.Xmpi_robust, -0.01, [0.5, 0.5, 0.5]);
            Graphics.show_convex_timeslice(x_nominal_seq(:, 1)+obj.Z, obj.N, 'g', 'FaceAlpha', .3);
            Graphics.show_trajectory_timeslice(x_nominal_seq, 'gs-', 'LineWidth', 1.2);
            Graphics.show_trajectory_timeslice(x_real_seq(:, 1:obj.N+1), 'b*-', 'LineWidth', 1.2);
            leg = legend('$X_c$', '$X_c\ominus Z$', '$X_f (= X_{MPI})$', '$X_f\ominus Z$', 'Tube', 'Nominal', 'Real');
            set(leg, 'Interpreter', 'latex')
            
            for i=2:obj.N+1 % show remaining tubes.
                    Graphics.show_convex_timeslice(x_nominal_seq(:, i)+obj.Z, obj.N-i+1, 'g', 'FaceAlpha', .3);
            end
            
            xlabel('x1');
            ylabel('x2');
            zlabel('time (minus)');
            xlim([obj.optcon.x_min(1), obj.optcon.x_max(1)]);
            ylim([obj.optcon.x_min(2), obj.optcon.x_max(2)]);
            zlim([-0.05, obj.N]);
            set( gca, 'FontSize',12); 
            view([10, 10])
            grid on;
        end
    end
end
