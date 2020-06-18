classdef Graphics<handle
% Library class for definition of different output and plot of simulation
%
% This is a function collection rather than a class
% So, All methods are static.
% plot trajectory and show convex set which are projected from n-dim
%
% printSolution =arguments:
%       iprint:     = 0  Print closed loop data(default)
%                   = 1  Print closed loop data and errors of the
%                        optimization method
%                   = 2  Print closed loop data and errors and warnings of
%                        the optimization method
%                   >= 5 Print closed loop data and errors and warnings of
%                        the optimization method as well as graphical
%                        output of closed loop state trajectories
%                   >=10 Print closed loop data and errors and warnings of
%                        the optimization method with error and warning
%                        description
%
%
%   printHeader:         Clarifying header for selective output of closed
%                        loop data, cf. printClosedloopData
%
%   printClosedloopData: Selective output of closed loop data
%
%   plotTrajectories:    Graphical output of the trajectories, requires
%                        iprint >= 4
%
    
    
    methods (Static)
        function show_convex(P, varargin)
            P_reduced = projectPolytope2Plane(P);
            switch numel(varargin)
                case 1
                    fill(P_reduced.V(:, 1), P_reduced.V(:, 2), varargin{1});
                case 3
                    fill(P_reduced.V(:, 1), P_reduced.V(:, 2), varargin{1}, varargin{2}, varargin{3});
            end
            hold on;
        end

        function show_convex_timeslice(P, z_axis, varargin)
            P_reduced = projectPolytope2Plane(P);
            num_vertex = size(P_reduced.V, 1);
            vertex_3d = [P_reduced.V, z_axis*ones(num_vertex, 1)];
            switch numel(varargin)
                case 1
                    fill3(vertex_3d(:, 1), vertex_3d(:, 2), vertex_3d(:, 3), varargin{1});
                case 3
                    fill3(vertex_3d(:, 1), vertex_3d(:, 2), vertex_3d(:, 3), varargin{1}, varargin{2}, varargin{3});
            end
            hold on;
        end
        
        function show_convex_sliding(S1, S2, h1, h2, varargin)
            % s1 and s2 must has a similarity.
           S1_reduced = projectPolytope2Plane(S1);
           S2_reduced = projectPolytope2Plane(S2);
           V1 = S1_reduced.V;
           V2 = S2_reduced.V;
           num_vert = size(V1, 1);
           for i = 1:num_vert
               plot3([V1(i, 1); V2(i, 1)], [V1(i, 2); V2(i, 2)], [h1, h2], varargin{1})  
           end
           hold on;
        end

        function show_trajectory(x_seq, varargin)
            switch numel(varargin)
                case 1
                    plot(x_seq(1, :), x_seq(2, :), varargin{1});
                case 3
                    plot(x_seq(1, :), x_seq(2, :), varargin{1}, varargin{2}, varargin{3});     
            end
            hold on;
        end
        
        function show_trajectory_timeslice( x_seq, varargin)
            num_xseq = size(x_seq, 2);
            x_seq_3d = [x_seq(1, :); x_seq(2, :); flip(0:(num_xseq-1))];
            switch numel(varargin)
                case 1
                    plot3(x_seq_3d(1, :), x_seq_3d(2, :), x_seq_3d(3, :), varargin{1});
                case 3
                    plot3(x_seq_3d(1, :), x_seq_3d(2, :), x_seq_3d(3, :), varargin{1}, varargin{2}, varargin{3});     
            end
            hold on;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Definition of output text format
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function printHeader()
            fprintf('   k  |      u(k)        x(1)        x(2)       x(3)        x(4)     Time\n');
            fprintf('-------------------------------------------------------------------------------------------------\n');
        end

        function printClosedloopData(mpciter, u, x, t_Elapsed , message)
            if size(x,1)==2 
            fprintf(' %3d  | %+11.6f %+11.6f %+11.6f                      %+6.3f %10s',  ...
                     mpciter, u(1,1), x(1), x(2), t_Elapsed, message);
            elseif size(x,1)==3
            fprintf(' %3d  | %+11.6f %+11.6f %+11.6f %+11.6f            %+6.3f %10s',  ...
                     mpciter, u(1,1), x(1), x(2), x(3),  t_Elapsed, message);
            elseif size(x,1)==4 
            fprintf(' %3d  | %+11.6f %+11.6f %+11.6f %+11.6f %+11.6f  %+6.3f %10s',  ...
                     mpciter, u(1,1), x(1), x(2), x(3), x(4), t_Elapsed, message);
                 
            end
        end
        
        function plotTrajectories(u, t_intermediate, x_intermediate)
            

            figure(2);
            subplot(221)
            plot(t_intermediate,x_intermediate(1,:),'-ok');
            title(['x_1 closed loop trajectory']);
            xlabel('n');
            ylabel('x_1(n)');
            grid on;
            hold on;
            subplot(222)
            plot(t_intermediate,x_intermediate(2,:),'-ok');
            title(['x_2 closed loop trajectory']);
            xlabel('n');
            ylabel(' x_2(n)');
            grid on;
            hold on;
            if size(x_intermediate,1)>=3   
                subplot(223)
                plot(t_intermediate,x_intermediate(3,:),'-ok');
                title([' x_3 closed loop trajectory']);
                xlabel('n');
                ylabel(' x_3(n)');
                grid on;
                hold on;
            end
            if size(x_intermediate, 1 )>=4               
                subplot(224)
                plot(t_intermediate,x_intermediate(4,:),'-ok');
                title([' x_4 closed loop trajectory']);
                xlabel('n');
                ylabel(' x_4(n)');
                grid on;
                hold on;        
            end


        end
        
        function plotTrajectoriesXY( x_intermediate, A , B)
            
            figure(1)
            plot(x_intermediate(A,:),x_intermediate(B,:),'or', ...
                 'MarkerFaceColor','r');
            axis square;
            title('x_2/x_4 closed loop trajectory');
            xlabel('x_2');
            ylabel('x_4');
            grid on;
            hold on;
            
        end
        
        function printSolution(mpciter, iprint, exitflag, t_Elapsed, u,t,x, message)
            if (mpciter == 1)
                Graphics.printHeader();
            end
            Graphics.printClosedloopData(mpciter, u(:,mpciter), x(:,mpciter), t_Elapsed, message);
            switch exitflag
                case -2
                if ( iprint >= 1 && iprint < 10 )
                    fprintf(' Error F\n');
                elseif ( iprint >= 10 )
                    fprintf(' Error: No feasible point was found\n')
                end
                case -1
                if ( iprint >= 1 && iprint < 10 )
                    fprintf(' Error OT\n');
                elseif ( iprint >= 10 )
                    fprintf([' Error: The output function terminated the',...
                             ' algorithm\n'])
                end
                case 0
                if ( iprint == 1 )
                    fprintf('\n');
                elseif ( iprint >= 2 && iprint < 10 )
                    fprintf(' Warning IT\n');
                elseif ( iprint >= 10 )
                    fprintf([' Warning: Number of iterations exceeded',...
                             ' options.MaxIter or number of function',...
                             ' evaluations exceeded options.FunEvals\n'])
                end
                case 1
                if ( iprint == 1 )
                    fprintf('\n');
                elseif ( iprint >= 2 && iprint < 10 )
                    fprintf(' \n');
                elseif ( iprint >= 10 )
                    fprintf([' First-order optimality measure was less',...
                             ' than options.TolFun, and maximum constraint',...
                             ' violation was less than options.TolCon\n'])
                end
                case 2
                if ( iprint == 1 )
                    fprintf('\n');
                elseif ( iprint >= 2 && iprint < 10 )
                    fprintf(' Warning TX\n');
                elseif ( iprint >= 10 )
                    fprintf(' Warning: Change in x was less than options.TolX\n')
                end
                case 3
                if ( iprint == 1 )
                    fprintf('\n');
                elseif ( iprint >= 2 && iprint < 10 )
                    fprintf(' Warning TJ\n');
                elseif ( iprint >= 10 )
                    fprintf([' Warning: Change in the objective function',...
                             ' value was less than options.TolFun\n'])
                end
                case 4
                if ( iprint == 1 )
                    fprintf('\n');
                elseif ( iprint >= 2 && iprint < 10 )
                    fprintf(' Warning S\n');
                elseif ( iprint >= 10 )
                    fprintf([' Warning: Magnitude of the search direction',...
                             ' was less than 2*options.TolX and constraint',...
                             ' violation was less than options.TolCon\n'])
                end
                case 5
                if ( iprint == 1 )
                    fprintf('\n');
                elseif ( iprint >= 2 && iprint < 10 )
                    fprintf(' Warning D\n');
                elseif ( iprint >= 10 )
                    fprintf([' Warning: Magnitude of directional derivative',...
                             ' in search direction was less than',...
                             ' 2*options.TolFun and maximum constraint',...
                             ' violation was less than options.TolCon\n'])
                end
            end
            if ( iprint >= 5 )
                Graphics.plotTrajectories( u, t, x)
            end
        end
        
        
        
    end
end

function P_projected = projectPolytope2Plane(P)
    vert = P.V;
    x_vert = round(vert(:, 1), 5);
    y_vert = round(vert(:, 2), 5);
    idx = convhull(x_vert, y_vert);
    P_projected = Polyhedron([x_vert(idx), y_vert(idx)]);
end
