classdef EventGenerator < handle
    % EventGenerator , a class for Event Generator
    % 
    %  smshariatzadeh@yahoo.com  
    %
    properties (SetAccess = private)
        % common notation for linear system
        
        nx;     % dim of state space 
        Xold;   % event generator memory that hold state of system durring the last event
        Sigma;  % event generator parameters
    end 

    methods (Access = public)

        function obj = EventGenerator(Sigma, nx )
            % EventGenerator Construct an instance of this class
            % choose method of event generation
            obj.Sigma = Sigma;
            obj.nx = nx;
            obj.Xold = zeros(size(nx, 1 )); % for generate event in the first time step
        end

        function [Event, x_error, snormX, normX]  = Check_event(obj, NewX)
             % Read new state and calculate some variable and check for event
             sx = size(NewX);
             if sx(1) ~= obj.nx
                error ('New state is not consistante with the system state space') 
             end    
             
             % calculation
             x_error = NewX-obj.Xold;  
             x_error = abs(norm(x_error));
             snormX = obj.Sigma*(norm(NewX));
             normX = (norm(NewX));
             
             % check event
             if (obj.Sigma*normX)<= abs(norm(x_error)) && normX >0.01

                 %event generated
                 %
                 Event = 1;
                 obj.Xold = NewX;  %save X in event generator memory for next step
             else
                 % event not triggered 
                 Event = 0;

             end
        end
    end
end

