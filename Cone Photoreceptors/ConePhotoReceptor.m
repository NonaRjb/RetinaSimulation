classdef ConePhotoReceptor < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        t
        dt
        buffer_size
        flag
        t_vec
        method
        % constants cone
		E = -72.0;      % [mV]
        tL = 10;        % [ms]
        gf_bar = 0.2;
        gi_bar = 0.7;
        tf = 67;        % [ms]
        Vf = -57.0;     % [mV]
        Veps = 4.0;     % [mV]
        alpha = 83.3;   % [sec^(-1)]
        k12_bar = 10;   % [sec^(-1)]
        k12_M = 250;    % [sec^(-1)]
        k23 = 17;       % [sec^(-1)]
        k34 = 1.3;      % [sec^(-1)]
        k32 = 0.03;     % [sec^(-1)]
        A = 6.76;       % k12/k21
        nu = 2.125;     % [sec^(-2)]
        K = 10;         % [sec]
        offset;
        % variables
        Y
    end
    
    methods
        function obj = ConePhotoReceptor(Y0, buffer_size, dt, method)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            if nargin==1
                buffer_size = 10000;
                dt = 1e-06;
                method = 'euler';
            elseif nargin == 2
                dt = 1e-06;
                method = 'euler';
            elseif nargin == 3
                method = 'euler';
            end
            
            obj.t = 1;
            obj.dt = dt;
            obj.buffer_size = buffer_size;
            obj.flag = 0;
            obj.method = method;
            
            obj.t_vec = zeros(buffer_size, 1); 
            obj.Y = zeros(10, buffer_size);

			% variables
			obj.Y(:, 1) = Y0;
            obj.offset = -Y0(1);
        end
		
		function update_time(obj)
            if obj.t+1 < obj.buffer_size
                obj.t = obj.t + 1;
            end
		end
        
        function [y, curr_t, c]  = solve(obj,iC)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
			
            if (obj.t+1 == obj.buffer_size) && (obj.flag == 1)
                k = obj.t+1;
            else
                k = obj.t;
            end
			% explicit functions
			
            Gi = obj.gi_bar/(1+obj.Y(8, k)/obj.K);
            F = obj.gf_bar/(1+exp((obj.Y(1, k)-obj.Vf)/obj.Veps));
            k12 = obj.k12_bar+obj.nu*obj.Y(9, k)*(obj.k12_M-obj.k12_bar)...
                /(obj.k12_M-obj.k12_bar+obj.nu*obj.Y(9, k));
            k21 = k12/obj.A;
			
            consts = [obj.E; obj.tL; obj.tf; obj.alpha; obj.k23; obj.k32;...
                obj.k34; obj.offset];
            
            vars = [Gi; F; iC; k12; k21];
            
            D = f(obj.Y(:,k), vars, consts);
            
            %%%% specify dt %%%%
%             
%             if abs(max(D, [], 'all')) > 300
%                 obj.dt = 1e-06;
%             elseif abs(max(D, [], 'all')) <= 300 && abs(max(D, [], 'all')) > 200
%                 obj.dt = 2*1e-06;
%             elseif abs(max(D, [], 'all')) <= 200 && abs(max(D, [], 'all')) > 100
%                 obj.dt = 5*1e-06;
%             elseif abs(max(D, [], 'all')) <= 100 && abs(max(D, [], 'all')) > 10
%                 obj.dt = 1e-05;
%             elseif abs(max(D, [], 'all')) <= 100 && abs(max(D, [], 'all')) > 10
%                 obj.dt = 5*1e-05;
%             else
%                 obj.dt = 1e-04;
%             end
            
%             maxD = abs(max(D, [], 'all'));
%             if maxD > 1e-03
%                 obj.dt = 1/maxD*1e-03;
%             else
%                 obj.dt = 0.1;
%             end
            %objdt = obj.dt;
            
            if strcmp(obj.method, 'rk4')
            %%%%% runge-kutta 4
                k1 = obj.dt * D;
                k2 = obj.dt * f(obj.Y(:, k)+k1/2, vars, consts);
                k3 = obj.dt * f(obj.Y(:, k)+k2/2, vars, consts);
                k4 = obj.dt * f(obj.Y(:, k)+k3, vars, consts);
                k_tot = 1/6*(k1+2*k2+2*k3+k4);
            elseif strcmp(obj.method, 'euler')
            %%%%% forward euler
                k_tot = obj.dt*D;
            else
                k_tot = obj.dt*D;
            end
            
            
            % values of variables at the next time step
            
            % we have already reached maximum buffer size so we should discard
            % primary values of the variables
            if (obj.t+1 == obj.buffer_size) && (obj.flag == 1)
                obj.Y = [obj.Y(:, 2:end), obj.Y(:, k)+k_tot];
                obj.t_vec = [obj.t_vec(2:end); obj.t_vec(end)+obj.dt];
                y = obj.Y(:, end);
                curr_t = obj.t_vec(end);
            % we have just reached maximum buffer size so e should change
            % flag value to 1
            elseif (obj.t+1 == obj.buffer_size) && (obj.flag == 0)
                obj.Y(:, k+1) = obj.Y(:, k)+k_tot;
                obj.t_vec(k+1) = obj.t_vec(k)+obj.dt;
                obj.flag = 1;
                y = obj.Y(:, k+1);
                curr_t = obj.t_vec(k+1);
            % we have not yet reached the buffer size
            else
                obj.Y(:, k+1) = obj.Y(:, k)+k_tot;
                obj.t_vec(k+1) = obj.t_vec(k)+obj.dt;
                y = obj.Y(:, k+1);
                curr_t = obj.t_vec(k+1);
            end
            
            
            c = [Gi; F; k12; k21];
            
        end
        
        function v = get_V(obj)
            v = obj.Y(1, :);
        end
        
        function s = get_S(obj)
            s = obj.Y(9, :);
        end
        
        function tVector = get_tvec(obj)
            tVector = obj.t_vec;
        end
    end
end


%%%% F %%%%
function D = f(Y, vars, consts)

% consts
E = consts(1);
tL = consts(2);
tf = consts(3);
alpha = consts(4);
k23 = consts(5);
k32 = consts(6);
k34 = consts(7);
offset = consts(8);

% vars
Gi = vars(1);
F = vars(2);
iC = vars(3);
k12 = vars(4);
k21 = vars(5);

D = zeros(10, 1);
D(1) = 1000*(E-Y(1)*(1+Y(2)+Gi))/tL;
D(2) = 1000*(F-Y(2))/tf;
D(3) = iC-alpha*Y(3);
D(4) = alpha*(Y(3)-Y(4));
D(5) = alpha*(Y(4)-Y(5));
D(6) = alpha*(Y(5)-Y(6));
D(7) = alpha*(Y(6)-Y(7));
D(8) = alpha*Y(7)-k12*Y(8)+k21*Y(9);
D(9) = k12*Y(8)-(k21+k23)*Y(9)+k32*Y(10);
D(10) = k23*Y(9)-(k32+k34)*Y(10);

end
