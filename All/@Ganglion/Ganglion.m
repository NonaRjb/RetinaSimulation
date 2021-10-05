classdef Ganglion < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        t
        dt
        buffer_size
        flag
        t_vec
        method
        % constants ganglion
        tscale = 1000; 
		C = 1;      % [uF/cm^2]
		gNa = 53; %50;   % [mS/cm^2]
		ENa = 35;   % [mV]
		gA = 36;    % [mS/cm^2]
		gCa = 2.2;  % [mS/cm^2]
		gK = 12;    % [mS/cm^2]
        EK = -75;   % [mV]
		gKCa = 0.05;% [mS/cm^2]
        gL = 0.05;  % [mS/cm^2]
        El = -62.0; % [mV]
		R = 8.31;
        T = 293;
        F = 96500;
        Cae = 1800;
        Cadiss = 1.0;
        r = 12.5;
        Cares = 0.1;
        tCa = 1.5; %50;
		% chemical synapse constants
        gmax = 2.56;              % [nS] Maximum synapse conductance
        Esyn = 0.0;               % [mV] Synapse's reversal potential
        tau = 0.01;               % [ms] Time constant
        Vslope = 20;              % [mV] Voltage sensitivity of the synapse
        Vth = -36.424516776897130;% [mV]
        % variables
        Y
    end
    
    methods
        function obj = Ganglion(Y0, buffer_size, dt, method)
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
            obj.Y = zeros(9, buffer_size);

			% variables
			obj.Y(:, 1) = Y0;
        end
		
		function update_time(obj)
            if obj.t+1 < obj.buffer_size
                obj.t = obj.t + 1;
            end
		end
        
        function [y, curr_t, c]  = solve(obj,Vpre)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
			
            if (obj.t+1 == obj.buffer_size) && (obj.flag == 1)
                k = obj.t+1;
            else
                k = obj.t;
            end
			% explicit functions
			%%% Na %%%
			am_Na = -0.6*(obj.Y(1, k)+30.0)/(exp(-0.1*(obj.Y(1, k)+30.0))-1.0);
            bm_Na = 20.0*exp(-(obj.Y(1, k)+55.0)/18.0);
            ah_Na = 0.4*exp(-(obj.Y(1, k)+50.0)/20.0);
            bh_Na = 6.0/(exp(-0.1*(obj.Y(1, k)+20.0))+1.0);
            iNa = obj.gNa*obj.Y(2, k)^3*obj.Y(3, k)*(obj.Y(1, k)-obj.ENa);
			%%% A %%%
			am_A = (-0.006*(obj.Y(1, k)+90.0))/(exp(-0.1*(obj.Y(1, k)+90.0))-1.0);
			bm_A = 0.1*exp(-(obj.Y(1, k)+30.0)/10.0);
            ah_A = 0.04*exp(-(obj.Y(1, k)+70.0)/20.0);
            bh_A = 0.6/(exp(-0.1*(obj.Y(1, k)+40.0))+1.0);
            iA = obj.gA*obj.Y(7, k)^3*obj.Y(8, k)*(obj.Y(1, k)-obj.EK);
            %%% K %%%
			an_K = (-0.02*(obj.Y(1, k)+40.0))/(exp(-0.1*(obj.Y(1, k)+40.0))-1.0);
			bn_K = 0.4*exp(-(obj.Y(1, k)+50.0)/80.0);
            iK = obj.gK*obj.Y(6, k)^4*(obj.Y(1, k)-obj.EK);
            %%% Ca %%%
			Eca = 1000.0*(obj.R*obj.T/(2*obj.F))*log(obj.Cae/obj.Y(5, k));
			an_Ca = -0.3*(obj.Y(1, k)+13.0)/(exp(-0.1*(obj.Y(1, k)+13.0))-1.0);
			bn_Ca = 10.0*exp(-(obj.Y(1, k)+38.0)/18.0);
			iCa = obj.gCa*obj.Y(4, k)^3*(obj.Y(1, k)-Eca);
			%%% K_Ca %%%
            Caidiss2 = (obj.Y(5, k)/obj.Cadiss)^2;
			iKCa = obj.gKCa*Caidiss2/(1.0+Caidiss2)*(obj.Y(1, k)-obj.EK);
			%%% L %%%
			iL = obj.gL*(obj.Y(1, k)-obj.El);
            %%% Isyn %%%
%             if Vpre > obj.Vth
%                 S_inf = tanh((Vpre-obj.Vth)/obj.Vslope);
%             else
%                 S_inf = 0;
%             end
            S_inf = tanh(abs(Vpre-obj.Vth)/obj.Vslope);
            Isyn = obj.gmax*obj.Y(9, k)*(obj.Y(1, k)-obj.Esyn);
			
            consts = [obj.C; obj.F; obj.r; obj.Cares; obj.tCa; obj.tau;...
                obj.Vth; obj.tscale];
            
            vars = [iNa; iA; iK; iCa; iKCa; iL; am_Na; bm_Na; ah_Na; bh_Na;...
                an_Ca; bn_Ca; an_K; bn_K; am_A; bm_A; ah_A; bh_A; Isyn;...
                S_inf; Vpre];
            
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
            
            
            c = [am_Na, bm_Na, ah_Na, bh_Na, iNa, am_A, bm_A, ah_A, ...
                bh_A, iA, an_K, bn_K, iK, Eca, an_Ca, bn_Ca, iCa,...
                Caidiss2, iKCa, iL, S_inf, Isyn];
            
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
        
        function set_Vth(obj, vth)
            obj.Vth = vth;
        end
    end
end


%%%% F %%%%
function D = f(Y, vars, consts)

% consts
C = consts(1);
F = consts(2);
r = consts(3);
Cares = consts(4);
tCa = consts(5);
tau = consts(6);
Vth = consts(7);
tscale = consts(8);

% vars
iNa = vars(1);
iA = vars(2);
iK = vars(3);
iCa = vars(4);
iKCa = vars(5);
iL = vars(6);
am_Na = vars(7);
bm_Na = vars(8);
ah_Na = vars(9);
bh_Na = vars(10);
an_Ca = vars(11);
bn_Ca = vars(12);
an_K = vars(13);
bn_K = vars(14);
am_A = vars(15);
bm_A = vars(16);
ah_A = vars(17);
bh_A = vars(18);
Isyn = vars(19);
S_inf = vars(20);
Vpre = vars(21);

D = zeros(9, 1);
D(1) = ((-(iNa+iA+iK+iCa+iKCa+iL)*20-Isyn)/(20*C))*tscale;
D(2) = (am_Na*(1.0-Y(2))-bm_Na*Y(2))*tscale;
D(3) = (ah_Na*(1.0-Y(3))-bh_Na*Y(3))*tscale;
D(4) = (an_Ca*(1.0-Y(4))-bn_Ca*Y(4))*tscale;
D(5) = (-(((10^4)*3.0*iCa)/(2.0*F*r))-((Y(5)-Cares)/tCa))*tscale;
D(6) = (an_K*(1.0-Y(6))-bn_K*Y(6))*tscale;
D(7) = (am_A*(1.0-Y(7))-bm_A*Y(7))*tscale;
D(8) = (ah_A*(1.0-Y(8))-bh_A*Y(8))*tscale;
if Vpre <= Vth
    D(9) = 0;
elseif abs(Vpre-Vth)>= 0
    D(9) = (S_inf-Y(9))/((1-S_inf)*tau);
else
    D(9) = 0;
end




end
