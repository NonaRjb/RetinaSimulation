classdef Amacrine < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        t
        dt
        buffer_size
        flag
        t_vec
        % constants ganglion
        tscale = 1000; 
		C = 20;     % [uF/cm^2]
		gNa = 4;    % [mS/cm^2]
		ENa = 40;   % [mV]
		gK = 0.4;   % [mS/cm^2]
        EK = -80;   % [mV]
        gL = 0.46;  % [mS/cm^2]
        El = -54.0; % [mV]
        celcius_Na_A = 22;
        celcius_K_A = 22;
		% chemical synapse constants
        gmax = 1.2;               % [nS] Maximum synapse conductance
        Esyn = -10.0;             % [mV] Synapse's reversal potential
        tau = 0.01;               % [ms] Time constant
        Vslope = 20;              % [mV] Voltage sensitivity of the synapse
        Vth = -36.424516776897130;% [mV]
        % variables
        Y
    end
    
    methods
        function obj = Amacrine(Y0, buffer_size, dt)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            if nargin==1
                buffer_size = 10000;
                dt = 1e-06;
            elseif nargin == 2
                dt = 1e-06;
            end
            
            obj.t = 1;
            obj.dt = dt;
            obj.buffer_size = buffer_size;
            obj.flag = 0;
            
            obj.t_vec = zeros(buffer_size, 1); 
            obj.Y = zeros(5, buffer_size);

			% variables
			obj.Y(:, 1) = Y0;
        end
		
		function update_time(obj)
            if obj.t+1 < obj.buffer_size
                obj.t = obj.t + 1;
            end
		end
        
        function [y, curr_t, c]  = solve(obj, Vpre)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
			
            if (obj.t+1 == obj.buffer_size) && (obj.flag == 1)
                k = obj.t+1;
            else
                k = obj.t;
            end
			% explicit functions
			%%% Na %%%
            [am_A, bm_A] = m(obj.Y(1, k), obj.celcius_Na_A);
            [ah_A, bh_A] = h(obj.Y(1, k), obj.celcius_Na_A);
			tau_m_Na = 1/(am_A+bm_A);
            m_inf_Na = am_A/(am_A+bm_A);
            tau_h_Na = 1/(ah_A/3.5+bh_A/3.5);
            h_inf_Na = ah_A/(ah_A+bh_A);
            iNa = obj.gNa*obj.Y(2, k)^3*obj.Y(3, k)*(obj.Y(1, k)-obj.ENa);
            %%% K %%%
            [an_A, bn_A] = n(V, obj.celcius_K_A);
            tau_n_K = 1/(an_A/5+bn_A/5);
            n_inf_K = an_A/(an_A+bn_A);
            iK = obj.gK*obj.Y(4, k)^4*(obj.Y(1, k)-obj.EK);
			%%% L %%%
			iL = obj.gL*(obj.Y(1, k)-obj.El);
            %%% Isyn %%%
%             if Vpre > obj.Vth
%                 S_inf = tanh((Vpre-obj.Vth)/obj.Vslope);
%             else
%                 S_inf = 0;
%             end
            S_inf = tanh(abs(Vpre-obj.Vth)/obj.Vslope);
            Isyn = obj.gmax*obj.Y(5, k)*(obj.Y(1, k)-obj.Esyn);
			
            consts = [obj.C; obj.tau; obj.Vth; obj.tscale];
            
            vars = [iNa; iK; iL; tau_m_Na; m_inf_Na; tau_h_Na; ...
                h_inf_Na; tau_n_K; n_inf_K; Isyn; S_inf; Vpre];
            
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
            
            %%%%% runge-kutta 4
%             k1 = obj.dt * D;
%             k2 = obj.dt * f(obj.Y(:, k)+k1/2, vars, consts);
%             k3 = obj.dt * f(obj.Y(:, k)+k2/2, vars, consts);
%             k4 = obj.dt * f(obj.Y(:, k)+k3, vars, consts);
%             k_tot = 1/6*(k1+2*k2+2*k3+k4);
            %%%%% forward euler
            k_tot = obj.dt*D;
            
            
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
            
            
            c = [am_A, bm_A, ah_A, bh_A, iNa, an_A, bn_A, iK, ...
                iL, S_inf, Isyn];
            
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

%%%% Q10_A %%%%
function fun = q10_A(temp)

fun = 3.0 ^ ((temp-6.3)/10.0);

end

%%%% Expm1_A %%%%
function fun = expm1_A(x, y)

if abs(x/y) < 1e-06
    fun = y*(1-x/(2*y));
else
    fun = x/(exp(x/y)-1);
end

end

%%%% m_A %%%%
function [am_A, bm_A] = m(V, temp_Na)

V = -V-65;
am_A = q10_A(temp_Na)*0.1*expm1_A((V+25),10);
bm_A = q10_A(temp_Na)*4*exp(V/18);

end

%%%% h_A %%%%
function [ah_A, bh_A] = h(V, temp_Na)

V = -V-65;
ah_A = q10_A(temp_Na)*0.07*exp(V/20);
bh_A = q10_A(temp_Na)*1/(exp(0.1*V+3)+1);

end

%%%% n_A %%%%
function [an_A, bn_A] = n(V, temp_K)

V = -V-65;
an_A = q10_A(temp_K)*0.01*expm1_A(V+10,10);
bn_A = q10_A(temp_K)*0.125*exp(V/80);

end

%%%% F %%%%
function D = f(Y, vars, consts)

% consts
C = consts(1);
tau = consts(2);
Vth = consts(3);
tscale = consts(4);

% vars
iNa = vars(1);
iK = vars(2);
iL = vars(3);
tau_m_Na = vars(4);
m_inf_Na = vars(5);
tau_h_Na = vars(6);
h_inf_Na = vars(7);
tau_n_K = vars(8);
n_inf_K = vars(9);
Isyn = vars(10);
S_inf = vars(11);
Vpre = vars(12);

D = zeros(5, 1);
D(1) = ((-(iNa+iA+iK+iCa+iKCa+iL)*20-Isyn)/(20*C))*tscale;
D(2) = (m_inf_Na-Y(2))/tau_m_Na*tscale;
D(3) = (h_inf_Na-Y(3))/tau_h_Na*tscale;
D(4) = (n_inf_K-Y(4))/tau_n_K*tscale;
if Vpre < Vth
    D(5) = 0;
elseif abs(Vpre-Vth)>= 0
    D(5) = (S_inf-Y(9))/((1-S_inf)*tau);
else
    D(5) = 0;
end




end
