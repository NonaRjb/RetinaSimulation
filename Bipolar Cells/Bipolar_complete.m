classdef Bipolar_complete < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        t
        dt
        buffer_size
        flag
        t_vec
        method
        % constants bipolar
		C = 0.01;   % [pF]
		gKv = 2.0;  % [nS]
		Ek = -58;    % [mV]
		gA = 50;    % [nS]
		gh = 0.975; % [nS]
		Eh = -17.7; % [mV]
		gCa = 1.1;  % [nS]
		Ca0 = 2500; % [uM]
		gKCa = 8.5; % [nS]
		gL = 0.23;  % [nS]
		El = -21.0;   % [mV]
		F = 9.649e04;   % [cmol^(-1)] Faraday constant
		DCa = 6e-08;    % [dm^2*sec^(-1)] Ca diffusion coefficient 
		Vs = 1.692e-13; % [dm^(-3)] Volume of submembrane area 
		Vd = 7.356e-13; % [dm^(-3)] Volume of the deep intracellular area
		Ssd = 3.94e-08;    % [dm^(-2)] Surface area of the submembrane and the deep intracellular area spherical boundary
		dsd = 5.8e-05;   % [dm] Distance between submembrane area and the deep intracellular area
		Cablmax_s = 300;  % [uM] Total low-affinity buffer concentration
		Cabhmax_s = 100;  % [uM] Total high-affinity buffer concentration
		Cablmax_d = 500;  % [uM] Total low-affinity buffer concentration
		Cabhmax_d = 300;  % [uM] Total high-affinity buffer concentration
		abl = 0.4;      % [sec^(-1)uM^(-1)] On rate constants for the binding of Ca to low-affinity buffer
		bbl = 0.2;      % [sec^(-1)uM^(-1)] Off rate constants for the binding of Ca to low-affinity buffer
		abh = 100;      % [sec^(-1)uM^(-1)] On rate constants for the binding of Ca to high-affinity buffer
		bbh = 90;       % [sec^(-1)uM^(-1)] Off rate constants for the binding of Ca to high-affinity buffer
		Jex = 9.0*40;        % [pA] Maximum Na-Ca exchanger current
		Jex2 = 9.5*40;     % [pA] Maximum Ca-ATPase exchanger current
		Camin = 0.01;   % [uM] Minimum intracellular Ca concentration for Ca extrusion
		% chemical synapse constants
        gmax = 2.56;    % [nS] Maximum synapse conductance
        Esyn = 0.0;       % [mV] Synapse's reversal potential
        tau = 0.01;       % [ms] Time constant
        Vslope = 20;    % [mV] Voltage sensitivity of the synapse
        Vth = -36.185963;      % [mV]
        % variables
        Y
    end
    
    methods
        function obj = Bipolar_complete(Y0, buffer_size, dt, method)
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
            obj.Y = zeros(19, buffer_size);

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
			%%% Kv %%%
			[am_Kv, bm_Kv, ah_Kv, bh_Kv] = Kv(obj.Y(1, k));
			iKv = obj.gKv*obj.Y(2, k)^3*obj.Y(3, k)*(obj.Y(1, k)-obj.Ek);
			%%% A %%%
			[am_A, bm_A, ah_A, bh_A] = A(obj.Y(1, k));
			iA = obj.gA*obj.Y(4, k)^3*obj.Y(5, k)*(obj.Y(1, k)-obj.Ek);
			%%% h %%%
			[ah, bh] = h(obj.Y(1, k));
			ih = obj.gh*(obj.Y(8, k)+obj.Y(9, k)+obj.Y(10, k))*(obj.Y(1, k)-obj.Eh);
			%%% Ca %%%
			Eca = 12.9*log(obj.Ca0/obj.Y(13, k));
			[am_Ca, bm_Ca] = Ca(obj.Y(1, k));
			hCa = exp(-(obj.Y(1, k)-50)/11)/(exp(-(obj.Y(1, k)-50)/11)+1.0);
			iCa = obj.gCa*obj.Y(11, k)^4*hCa*(obj.Y(1, k)-Eca);
			%%% K_Ca %%%
			[am_KCa, bm_KCa] = K_Ca(obj.Y(1, k));
			mKc1 = obj.Y(13, k)/(obj.Y(13, k)+0.2);
			iKCa = obj.gKCa*obj.Y(12, k)^2*mKc1*(obj.Y(1, k)-obj.Ek);
			%%% L %%%
			iL = obj.gL*(obj.Y(1, k)-obj.El);
			%%% Iex %%%
			Iex = obj.Jex*(obj.Y(13, k)-obj.Camin)/(obj.Y(13, k)...
                -obj.Camin+2.3)*exp(-(obj.Y(1, k)+14.0)/70.0);
			Iex2 = obj.Jex2*(obj.Y(13, k)-obj.Camin)/(obj.Y(13, k)-obj.Camin+0.5);
            %%% Isyn %%%
%             if Vpre > obj.Vth
%                 S_inf = tanh((Vpre-obj.Vth)/obj.Vslope);
%             else
%                 S_inf = 0;
%             end
            S_inf = tanh(abs(Vpre-obj.Vth)/obj.Vslope);
            %Isyn = obj.gmax*obj.Y(19, k)*(obj.Y(1, k)-obj.Esyn);
            Isyn = Vpre;
			
            consts = [obj.C; obj.F; obj.Vs; obj.Vd; obj.DCa; obj.Ssd;... 
                obj.dsd; obj.bbl; obj.bbh; obj.abl; obj.abh; obj.Cablmax_s;...
                obj.Cabhmax_s; obj.Cabhmax_s; obj.Cabhmax_s; obj.tau; obj.Vth];
            
            vars = [iKv; iA; ih; iCa; iKCa; iL; am_Kv; bm_Kv; ah_Kv; bh_Kv;...
                am_A; bm_A; ah_A; bh_A; ah; bh; am_Ca; bm_Ca; am_KCa;...
                bm_KCa; Iex; Iex2; Isyn; S_inf; Vpre];
            
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
            
            
            c = [am_Kv, bm_Kv, ah_Kv, bh_Kv, iKv, am_A, bm_A, ah_A, ...
                bh_A, iA, ah, bh, ih, Eca, am_Ca, bm_Ca, hCa, iCa,...
                am_KCa, bm_KCa, mKc1, iKCa, iL, Iex, Iex2, S_inf, Isyn];
            
        end
        
        function v = get_V(obj)
            v = obj.Y(1, :);
        end
        
        function s = get_S(obj)
            s = obj.Y(19, :);
        end
        
        function cas = get_Cas(obj)
            cas = obj.Y(13, :);
        end
        
        function tVector = get_tvec(obj)
            tVector = obj.t_vec;
        end
        
        function set_Vth(obj, vth)
            obj.Vth = vth;
        end
    end
end


%%%%% Kv %%%%%
function [am_Kv, bm_Kv, ah_Kv, bh_Kv] = Kv(V)

am_Kv = 400/(exp(-(V-15)/36)+1.0);
bm_Kv = 1.0 * exp(-V/13);
ah_Kv = 0.0003*exp(-V/7);
bh_Kv = 80/(exp((V+115)/15)+1.0)+0.02;

end
%%%%%%%%%%%%%%%%

%%%%% A %%%%%
function [am_A, bm_A, ah_A, bh_A] = A(V)

am_A = 2400/(exp(-(V-50.0)/28.0)+1.0);
bm_A = 12*exp(-V/10.0);
ah_A = 0.045*exp(-V/13);
bh_A = 75/(exp(-(V+30)/15)+1.0);

end
%%%%%%%%%%%%%%%%

%%%%% h %%%%%
function [ah, bh] = h(V)

ah = 3/(exp((V+110)/15)+1);
bh = 1.5/(exp(-(V+115)/15)+1);

end
%%%%%%%%%%%%%%%%


%%%%% Ca %%%%%
function [am_Ca, bm_Ca] = Ca(V)

am_Ca = 12000*(120-V)/(exp((120-V)/25)-1.0);
bm_Ca = 40000/(exp((V+68)/25)+1.0);

end
%%%%%%%%%%%%%%%%


%%%%% K_Ca %%%%%
function [am_KCa, bm_KCa] = K_Ca(V)

am_KCa = 100*(230-V)/(exp((230-V)/52)-1.0);
bm_KCa = 120*exp(-V/95);

end
%%%%%%%%%%%%%%%%

%%%% F %%%%
function D = f(Y, vars, consts)

% consts
C = consts(1);
F = consts(2);
Vs = consts(3);
Vd = consts(4);
DCa = consts(5);
Ssd = consts(6);
dsd = consts(7);
bbl = consts(8);
bbh = consts(9);
abl = consts(10);
abh = consts(11);
Cablmax_s = consts(12);
Cabhmax_s = consts(13);
Cablmax_d = consts(14);
Cabhmax_d = consts(15);
tau = consts(16);
Vth = consts(17);

% vars
iKv = vars(1);
iA = vars(2);
ih = vars(3);
iCa = vars(4);
iKCa = vars(5);
iL = vars(6);
am_Kv = vars(7);
bm_Kv = vars(8);
ah_Kv = vars(9);
bh_Kv = vars(10);
am_A = vars(11);
bm_A = vars(12);
ah_A = vars(13);
bh_A = vars(14);
ah = vars(15);
bh = vars(16);
am_Ca = vars(17);
bm_Ca = vars(18);
am_KCa = vars(19);
bm_KCa = vars(20);
Iex = vars(21);
Iex2 = vars(22);
Isyn = vars(23);
S_inf = vars(24);
Vpre = vars(25);


D = zeros(19, 1);
D(1) = 1/C*-(iKv+ih+iCa+iKCa+iL-Isyn);
D(2) = am_Kv*(1.0-Y(2))-bm_Kv*Y(2);
D(3) = ah_Kv*(1.0-Y(3))-bh_Kv*Y(3);
D(4) = am_A*(1.0-Y(4))-bm_A*Y(4);
D(5) = ah_A*(1.0-Y(5))-bh_A*Y(5);
D(6) = -4.0*ah*Y(6)+bh*Y(7);
D(7) = 4.0*ah*Y(6)-(3.0*ah+bh)*Y(7)+2.0*bh*Y(8);
D(8) = 3.0*ah*Y(7)-2.0*(ah+bh)*Y(8)+3.0*bh*Y(9);
D(9) = 2.0*ah*Y(8)-(ah+3.0*bh)*Y(9)+4.0*bh*Y(10);
D(10) = ah*Y(9)-4.0*bh*Y(10);
D(11) = am_Ca*(1.0-Y(11))-bm_Ca*Y(11);
D(12) = am_KCa*(1.0-Y(12))-bm_KCa*Y(12);
D(13) = -(iCa+Iex+Iex2)/(2*F*Vs)*10^(-6)-DCa*Ssd/(Vs*dsd)*(Y(13)-Y(14))+...
    bbl*Y(15)-abl*Y(13)*(Cablmax_s-Y(15))+bbh*Y(16)-abh*Y(13)*(Cabhmax_s-Y(16));
D(14) = DCa*Ssd/(Vd*dsd)*(Y(13)-Y(14))+bbl*Y(17)-abl*Y(14)*(Cablmax_d-Y(17))+...
    bbh*Y(18)-abh*Y(14)*(Cabhmax_d-Y(18));
D(15) = abl*Y(13)*(Cablmax_s-Y(15))-bbl*Y(15);
D(16) = abh*Y(13)*(Cabhmax_s-Y(16))-bbh*Y(16);
D(17) = abl*Y(14)*(Cablmax_d-Y(17))-bbl*Y(17);
D(18) = abh*Y(14)*(Cabhmax_d-Y(18))-bbh*Y(18);
if Vpre >= Vth
    D(19) = 0;
elseif abs(Vpre-Vth)>= 0
    D(19) = (S_inf-Y(19))/((1-S_inf)*tau);
else
    D(19) = 0;
end



end
