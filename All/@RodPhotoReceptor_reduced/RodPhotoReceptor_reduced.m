classdef RodPhotoReceptor_reduced < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        t
        dt
        buffer_size
        flag
        t_vec
        method
        % constants
		Cm = 0.02;      % [nF]
		gKv = 2.0;      % [nS]
		Ek = -74;       % [mV]
		gCa = 0.7;      % [nS]
		Ca0 = 1600;     % [uM]
		gClCa = 2.0;    % [nS]
		Eclca = -20;    % [mV]
		gKCa = 5.0;     % [nS]
		gh = 3.0;       % [nS]
		Eh = -32;       % [mV]
		gL = 0.35;      % [nS]
		EL = -77;       % [mV]
		F = 9.648*10^4;         % [cmol^(-1)] Faraday const.
		V1 = 3.812*10^(-13);    % [dm^3] Volume of submembrane area
		V2 = 5.236*10^(-13);    % [dm^3] Volume of deep intracellular area
		DCa = 6*10^(-8);        % [d(m^2)(s^(-1))] Ca diffusion coefficient 
		S1 = 3.142*10^(-8);     % [dm^2] Surface area of the submembrane and the deep intracellular area spherical boundary
		delta = 5.9*10^(-5);    % [dm] Distance between submembrane area and the deep intracellular area
		Lb1 = 0.4;              % [s^(-1)uM^(-1)] On rate constant for the binding of Ca to low-affinity buffer
		Lb2 = 0.2;              % [s^(-1)uM^(-1)] Off rate constant for the binding of Ca to low-affinity buffer
		Hb1 = 100;              % [s^(-1)uM^(-1)]  On rate constant for the binding of Ca to high-affinity buffer
		Hb2 = 90;               % [s^(-1)uM^(-1)]  off rate constants for the binding of Ca to high-affinity buffer
		BL = 500;               % [uM] Total low-affinity buffer concentration
		BH = 300;               % [uM] Total high-affinity buffer concentration
		jex = 20;               % [pA] Maximum Na-Ca exchanger current
		jex2 = 20;              % [pA] Maximum Ca-ATPase exchanger current
		Cae = 0.01;             % [uM] Minimum intracellular Ca2+ concentration
		a1 = 50.0;                % [s^(-1)] Rate constant of Rh* inactivation
		a2 = 0.0003;            % [s^(-1)] Rate constant of the reaction Rhi -> Rh*
		a3 = 0.03;              % [s^(-1)] Rate constant of the decay of inactive rhodopsin
		e = 0.5;                % [s^(-1)uM^(-1)] Rate constant of T* activation
		Ttot = 1000;            % [uM] Total transducin
		b1 = 2.5;               % [s^(-1)] Rate constant of T* inactivation
		tau1 = 0.2;             % [s^(-1)uM^(-1)] Rate constant of PDE activation
		tau2 = 5;               % [s^(-1)] Rate constant of PDE inactivation
		PDEtot = 100;           % [uM] Phosphodiasterase
		gammaCa = 50;           % [s^(-1)] Rate constant of Ca2+ extrusion in the absence of Ca2+ buffers mediated by the Na+?Ca2+ exchanger
		Ca0p = 0.1;             % [uM] Intracellular Ca2+ concentration at the steady state
		b = 0.25;               % [uMs^(-1)pA^(-1)] Proportionality constant between Ca2+ influx and photocurrent
		k1 = 0.2;               % [uM^(-1)] on rate constants for the binding of Ca2+ to the buffer
		k2 = 0.8;               % [uM^(-1)] off rate constants for the binding of Ca2+ to the buffer
		eT = 500;               % [uM] Low-affinity Ca2+ buffer concentration
		Vbar = 0.4;             % [s^(-1)] Cyclic GMP hydrolysis in dark
		Kc = 0.1;               % [uM] Maximal activity of guanylate cyclase
		Amax = 65.6;            % [uMs^(-1)] Maximal activity of guanylate cyclase
		sigma = 1.0;            % [s^(-1)uM^(-1)] Proportionality constant
		jmax = 5040;            % [pA] Maximal cyclic GMP gated current inexcised patches
		% variables
        Y
    end
    
    methods
        function obj = RodPhotoReceptor_reduced(Y0, buffer_size, dt, method)
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
            obj.Y = zeros(21, buffer_size);

			% variables
			obj.Y(:, 1) = Y0;
        end
		
		function update_time(obj)
            if obj.t+1 < obj.buffer_size
                obj.t = obj.t + 1;
            end
		end
        
        function [y, curr_t, c]  = solve(obj,jhv)
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
            %m_Kv = am_Kv/(am_Kv+bm_Kv);
			iKv = obj.gKv*(am_Kv/(am_Kv+bm_Kv))^3*obj.Y(2, k)*(obj.Y(1,k)-obj.Ek);
			%%% Ca %%%
			Eca = -12.5*log(obj.Y(9, k)/obj.Ca0);
			[am_Ca, bm_Ca] = Ca(obj.Y(1, k));
			hCa = exp((40-obj.Y(1, k))/18)/(1+exp((40-obj.Y(1, k))/18));
			iCa = obj.gCa*obj.Y(3, k)^4*hCa*(obj.Y(1, k)-Eca);
			%%% Cl_Ca %%%
			mCl_Ca = 1/(1+exp((0.37-obj.Y(9, k))/0.09));
			iCl = obj.gClCa*mCl_Ca*(obj.Y(1, k)-obj.Eclca);
			%%% K_Ca %%%
			[am_KCa, bm_KCa] = K_Ca(obj.Y(1, k));
			mKCas = obj.Y(9, k)/(obj.Y(9, k)+0.3);
            % m_KCa = am_KCa/(am_KCa+bm_KCa);
			iKCa = obj.gKCa*(am_KCa/(am_KCa+bm_KCa))^2*mKCas*(obj.Y(1, k)-obj.Ek);
			%%% h %%%
			[ah, bh] = h(obj.Y(1, k));
			ih = obj.gh*(obj.Y(6, k)+obj.Y(7, k)+obj.Y(8, k))*(obj.Y(1, k)-obj.Eh);
			%%% L %%%
			iL =  obj.gL*(obj.Y(1, k)-obj.EL);
			%%% Cas %%%
			iex = obj.jex*exp(-(obj.Y(1, k)+14)/70)*(obj.Y(9, k)-obj.Cae)/(obj.Y(9, k)-obj.Cae+2.3);
			iex2 = obj.jex2*(obj.Y(9, k)-obj.Cae)/(obj.Y(9, k)-obj.Cae+0.5);
			%%% photo %%%
			j = obj.jmax*(obj.Y(21, k))^3/(obj.Y(21, k)^3+10^3);
			iPhoto = -j*(1.0-exp((obj.Y(1, k)-8.5)/17.0));
			
            consts = [obj.Cm; obj.F; obj.V1; obj.V2; obj.DCa; obj.S1;... 
                obj.delta; obj.Lb1; obj.Lb2; obj.Hb1; obj.Hb2; obj.BL;...
                obj.BH; obj.a1; obj.a2; obj.a3; obj.e; obj.Ttot; ...
                obj.b1; obj.tau1; obj.tau2; obj.PDEtot; obj.b; obj.gammaCa;...
                obj.Ca0p; obj.k1; obj.k2; obj.eT; obj.Amax; obj.Kc;...
                obj.Vbar; obj.sigma];
            
            vars = [iKv; iCa; iCl; iKCa; ih; iL; iPhoto; iex; iex2;...
                ah_Kv; bh_Kv; am_Ca; bm_Ca; ah; bh; jhv; j];
            
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
            
            
            c = [am_Kv, bm_Kv, ah_Kv, bh_Kv, iKv, Eca, am_Ca, bm_Ca, hCa, iCa,...
                mCl_Ca, iCl, am_KCa, bm_KCa, mKCas, iKCa, ah, bh, ih, iL, ...
                iex, iex2, j, iPhoto];
            
        end
        
        function v = get_V(obj)
            v = obj.Y(1, :);
        end
        
        function cgmp = get_cGMP(obj)
            cgmp = obj.Y(21, :);
        end
        
        function cas = get_Cas(obj)
            cas = obj.Y(9, :);
        end
        
        function tVector = get_tvec(obj)
            tVector = obj.t_vec;
        end
    end
end


%%%%% Kv %%%%%
function [am_Kv, bm_Kv, ah_Kv, bh_Kv] = Kv(V)

am_Kv = 5*(100-V)/(exp((100-V)/42)-1);
bm_Kv = 9*exp(-(V-20)/40);
ah_Kv = 0.15*exp(-V/22);
bh_Kv = 0.4125/(exp((10-V)/7)+1);

end
%%%%%%%%%%%%%%%%

%%%%% Ca %%%%%
function [am_Ca, bm_Ca] = Ca(V)

am_Ca = 3*(80-V)/(exp((80-V)/25)-1);
bm_Ca = 10/(1+exp((V+38)/7));

end
%%%%%%%%%%%%%%%%

%%%%% K_Ca %%%%%
function [am_KCa, bm_KCa] = K_Ca(V)

am_KCa = 15*(80-V)/(exp((80-V)/40)-1);
bm_KCa = 20*exp(-V/35);

end
%%%%%%%%%%%%%%%%

%%%%% h %%%%%
function [ah, bh] = h(V)

ah = 8/(exp((V+78)/14)+1);
bh = 18/exp(-(V+8)/19+1);

end
%%%%%%%%%%%%%%%%

%%%% F %%%%
function D = f(Y, vars, consts)

% consts
Cm = consts(1);
F = consts(2);
V1 = consts(3);
V2 = consts(4);
DCa = consts(5);
S1 = consts(6);
delta = consts(7);
Lb1 = consts(8);
Lb2 = consts(9);
Hb1 = consts(10);
Hb2 = consts(11);
BL = consts(12);
BH = consts(13);
a1 = consts(14);
a2 = consts(15);
a3 = consts(16);
e = consts(17);
Ttot = consts(18);
b1 = consts(19);
tau1 = consts(20);
tau2 = consts(21);
PDEtot = consts(22);
b = consts(23);
gammaCa = consts(24);
Ca0p = consts(25);
k1 = consts(26);
k2 = consts(27);
eT = consts(28);
Amax = consts(29);
Kc = consts(30);
Vbar = consts(31);
sigma = consts(32);

% vars
iKv = vars(1);
iCa = vars(2);
iCl = vars(3);
iKCa = vars(4);
ih = vars(5);
iL = vars(6);
iPhoto = vars(7);
iex = vars(8);
iex2 = vars(9);
%am_Kv = vars(10);
%bm_Kv = vars(11);
ah_Kv = vars(10);
bh_Kv = vars(11);
am_Ca = vars(12);
bm_Ca = vars(13);
%am_KCa = vars(16);
%bm_KCa = vars(17);
ah = vars(14);
bh = vars(15);
jhv = vars(16);
j = vars(17);


D = zeros(21, 1);
D(1) = 1/Cm*-(iKv+iCa+iCl+iKCa+ih+iL+iPhoto+iex+iex2);
%D(2) = am_Kv*(1-Y(2))-bm_Kv*Y(2);
D(2) = ah_Kv*(1-Y(2))-bh_Kv*Y(2);
D(3) = am_Ca*(1-Y(3))-bm_Ca*Y(3);
%D(5) = am_KCa*(1-Y(5))-bm_KCa*Y(5);
D(4) = -4*ah*Y(4)+bh*Y(5);
D(5) = 4*ah*Y(4)-3*(ah+bh)*Y(5)+2*bh*Y(6);
D(6) = 3*ah*Y(5)-(2*ah+2*bh)*Y(6)+3*bh*Y(7);
D(7) = 2*ah*Y(6)-(ah+3*bh)*Y(7)+4*bh*Y(8);
D(8) = ah*Y(7)-4*bh*Y(8);
D(9) = -(iCa+iex+iex2)/(2*F*V1)*10^(-6)-DCa*S1/(delta*V1)*(Y(9)-Y(10))...
    -Lb1*Y(9)*(BL-Y(11))+Lb2*Y(11)-Hb1*Y(9)*(BH-Y(12))+Hb2*Y(12);
D(10) = DCa*S1/(delta*V2)*(Y(9)-Y(10))-Lb1*Y(10)*(BL-Y(13))+...
    Lb2*Y(13)-Hb1*Y(10)*(BH-Y(14))+Hb2*Y(14);
D(11) = Lb1*Y(9)*(BL-Y(11))-Lb2*Y(11);
D(12) = Hb1*Y(9)*(BH-Y(12))-Hb2*Y(12);
D(13) = Lb1*Y(10)*(BL-Y(13))-Lb2*Y(13);
D(14) = Hb1*Y(10)*(BH-Y(14))-Hb2*Y(14);
D(15) = jhv-a1*Y(15)+a2*Y(16);
D(16) = a1*Y(15)-(a2+a3)*Y(16);
D(17) = e*Y(15)*(Ttot-Y(17))-b1*Y(17)-tau1*Y(17)*(PDEtot-Y(18))+tau2*Y(18);
D(18) = tau1*Y(17)*(PDEtot-Y(18))-tau2*Y(18);
D(19) = b*j-gammaCa*(Y(19)-Ca0p)-k1*(eT-Y(20))*Y(19)+k2*Y(20);
D(20) = k1*(eT-Y(20))*Y(19)-k2*Y(20);
D(21) = Amax/(1+(Y(19)/Kc)^4)-Y(21)*(Vbar+sigma*Y(18));


end
