classdef RodPhotoReceptor < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        t
        dt
        buffer_size
        flag
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
		BH = 200;               % [uM] Total high-affinity buffer concentration
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
        function obj = RodPhotoReceptor(Y0, buffer_size, dt)
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
                       
            obj.Y = zeros(23, buffer_size);

			% variables
			obj.Y(:, 1) = Y0;
        end
		
		function update_time(obj)
            if obj.t+1 < obj.buffer_size
                obj.t = obj.t + 1;
            end
		end
        
        function [y, c]  = solve(obj,jhv)
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
			iKv = obj.gKv*obj.Y(2, k)^3*obj.Y(3, k)*(obj.Y(1,k)-obj.Ek);
			%%% Ca %%%
			Eca = -12.5*log(obj.Y(11, k)/obj.Ca0);
            %Eca = -12.5;
			[am_Ca, bm_Ca] = Ca(obj.Y(1, k));
			hCa = exp((40-obj.Y(1, k))/18)/(1+exp((40-obj.Y(1, k))/18));
			iCa = obj.gCa*obj.Y(4, k)^4*hCa*(obj.Y(1, k)-Eca);
			%%% Cl_Ca %%%
			mCl_Ca = 1/(1+exp((0.37-obj.Y(11, k))/0.09));
			iCl = obj.gClCa*mCl_Ca*(obj.Y(1, k)-obj.Eclca);
			%%% K_Ca %%%
			[am_KCa, bm_KCa] = K_Ca(obj.Y(1, k));
			mKCas = obj.Y(11, k)/(obj.Y(11, k)+0.3);
			iKCa = obj.gKCa*obj.Y(5, k)^2*mKCas*(obj.Y(1, k)-obj.Ek);
			%%% h %%%
			[ah, bh] = h(obj.Y(1, k));
			ih = obj.gh*(obj.Y(8, k)+obj.Y(9, k)+obj.Y(10, k))*(obj.Y(1, k)-obj.Eh);
			%%% L %%%
			iL =  obj.gL*(obj.Y(1, k)-obj.EL);
			%%% Cas %%%
			iex = obj.jex*exp(-(obj.Y(1, k)+14)/70)*(obj.Y(11, k)-obj.Cae)/(obj.Y(11, k)-obj.Cae+2.3);
			iex2 = obj.jex2*(obj.Y(11, k)-obj.Cae)/(obj.Y(11, k)-obj.Cae+0.5);
			%%% photo %%%
			j = obj.jmax*(obj.Y(23, k))^3/(obj.Y(23, k)^3+10^3);
			iPhoto = -j*(1.0-exp((obj.Y(1, k)-8.5)/17.0));
			
            
            D = zeros(23, 1);
            D(1) = 1/obj.Cm*-(iKv+iCa+iCl+iKCa+ih+iL+iPhoto+iex+iex2);
            D(2) = am_Kv*(1-obj.Y(2, k))-bm_Kv*obj.Y(2, k);
            D(3) = ah_Kv*(1-obj.Y(3, k))-bh_Kv*obj.Y(3, k);
            D(4) = am_Ca*(1-obj.Y(4, k))-bm_Ca*obj.Y(4, k);
            D(5) = am_KCa*(1-obj.Y(5, k))-bm_KCa*obj.Y(5, k);
            D(6) = -4*ah*obj.Y(6, k)+bh*obj.Y(7, k);
            D(7) = 4*ah*obj.Y(6, k)-3*(ah+bh)*obj.Y(7, k)+2*bh*obj.Y(8, k);
            D(8) = 3*ah*obj.Y(7, k)-(2*ah+2*bh)*obj.Y(8, k)+3*bh*obj.Y(9, k);
            D(9) = 2*ah*obj.Y(8, k)-(ah+3*bh)*obj.Y(9, k)+4*bh*obj.Y(10, k);
            D(10) = ah*obj.Y(9, k)-4*bh*obj.Y(10, k);
            D(11) = -(iCa+iex+iex2)/(2*obj.F*obj.V1)*10^(-6)...
                -obj.DCa*obj.S1/(obj.delta*obj.V1)*(obj.Y(11, k)-obj.Y(12, k))...
                -obj.Lb1*obj.Y(11, k)*(obj.BL-obj.Y(13, k))+obj.Lb2*obj.Y(13, k)...
                -obj.Hb1*obj.Y(11, k)*(obj.BH-obj.Y(14, k))+obj.Hb2*obj.Y(14, k);
            D(12) = obj.DCa*obj.S1/(obj.delta*obj.V2)*(obj.Y(11, k)-...
                obj.Y(12, k))-obj.Lb1*obj.Y(12, k)*(obj.BL-obj.Y(15, k))+...
                obj.Lb2*obj.Y(15, k)-obj.Hb1*obj.Y(12, k)*(obj.BH-obj.Y(16, k))...
                +obj.Hb2*obj.Y(16, k);
            D(13) = obj.Lb1*obj.Y(11, k)*(obj.BL-obj.Y(13, k))-obj.Lb2*obj.Y(13, k);
            D(14) = obj.Hb1*obj.Y(11, k)*(obj.BH-obj.Y(14, k))-obj.Hb2*obj.Y(14, k);
            D(15) = obj.Lb1*obj.Y(12, k)*(obj.BL-obj.Y(15, k))-obj.Lb2*obj.Y(15, k);
            D(16) = obj.Hb1*obj.Y(12, k)*(obj.BH-obj.Y(16, k))-obj.Hb2*obj.Y(16, k);
            D(17) = jhv-obj.a1*obj.Y(17, k)+obj.a2*obj.Y(18, k);
            D(18) = obj.a1*obj.Y(17, k)-(obj.a2+obj.a3)*obj.Y(18, k);
            D(19) = obj.e*obj.Y(17, k)*(obj.Ttot-obj.Y(19, k))-obj.b1*obj.Y(19, k)...
                -obj.tau1*obj.Y(19, k)*(obj.PDEtot-obj.Y(20, k))+obj.tau2*obj.Y(20, k);
            D(20) = obj.tau1*obj.Y(19, k)*(obj.PDEtot-obj.Y(20, k))-obj.tau2*obj.Y(20, k);
            D(21) = obj.b*j-obj.gammaCa*(obj.Y(21, k)-obj.Ca0p)...
                -obj.k1*(obj.eT-obj.Y(22, k))*obj.Y(21, k)+obj.k2*obj.Y(22, k);
            D(22) = obj.k1*(obj.eT-obj.Y(22, k))*obj.Y(21, k)-obj.k2*obj.Y(22, k);
            D(23) = obj.Amax/(1+(obj.Y(21, k)/obj.Kc)^4)-obj.Y(23, k)*...
                (obj.Vbar+obj.sigma*obj.Y(20, k));
           
            % values of variables at the next time step
            
            % we have already reached maximum buffer size so we should discard
            % primary values of the variables
            if (obj.t+1 == obj.buffer_size) && (obj.flag == 1)
                obj.Y = [obj.Y(:, 2:end), obj.Y(:, k)+obj.dt*D];
                y = obj.Y(:, end);
            % we have just reached maximum buffer size so e should change
            % flag value to 1
            elseif (obj.t+1 == obj.buffer_size) && (obj.flag == 0)
                obj.Y(:, k+1) = obj.Y(:, k)+obj.dt*D;
                obj.flag = 1;
                y = obj.Y(:, k+1);
            % we have not yet reached the buffer size
            else
                obj.Y(:, k+1) = obj.Y(:, k)+obj.dt*D;
                y = obj.Y(:, k+1);
            end
            
            c = [am_Kv, bm_Kv, ah_Kv, bh_Kv, iKv, Eca, am_Ca, bm_Ca, hCa, iCa,...
                mCl_Ca, iCl, am_KCa, bm_KCa, mKCas, iKCa, ah, bh, ih, iL, ...
                iex, iex2, j, iPhoto];
            
        end
        
        function v = get_V(obj)
            v = obj.Y(1, :);
        end
        
        function cgmp = get_cGMP(obj)
            cgmp = obj.Y(23, :);
        end
        
        function cas = get_Cas(obj)
            cas = obj.Y(11, :);
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

