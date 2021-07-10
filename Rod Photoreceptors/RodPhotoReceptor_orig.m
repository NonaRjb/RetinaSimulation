classdef RodPhotoReceptor < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        t
        dt
        buffer_size
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
        V
        mKv
        hKv
        mCa
        mKCa
        C1
        C2
        O1
        O2
        O3
        Cas
        Caf
        Cabls
        Cabhs
        Cablf
        Cabhf
        Rh
        Rhi
        Tr
        PDE
        Ca2
        Cab
        cGMP
    end
    
    methods
        function obj = RodPhotoReceptor(Y, buffer_size, dt)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            if nargin==1
                buffer_size = 10000;
                dt = 1e-04;
            elseif nargin == 2
                dt = 1e-05;
            end
            
            obj.t = 1;
            obj.dt = dt;
            obj.buffer_size = buffer_size;
            
            obj.V = zeros(buffer_size, 1);
			obj.mKv = zeros(buffer_size, 1); 
			obj.hKv = zeros(buffer_size, 1);
			obj.mCa = zeros(buffer_size, 1);
			obj.mKCa = zeros(buffer_size, 1);
			obj.C1 = zeros(buffer_size, 1);
			obj.C2 = zeros(buffer_size, 1);
			obj.O1 = zeros(buffer_size, 1);
			obj.O2 = zeros(buffer_size, 1);
			obj.O3 = zeros(buffer_size, 1);
			obj.Cas = zeros(buffer_size, 1);
			obj.Caf = zeros(buffer_size, 1);
			obj.Cabls = zeros(buffer_size, 1);
			obj.Cabhs = zeros(buffer_size, 1);
			obj.Cablf = zeros(buffer_size, 1);
			obj.Cabhf = zeros(buffer_size, 1);
			obj.Rh = zeros(buffer_size, 1);
			obj.Rhi = zeros(buffer_size, 1);
			obj.Tr = zeros(buffer_size, 1);
			obj.PDE = zeros(buffer_size, 1);
			obj.Ca2 = zeros(buffer_size, 1);
			obj.Cab = zeros(buffer_size, 1);
			obj.cGMP = zeros(buffer_size, 1);

			% variables
			obj.V(1) = Y(1);
			obj.mKv(1) = Y(2); 
			obj.hKv(1) = Y(3);
			obj.mCa(1) = Y(4);
			obj.mKCa(1) = Y(5);
			obj.C1(1) = Y(6);
			obj.C2(1) = Y(7);
			obj.O1(1) = Y(8);
			obj.O2(1) = Y(9);
			obj.O3(1) = Y(10);
			obj.Cas(1) = Y(11);
			obj.Caf(1) = Y(12);
			obj.Cabls(1) = Y(13);
			obj.Cabhs(1) = Y(14);
			obj.Cablf(1) = Y(15);
			obj.Cabhf(1) = Y(16);
			obj.Rh(1) = Y(17);
			obj.Rhi(1) = Y(18);
			obj.Tr(1) = Y(19);
			obj.PDE(1) = Y(20);
			obj.Ca2(1) = Y(21);
			obj.Cab(1) = Y(22);
			obj.cGMP(1) = Y(23);
        end
		
		function update_time(obj)
			if obj.t+1 == obj.buffer_size
                k = obj.t+1;
				obj.V(1:k-1) = obj.V(2:k);
				obj.mKv(1:k-1) = obj.mKv(2:k); 
				obj.hKv(1:k-1) = obj.hKv(2:k);
				obj.mCa(1:k-1) = obj.mCa(2:k);
				obj.mKCa(1:k-1) = obj.mKCa(2:k);
				obj.C1(1:k-1) = obj.C1(2:k);
				obj.C2(1:k-1) = obj.C2(2:k);
				obj.O1(1:k-1) = obj.O1(2:k);
				obj.O2(1:k-1) = obj.O2(2:k);
				obj.O3(1:k-1) = obj.O3(2:k);
				obj.Cas(1:k-1) = obj.Cas(2:k);
				obj.Caf(1:k-1) = obj.Caf(2:k);
				obj.Cabls(1:k-1) = obj.Cabls(2:k);
				obj.Cabhs(1:k-1) = obj.Cabhs(2:k);
				obj.Cablf(1:k-1) = obj.Cablf(2:k);
				obj.Cabhf(1:k-1) = obj.Cabhf(2:k);
				obj.Rh(1:k-1) = obj.Rh(2:k);
				obj.Rhi(1:k-1) = obj.Rhi(2:k);
				obj.Tr(1:k-1) = obj.Tr(2:k);
				obj.PDE(1:k-1) = obj.PDE(2:k);
				obj.Ca2(1:k-1) = obj.Ca2(2:k);
				obj.Cab(1:k-1) = obj.Cab(2:k);
				obj.cGMP(1:k-1) = obj.cGMP(2:k);
            elseif obj.t+1 < obj.buffer_size
                obj.t = obj.t + 1;
            end
		end
        
        function [Y, c]  = solve(obj,jhv)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
			
            k = obj.t;
			% explicit functions
			%%% Kv %%%
			[am_Kv, bm_Kv, ah_Kv, bh_Kv] = Kv(obj.V(k));
			iKv = obj.gKv*obj.mKv(k)^3*obj.hKv(k)*(obj.V(k)-obj.Ek);
			%%% Ca %%%
			Eca = -12.5*log(obj.Cas(k)/obj.Ca0);
            %Eca = -12.5;
			[am_Ca, bm_Ca] = Ca(obj.V(k));
			hCa = exp((40-obj.V(k))/18)/(1+exp((40-obj.V(k))/18));
			iCa = obj.gCa*obj.mCa(k)^4*hCa*(obj.V(k)-Eca);
			%%% Cl_Ca %%%
			mCl_Ca = 1/(1+exp((0.37-obj.Cas(k))/0.09));
			iCl = obj.gClCa*mCl_Ca*(obj.V(k)-obj.Eclca);
			%%% K_Ca %%%
			[am_KCa, bm_KCa] = K_Ca(obj.V(k));
			mKCas = obj.Cas(k)/(obj.Cas(k)+0.3);
			iKCa = obj.gKCa*obj.mKCa(k)^2*mKCas*(obj.V(k)-obj.Ek);
			%%% h %%%
			[ah, bh] = h(obj.V(k));
			ih = obj.gh*(obj.O1(k)+obj.O2(k)+obj.O3(k))*(obj.V(k)-obj.Eh);
			%%% L %%%
			iL =  obj.gL*(obj.V(k)-obj.EL);
			%%% Cas %%%
			iex = obj.jex*exp(-(obj.V(k)+14)/70)*(obj.Cas(k)-obj.Cae)/(obj.Cas(k)-obj.Cae+2.3);
			iex2 = obj.jex2*(obj.Cas(k)-obj.Cae)/(obj.Cas(k)-obj.Cae+0.5);
			%%% photo %%%
			j = obj.jmax*(obj.cGMP(k))^3/(obj.cGMP(k)^3+10^3);
			iPhoto = -j*(1.0-exp((obj.V(k)-8.5)/17.0));
			
			% values of variables at the next time step
            obj.V(k+1) = obj.V(k)+obj.dt/obj.Cm*-(iKv+iCa+iCl+iKCa+ih+iL+iPhoto+iex+iex2);
            obj.mKv(k+1) = obj.mKv(k)+obj.dt*(am_Kv*(1-obj.mKv(k))-bm_Kv*obj.mKv(k));
            obj.hKv(k+1) = obj.hKv(k)+obj.dt*(ah_Kv*(1-obj.hKv(k))-bh_Kv*obj.hKv(k));
			obj.mCa(k+1) = obj.mCa(k)+obj.dt*(am_Ca*(1-obj.mCa(k))-bm_Ca*obj.mCa(k));
			obj.mKCa(k+1) = obj.mKCa(k)+obj.dt*(am_KCa*(1-obj.mKCa(k))-bm_KCa*obj.mKCa(k));
			obj.C1(k+1) = obj.C1(k)+obj.dt*(-4*ah*obj.C1(k)+bh*obj.C2(k));
			obj.C2(k+1) = obj.C2(k)+obj.dt*(4*ah*obj.C1(k)-3*(ah+bh)*obj.C2(k)+2*bh*obj.O1(k));
			obj.O1(k+1) = obj.O1(k)+obj.dt*(3*ah*obj.C2(k)-(2*ah+2*bh)*obj.O1(k)+3*bh*obj.O2(k));
			obj.O2(k+1) = obj.O2(k)+obj.dt*(2*ah*obj.O1(k)-(ah+3*bh)*obj.O2(k)+4*bh*obj.O3(k));
			obj.O3(k+1) = obj.O3(k)+obj.dt*(ah*obj.O2(k)-4*bh*obj.O3(k));
			obj.Cas(k+1) = obj.Cas(k)+obj.dt*(-(iCa+iex+iex2)/(2*obj.F*obj.V1)*10^(-6)...
                -obj.DCa*obj.S1/(obj.delta*obj.V1)*(obj.Cas(k)-obj.Caf(k))...
                -obj.Lb1*obj.Cas(k)*(obj.BL-obj.Cabls(k))+obj.Lb2*obj.Cabls(k)...
                -obj.Hb1*obj.Cas(k)*(obj.BH-obj.Cabhs(k))+obj.Hb2*obj.Cabhs(k));
			obj.Caf(k+1) = obj.Caf(k)+obj.dt*(obj.DCa*obj.S1/(obj.delta*obj.V2)*(obj.Cas(k)-...
                obj.Caf(k))-obj.Lb1*obj.Caf(k)*(obj.BL-obj.Cablf(k))+...
                obj.Lb2*obj.Cablf(k)-obj.Hb1*obj.Caf(k)*(obj.BH-obj.Cabhf(k))...
                +obj.Hb2*obj.Cabhf(k));
			obj.Cabls(k+1) = obj.Cabls(k)+obj.dt*(obj.Lb1*obj.Cas(k)*...
                (obj.BL-obj.Cabls(k))-obj.Lb2*obj.Cabls(k));
			obj.Cabhs(k+1) = obj.Cabhs(k)+obj.dt*(obj.Hb1*obj.Cas(k)*...
                (obj.BH-obj.Cabhs(k))-obj.Hb2*obj.Cabhs(k));
			obj.Cablf(k+1) = obj.Cablf(k)+obj.dt*(obj.Lb1*obj.Caf(k)*...
                (obj.BL-obj.Cablf(k))-obj.Lb2*obj.Cablf(k));
			obj.Cabhf(k+1) = obj.Cabhf(k)+obj.dt*(obj.Hb1*obj.Caf(k)*...
                (obj.BH-obj.Cabhf(k))-obj.Hb2*obj.Cabhf(k));
			obj.Rh(k+1) = obj.Rh(k)+obj.dt*(jhv-obj.a1*obj.Rh(k)+obj.a2*obj.Rhi(k));
			obj.Rhi(k+1) = obj.Rhi(k)+obj.dt*(obj.a1*obj.Rh(k)-(obj.a2+obj.a3)*obj.Rhi(k));
			obj.Tr(k+1) = obj.Tr(k)+obj.dt*(obj.e*obj.Rh(k)*(obj.Ttot-obj.Tr(k))-...
                obj.b1*obj.Tr(k)-obj.tau1*obj.Tr(k)*(obj.PDEtot-obj.PDE(k))+obj.tau2*obj.PDE(k));
			obj.PDE(k+1) = obj.PDE(k)+obj.dt*(obj.tau1*obj.Tr(k)*...
                (obj.PDEtot-obj.PDE(k))-obj.tau2*obj.PDE(k));
			obj.Ca2(k+1) = obj.Ca2(k)+obj.dt*(obj.b*j-obj.gammaCa*(obj.Ca2(k)-obj.Ca0p)...
                -obj.k1*(obj.eT-obj.Cab(k))*obj.Ca2(k)+obj.k2*obj.Cab(k));
			obj.Cab(k+1) = obj.Cab(k)+obj.dt*(obj.k1*(obj.eT-obj.Cab(k))*...
                obj.Ca2(k)-obj.k2*obj.Cab(k));
			obj.cGMP(k+1) = obj.cGMP(k)+obj.dt*(obj.Amax/(1+(obj.Ca2(k)/obj.Kc)^4)...
                -obj.cGMP(k)*(obj.Vbar+obj.sigma*obj.PDE(k)));
            
            Y = [obj.V(k+1), obj.mKv(k+1), obj. hKv(k+1), obj.mCa(k+1), ...
                obj.mKCa(k+1), obj.C1(k+1), obj.C2(k+1), obj.O1(k+1), ...
                obj.O2(k+1), obj.O3(k+1), obj.Cas(k+1), obj.Caf(k+1), ...
                obj.Cabls(k+1), obj.Cabhs(k+1), obj.Cablf(k+1), ...
                obj.Cabhf(k+1), obj.Rh(k+1), obj.Rhi(k+1), obj.Tr(k+1), ...
                obj.PDE(k+1), obj.Ca2(k+1), obj.Cab(k+1), obj.cGMP(k+1)];
            c = [am_Kv, bm_Kv, ah_Kv, bh_Kv, iKv, Eca, am_Ca, bm_Ca, hCa, iCa,...
                mCl_Ca, iCl, am_KCa, bm_KCa, mKCas, iKCa, ah, bh, ih, iL, ...
                iex, iex2, j, iPhoto];
            
        end
        
        function v = get_V(obj)
            v = obj.V;
        end
        
        function cgmp = get_cGMP(obj)
            cgmp = obj.cGMP;
        end
        
        function cas = get_Cas(obj)
            cas = obj.Cas;
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


