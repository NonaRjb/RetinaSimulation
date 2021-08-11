function dY = odefuncs_retina(t, Y, jhvt, jhv)

% input current of light photons 
jhv = interp1(jhvt,jhv,t); % Interpolate the data set (gt,g) at time t

% constants rod
rod_Cm = 0.02;      % [nF]
rod_gKv = 2.0;      % [nS]
rod_Ek = -74.0;       % [mV]
rod_gCa = 0.7;      % [nS]
rod_Ca0 = 1600;     % [uM]
rod_gClCa = 2.0;    % [nS]
rod_Eclca = -20.0;    % [mV]
rod_gKCa = 5.0;     % [nS]
rod_gh = 3.0;       % [nS]
rod_Eh = -32.0;       % [mV]
rod_gL = 0.35;      % [nS]
rod_EL = -77.0;       % [mV]
rod_F = 9.64846*10^4;         % [cmol^(-1)] Faraday const.
rod_V1 = 3.812*10^(-13);    % [dm^3] Volume of submembrane area
rod_V2 = 5.236*10^(-13);    % [dm^3] Volume of deep intracellular area
rod_DCa = 6*10^(-8);        % [d(m^2)(s^(-1))] Ca diffusion coefficient 
rod_S1 = 3.142*10^(-8);     % [dm^2] Surface area of the submembrane and the deep intracellular area spherical boundary
rod_delta = 5.9*10^(-5);    % [dm] Distance between submembrane area and the deep intracellular area
rod_Lb1 = 0.4;              % [s^(-1)uM^(-1)] On rate constant for the binding of Ca to low-affinity buffer
rod_Lb2 = 0.2;              % [s^(-1)uM^(-1)] Off rate constant for the binding of Ca to low-affinity buffer
rod_Hb1 = 100;              % [s^(-1)uM^(-1)]  On rate constant for the binding of Ca to high-affinity buffer
rod_Hb2 = 90;               % [s^(-1)uM^(-1)]  off rate constants for the binding of Ca to high-affinity buffer
rod_BL = 500;               % [uM] Total low-affinity buffer concentration
rod_BH = 200;               % [uM] Total high-affinity buffer concentration
rod_jex = 20;               % [pA] Maximum Na-Ca exchanger current
rod_jex2 = 20;              % [pA] Maximum Ca-ATPase exchanger current
rod_Cae = 0.01;             % [uM] Minimum intracellular Ca2+ concentration
rod_a1 = 50.0;                % [s^(-1)] Rate constant of Rh* inactivation
rod_a2 = 0.0003;            % [s^(-1)] Rate constant of the reaction Rhi -> Rh*
rod_a3 = 0.03;              % [s^(-1)] Rate constant of the decay of inactive rhodopsin
rod_e = 0.5;                % [s^(-1)uM^(-1)] Rate constant of T* activation
rod_Ttot = 1000;            % [uM] Total transducin
rod_b1 = 2.5;               % [s^(-1)] Rate constant of T* inactivation
rod_tau1 = 0.2;             % [s^(-1)uM^(-1)] Rate constant of PDE activation
rod_tau2 = 5;               % [s^(-1)] Rate constant of PDE inactivation
rod_PDEtot = 100;           % [uM] Phosphodiasterase
rod_gammaCa = 50;           % [s^(-1)] Rate constant of Ca2+ extrusion in the absence of Ca2+ buffers mediated by the Na+?Ca2+ exchanger
rod_Ca0p = 0.1;             % [uM] Intracellular Ca2+ concentration at the steady state
rod_b = 0.25;               % [uMs^(-1)pA^(-1)] Proportionality constant between Ca2+ influx and photocurrent
rod_k1 = 0.2;               % [uM^(-1)] on rate constants for the binding of Ca2+ to the buffer
rod_k2 = 0.8;               % [uM^(-1)] off rate constants for the binding of Ca2+ to the buffer
rod_eT = 500;               % [uM] Low-affinity Ca2+ buffer concentration
rod_Vbar = 0.4;             % [s^(-1)] Cyclic GMP hydrolysis in dark
rod_Kc = 0.1;               % [uM] Maximal activity of guanylate cyclase
rod_Amax = 65.6;            % [uMs^(-1)] Maximal activity of guanylate cyclase
rod_sigma = 1.0;            % [s^(-1)uM^(-1)] Proportionality constant
rod_jmax = 5040;            % [pA] Maximal cyclic GMP gated current inexcised patches

% constants rod-bipolar synapse
rbp_syn_gmax = 2.56;    % [nS] Maximum synapse conductance
rbp_syn_Esyn = 0.0;       % [mV] Synapse's reversal potential
rbp_syn_tau = 0.01;       % [ms] Time constant
rbp_syn_Vslope = 20;    % [mV] Voltage sensitivity of the synapse
rbp_syn_Vth = -36.185963;      % [mV]

% constants bipolar
bp_C = 0.01;   % [pF]
bp_gKv = 2.0;  % [nS]
bp_Ek = -58;    % [mV]
bp_gA = 50;    % [nS]
bp_gh = 0.975; % [nS]
bp_Eh = -17.7; % [mV]
bp_gCa = 1.1;  % [nS]
bp_Ca0 = 2500; % [uM]
bp_gKCa = 8.5; % [nS]
bp_gL = 0.23;  % [nS]
bp_El = -21.0;   % [mV]
bp_F = 9.649e04;   % [cmol^(-1)] Faraday constant
bp_DCa = 6e-08;    % [dm^2*sec^(-1)] Ca diffusion coefficient 
bp_Vs = 1.692e-13; % [dm^(-3)] Volume of submembrane area 
bp_Vd = 7.356e-13; % [dm^(-3)] Volume of the deep intracellular area
bp_Ssd = 4e-08;    % [dm^(-2)] Surface area of the submembrane and the deep intracellular area spherical boundary
bp_dsd = 5.9e-05;   % [dm] Distance between submembrane area and the deep intracellular area
bp_Cablmax_s = 300;  % [uM] Total low-affinity buffer concentration
bp_Cabhmax_s = 100;  % [uM] Total high-affinity buffer concentration
bp_Cablmax_d = 500;  % [uM] Total low-affinity buffer concentration
bp_Cabhmax_d = 300;  % [uM] Total high-affinity buffer concentration
bp_abl = 0.4;      % [sec^(-1)uM^(-1)] On rate constants for the binding of Ca to low-affinity buffer
bp_bbl = 0.2;      % [sec^(-1)uM^(-1)] Off rate constants for the binding of Ca to low-affinity buffer
bp_abh = 100;      % [sec^(-1)uM^(-1)] On rate constants for the binding of Ca to high-affinity buffer
bp_bbh = 90;       % [sec^(-1)uM^(-1)] Off rate constants for the binding of Ca to high-affinity buffer
bp_Jex = 9.0*40;        % [pA] Maximum Na-Ca exchanger current
bp_Jex2 = 9.5*40;     % [pA] Maximum Ca-ATPase exchanger current
bp_Camin = 0.01;   % [uM] Minimum intracellular Ca concentration for Ca extrusion

% variables rod
rod_V = Y(1);
rod_mKv = Y(2); 
rod_hKv = Y(3);
rod_mCa = Y(4);
rod_mKCa = Y(5);
rod_C1 = Y(6);
rod_C2 = Y(7);
rod_O1 = Y(8);
rod_O2 = Y(9);
rod_O3 = Y(10);
rod_Cas = Y(11);
rod_Caf = Y(12);
rod_Cabls = Y(13);
rod_Cabhs = Y(14);
rod_Cablf = Y(15);
rod_Cabhf = Y(16);
rod_Rh = Y(17);
rod_Rhi = Y(18);
rod_Tr = Y(19);
rod_PDE = Y(20);
rod_Ca2 = Y(21);
rod_Cab = Y(22);
rod_cGMP = Y(23);

% variables rod-bipolar synapse
rbp_syn_S = Y(24);

% variables bipolar
bp_V = Y(25);
bp_mKv = Y(26);
bp_hKv = Y(27);
bp_mA = Y(28);
bp_hA = Y(29);
bp_C1 = Y(30);
bp_C2 = Y(31);
bp_O1 = Y(32);
bp_O2 = Y(33);
bp_O3 = Y(34);
bp_mCa = Y(35);
bp_mKCa = Y(36);
bp_Cas = Y(37);
bp_Cad = Y(38);
bp_Cabls = Y(39);
bp_Cabhs = Y(40);
bp_Cabld = Y(41);
bp_Cabhd = Y(42);

% explicit functions rod
%%% Kv %%%
[rod_am_Kv, rod_bm_Kv, rod_ah_Kv, rod_bh_Kv] = rod_Kv(rod_V);
rod_iKv = rod_gKv*rod_mKv^3*rod_hKv*(rod_V-rod_Ek);
%%% Ca %%%
rod_Eca = -12.5*log(rod_Cas/rod_Ca0);
[rod_am_Ca, rod_bm_Ca] = rod_Ca(rod_V);
rod_hCa = exp((40-rod_V)/18)/(1+exp((40-rod_V)/18));
rod_iCa = rod_gCa*rod_mCa^4*rod_hCa*(rod_V-rod_Eca);
%%% Cl_Ca %%%
rod_mCl_Ca = 1.0/(1+exp((0.37-rod_Cas)/0.09));
rod_iCl = rod_gClCa*rod_mCl_Ca*(rod_V-rod_Eclca);
%%% K_Ca %%%
[rod_am_KCa, rod_bm_KCa] = rod_K_Ca(rod_V);
rod_mKCas = rod_Cas/(rod_Cas+0.3);
rod_iKCa = rod_gKCa*rod_mKCa^2*rod_mKCas*(rod_V-rod_Ek);
%%% h %%%
[rod_ah, rod_bh] = rod_h(rod_V);
rod_ih = rod_gh*(rod_O1+rod_O2+rod_O3)*(rod_V-rod_Eh);
%%% L %%%
rod_iL = rod_gL*(rod_V-rod_EL);
%%% Cas %%%
rod_iex = rod_jex*exp(-(rod_V+14)/70)*(rod_Cas-rod_Cae)/(rod_Cas-rod_Cae+2.3);
rod_iex2 = rod_jex2*(rod_Cas-rod_Cae)/(rod_Cas-rod_Cae+0.5);
%%% photo %%%
rod_j = rod_jmax*(rod_cGMP)^3/(rod_cGMP^3+10^3);
rod_iPhoto = -rod_j*(1.0-exp((rod_V-8.5)/17.0));


% explicit functions bipolar
%%% Kv %%%
[bp_am_Kv, bp_bm_Kv, bp_ah_Kv, bp_bh_Kv] = bp_Kv(bp_V);
bp_iKv = bp_gKv*bp_mKv^3*bp_hKv*(bp_V-bp_Ek);
%%% A %%%
[bp_am_A, bp_bm_A, bp_ah_A, bp_bh_A] = bp_A(bp_V);
bp_iA = bp_gA*bp_mA^3*bp_hA*(bp_V-bp_Ek);
%%% h %%%
[bp_ah, bp_bh] = bp_h(bp_V);
bp_ih = bp_gh*(bp_O1+bp_O2+bp_O3)*(bp_V-bp_Eh);
%%% Ca %%%
bp_Eca = 12.9*log(bp_Ca0/bp_Cas);
[bp_am_Ca, bp_bm_Ca] = bp_Ca(bp_V);
bp_hCa = exp(-(bp_V-50)/11)/(exp(-(bp_V-50)/11)+1.0);
bp_iCa = bp_gCa*bp_mCa^4*bp_hCa*(bp_V-bp_Eca);
%%% K_Ca %%%
[bp_am_KCa, bp_bm_KCa] = bp_K_Ca(bp_V);
bp_mKc1 = bp_Cas/(bp_Cas+0.2);
bp_iKCa = bp_gKCa*bp_mKCa^2*bp_mKc1*(bp_V-bp_Ek);
%%% L %%%
bp_iL = bp_gL*(bp_V-bp_El);
%%% Iex %%%
bp_Iex = bp_Jex*(bp_Cas-bp_Camin)/(bp_Cas-bp_Camin+2.3)*exp(-(bp_V+14.0)/70.0);
bp_Iex2 = bp_Jex2*(bp_Cas-bp_Camin)/(bp_Cas-bp_Camin+0.5);
%%% Isyn %%%
% if rod_V > rbp_syn_Vth
%    rbp_syn_S_inf = tanh((rod_V-rbp_syn_Vth)/rbp_syn_Vslope);
% else
%    rbp_syn_S_inf = 0;
% end
rbp_syn_S_inf = tanh(abs(rod_V-rbp_syn_Vth)/rbp_syn_Vslope);
rbp_Isyn = rbp_syn_gmax*rbp_syn_S*(bp_V-rbp_syn_Esyn);

dY = zeros(42,1);

% differential equations rod
dY(1) = 1/rod_Cm*-(rod_iKv+rod_iCa+rod_iCl+rod_iKCa+rod_ih+rod_iL+...
    rod_iPhoto+rod_iex+rod_iex2);
dY(2) = rod_am_Kv*(1-rod_mKv)-rod_bm_Kv*rod_mKv;
dY(3) = rod_ah_Kv*(1-rod_hKv)-rod_bh_Kv*rod_hKv;
dY(4) = rod_am_Ca*(1-rod_mCa)-rod_bm_Ca*rod_mCa;
dY(5) = rod_am_KCa*(1-rod_mKCa)-rod_bm_KCa*rod_mKCa;
dY(6) = -4*rod_ah*rod_C1+rod_bh*rod_C2;
dY(7) = 4*rod_ah*rod_C1-(3*rod_ah+rod_bh)*rod_C2+2*rod_bh*rod_O1;
dY(8) = 3*rod_ah*rod_C2-(2*rod_ah+2*rod_bh)*rod_O1+3*rod_bh*rod_O2;
dY(9) = 2*rod_ah*rod_O1-(rod_ah+3*rod_bh)*rod_O2+4*rod_bh*rod_O3;
dY(10) = rod_ah*rod_O2-4*rod_bh*rod_O3;
dY(11) = -(rod_iCa+rod_iex+rod_iex2)/(2.0*rod_F*rod_V1)*10^(-6)-...
    rod_DCa*rod_S1/(rod_delta*rod_V1)*(rod_Cas-rod_Caf)-rod_Lb1*rod_Cas*...
    (rod_BL-rod_Cabls)+rod_Lb2*rod_Cabls-rod_Hb1*rod_Cas*...
    (rod_BH-rod_Cabhs)+rod_Hb2*rod_Cabhs;
dY(12) = rod_DCa*rod_S1/(rod_delta*rod_V2)*(rod_Cas-rod_Caf)-...
    rod_Lb1*rod_Caf*(rod_BL-rod_Cablf)+rod_Lb2*rod_Cablf-...
    rod_Hb1*rod_Caf*(rod_BH-rod_Cabhf)+rod_Hb2*rod_Cabhf;
dY(13) = rod_Lb1*rod_Cas*(rod_BL-rod_Cabls)-rod_Lb2*rod_Cabls;
dY(14) = rod_Hb1*rod_Cas*(rod_BH-rod_Cabhs)-rod_Hb2*rod_Cabhs;
dY(15) = rod_Lb1*rod_Caf*(rod_BL-rod_Cablf)-rod_Lb2*rod_Cablf;
dY(16) = rod_Hb1*rod_Caf*(rod_BH-rod_Cabhf)-rod_Hb2*rod_Cabhf;
dY(17) = jhv-rod_a1*rod_Rh+rod_a2*rod_Rhi;
dY(18) = rod_a1*rod_Rh-(rod_a2+rod_a3)*rod_Rhi;
dY(19) = rod_e*rod_Rh*(rod_Ttot-rod_Tr)-rod_b1*rod_Tr-rod_tau1*rod_Tr*...
    (rod_PDEtot-rod_PDE)+rod_tau2*rod_PDE;
dY(20) = rod_tau1*rod_Tr*(rod_PDEtot-rod_PDE)-rod_tau2*rod_PDE;
dY(21) = rod_b*rod_j-rod_gammaCa*(rod_Ca2-rod_Ca0p)-rod_k1*...
    (rod_eT-rod_Cab)*rod_Ca2+rod_k2*rod_Cab;
dY(22) = rod_k1*(rod_eT-rod_Cab)*rod_Ca2-rod_k2*rod_Cab;
dY(23) = rod_Amax/(1.0+rod_Ca2^4/rod_Kc^4)-rod_cGMP*(rod_Vbar+rod_sigma*rod_PDE);

% differential equations rod-bipolar synapse
dY(24) = 0;
if rod_V >= rbp_syn_Vth
	dY(24) = 0;
elseif abs(rod_V-rbp_syn_Vth) >= 0
	dY(24) = (rbp_syn_S_inf-rbp_syn_S)/((1-rbp_syn_S_inf)*rbp_syn_tau);
else
	dY(24) = 0;
end
% dY(24) = (rbp_syn_S_inf-rbp_syn_S)/((1-rbp_syn_S_inf)*rbp_syn_tau);

% differential equations bipolar
dY(25) = 1/bp_C*-(bp_iKv+bp_iA+bp_ih+bp_iCa+bp_iKCa+bp_iL+rbp_Isyn);
dY(26) = bp_am_Kv*(1.0-bp_mKv)-bp_bm_Kv*bp_mKv;
dY(27) = bp_ah_Kv*(1.0-bp_hKv)-bp_bh_Kv*bp_hKv;
dY(28) = bp_am_A*(1.0-bp_mA)-bp_bm_A*bp_mA;
dY(29) = bp_ah_A*(1.0-bp_hA)-bp_bh_A*bp_hA;
dY(30) = -4.0*bp_ah*bp_C1+bp_bh*bp_C2;
dY(31) = 4.0*bp_ah*bp_C1-(3*bp_ah+bp_bh)*bp_C2+2.0*bp_bh*bp_O1;
dY(32) = 3.0*bp_ah*bp_C2-2.0*(bp_ah+bp_bh)*bp_O1+3.0*bp_bh*bp_O2;
dY(33) = 2.0*bp_ah*bp_O1-(bp_ah+3*bp_bh)*bp_O2+4*bp_bh*bp_O3;
dY(34) = bp_ah*bp_O2-4*bp_bh*bp_O3;
dY(35) = bp_am_Ca*(1.0-bp_mCa)-bp_bm_Ca*bp_mCa;
dY(36) = bp_am_KCa*(1.0-bp_mKCa)-bp_bm_KCa*bp_mKCa;
dY(37) = -(bp_iCa+bp_Iex+bp_Iex2)/(2*bp_F*bp_Vs)*10^(-6)-bp_DCa*bp_Ssd/(bp_Vs*bp_dsd)*...
    (bp_Cas-bp_Cad)+bp_bbl*bp_Cabls-bp_abl*bp_Cas*(bp_Cablmax_s-bp_Cabls)+...
    bp_bbh*bp_Cabhs-bp_abh*bp_Cas*(bp_Cabhmax_s-bp_Cabhs);
dY(38) = bp_DCa*bp_Ssd/(bp_Vd*bp_dsd)*(bp_Cas-bp_Cad)+...
    bp_bbl*bp_Cabld-bp_abl*bp_Cad*(bp_Cablmax_d-bp_Cabld)+...
    bp_bbh*bp_Cabhd-bp_abh*bp_Cad*(bp_Cabhmax_d-bp_Cabhd);
dY(39) = bp_abl*bp_Cas*(bp_Cablmax_s-bp_Cabls)-bp_bbl*bp_Cabls;
dY(40) = bp_abh*bp_Cas*(bp_Cabhmax_s-bp_Cabhs)-bp_bbh*bp_Cabhs;
dY(41) = bp_abl*bp_Cad*(bp_Cablmax_d-bp_Cabld)-bp_bbl*bp_Cabld;
dY(42) = bp_abh*bp_Cad*(bp_Cabhmax_d-bp_Cabhd)-bp_bbh*bp_Cabhd;


end

%%%%% Rod Kv %%%%%
function [am_Kv, bm_Kv, ah_Kv, bh_Kv] = rod_Kv(V)

am_Kv = 5*(100-V)/(exp((100-V)/42)-1);
bm_Kv = 9*exp(-(V-20)/40);
ah_Kv = 0.15*exp(-V/22);
bh_Kv = 0.4125/(exp((10-V)/7)+1);

end
%%%%%%%%%%%%%%%%

%%%%% Rod Ca %%%%%
function [am_Ca, bm_Ca] = rod_Ca(V)

am_Ca = 3*(80-V)/(exp((80-V)/25)-1);
bm_Ca = 10/(1+exp((V+38)/7));

end
%%%%%%%%%%%%%%%%

%%%%% Rod K_Ca %%%%%
function [am_KCa, bm_KCa] = rod_K_Ca(V)

am_KCa = 15*(80-V)/(exp((80-V)/40)-1);
bm_KCa = 20*exp(-V/35);

end
%%%%%%%%%%%%%%%%

%%%%% Rod h %%%%%
function [ah, bh] = rod_h(V)

ah = 8/(exp((V+78)/14)+1);
bh = 18/exp(-(V+8)/19+1);

end
%%%%%%%%%%%%%%%%


%%%%% Bipolar Kv %%%%%
function [am_Kv, bm_Kv, ah_Kv, bh_Kv] = bp_Kv(V)

am_Kv = 400/(exp(-(V-15)/36)+1.0);
bm_Kv = 1.0*exp(-V/13);
ah_Kv = 0.0003*exp(-V/7);
bh_Kv = 80/(exp((V+115)/15)+1.0)+0.02;

end
%%%%%%%%%%%%%%%%

%%%%% Bipolar A %%%%%
function [am_A, bm_A, ah_A, bh_A] = bp_A(V)

am_A = 2400/(exp(-(V-50.0)/28.0)+1.0);
bm_A = 12*exp(-V/10.0);
ah_A = 0.045*exp(-V/13);
bh_A = 75/(exp(-(V+30)/15)+1.0);

end
%%%%%%%%%%%%%%%%

%%%%% Bipolar h %%%%%
function [ah, bh] = bp_h(V)

ah = 3/(exp((V+110)/15)+1.0);
bh = 1.5/(exp(-(V+115)/15)+1.0);

end
%%%%%%%%%%%%%%%%


%%%%% Bipolar Ca %%%%%
function [am_Ca, bm_Ca] = bp_Ca(V)

am_Ca = 12000*(120-V)/(exp((120-V)/25)-1.0);
bm_Ca = 40000/(exp((V+68)/25)+1.0);

end
%%%%%%%%%%%%%%%%


%%%%% Bipolar K_Ca %%%%%
function [am_KCa, bm_KCa] = bp_K_Ca(V)

am_KCa = 100*(230-V)/(exp((230-V)/52)-1.0);
bm_KCa = 120*exp(-V/95);

end
%%%%%%%%%%%%%%%%
