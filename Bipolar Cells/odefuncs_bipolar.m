function dY = odefuncs_bipolar(t, Y)

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
bp_dsd = 5.9e05;   % [dm] Distance between submembrane area and the deep intracellular area
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

% variables bipolar
bp_V = Y(1);
bp_mKv = Y(2);
bp_hKv = Y(3);
bp_mA = Y(4);
bp_hA = Y(5);
bp_C1 = Y(6);
bp_C2 = Y(7);
bp_O1 = Y(8);
bp_O2 = Y(9);
bp_O3 = Y(10);
bp_mCa = Y(11);
bp_mKCa = Y(12);
bp_Cas = Y(13);
bp_Cad = Y(14);
bp_Cabls = Y(15);
bp_Cabhs = Y(16);
bp_Cabld = Y(17);
bp_Cabhd = Y(18);

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
% bp_Eca = 12.9*log(bp_Ca0/0.011561930331728);
[bp_am_Ca, bp_bm_Ca] = bp_Ca(bp_V);
bp_hCa = exp(-(bp_V-50)/11)/(exp(-(bp_V-50)/11)+1.0);
% mca_4 = cast(bp_mCa, 'double')^cast(4, 'double');
bp_iCa = bp_gCa*bp_mCa^4*bp_hCa*(bp_V-bp_Eca);
%%% K_Ca %%%
[bp_am_KCa, bp_bm_KCa] = bp_K_Ca(bp_V);
bp_mKc1 = bp_Cas/(bp_Cas+0.2);
% bp_mKc1 = 0.011561930331728/(0.011561930331728+0.2);
bp_iKCa = bp_gKCa*bp_mKCa^2*bp_mKc1*(bp_V-bp_Ek);
%%% L %%%
bp_iL = bp_gL*(bp_V-bp_El);
%%% Iex %%%
bp_Iex = bp_Jex*(bp_Cas-bp_Camin)/(bp_Cas-bp_Camin+2.3)*exp(-(bp_V+14.0)/70.0);
bp_Iex2 = bp_Jex2*(bp_Cas-bp_Camin)/(bp_Cas-bp_Camin+0.5);
%%% Isyn %%%
%if rod_V > rbp_syn_Vth
%    rbp_syn_S_inf = tanh((rod_V-rbp_syn_Vth)/rbp_syn_Vslope);
%else
%    rbp_syn_S_inf = 0;
%end
% rbp_syn_S_inf = tanh(abs(rod_V-rbp_syn_Vth)/rbp_syn_Vslope);
% rbp_Isyn = rbp_syn_gmax*rbp_syn_S*(bp_V-rbp_syn_Esyn);

dY = zeros(18, 1);
% differential equations bipolar
dY(1) = 1/bp_C*-(bp_iKv+bp_iA+bp_ih+bp_iCa+bp_iKCa+bp_iL);
dY(2) = bp_am_Kv*(1.0-bp_mKv)-bp_bm_Kv*bp_mKv;
dY(3) = bp_ah_Kv*(1.0-bp_hKv)-bp_bh_Kv*bp_hKv;
dY(4) = bp_am_A*(1.0-bp_mA)-bp_bm_A*bp_mA;
dY(5) = bp_ah_A*(1.0-bp_hA)-bp_bh_A*bp_hA;
dY(6) = -4.0*bp_ah*bp_C1+bp_bh*bp_C2;
dY(7) = 4.0*bp_ah*bp_C1-(3*bp_ah+bp_bh)*bp_C2+2.0*bp_bh*bp_O1;
dY(8) = 3.0*bp_ah*bp_C2-2.0*(bp_ah+bp_bh)*bp_O1+3.0*bp_bh*bp_O2;
dY(9) = 2.0*bp_ah*bp_O1-(bp_ah+3*bp_bh)*bp_O2+4*bp_bh*bp_O3;
dY(10) = bp_ah*bp_O2-4*bp_bh*bp_O3;
dY(11) = bp_am_Ca*(1.0-bp_mCa)-bp_bm_Ca*bp_mCa;
dY(12) = bp_am_KCa*(1.0-bp_mKCa)-bp_bm_KCa*bp_mKCa;
dY(13) = -(bp_iCa+bp_Iex+bp_Iex2)/(2*bp_F*bp_Vs)*10^(-6)-bp_DCa*bp_Ssd/(bp_Vs*bp_dsd)*...
    (bp_Cas-bp_Cad)+bp_bbl*bp_Cabls-bp_abl*bp_Cas*(bp_Cablmax_s-bp_Cabls)+...
    bp_bbh*bp_Cabhs-bp_abh*bp_Cas*(bp_Cabhmax_s-bp_Cabhs);
dY(14) = bp_DCa*bp_Ssd/(bp_Vd*bp_dsd)*(bp_Cas-bp_Cad)+...
    bp_bbl*bp_Cabld-bp_abl*bp_Cad*(bp_Cablmax_d-bp_Cabld)+...
    bp_bbh*bp_Cabhd-bp_abh*bp_Cad*(bp_Cabhmax_d-bp_Cabhd);
dY(15) = bp_abl*bp_Cas*(bp_Cablmax_s-bp_Cabls)-bp_bbl*bp_Cabls;
dY(16) = bp_abh*bp_Cas*(bp_Cabhmax_s-bp_Cabhs)-bp_bbh*bp_Cabhs;
dY(17) = bp_abl*bp_Cad*(bp_Cablmax_d-bp_Cabld)-bp_bbl*bp_Cabld;
dY(18) = bp_abh*bp_Cad*(bp_Cabhmax_d-bp_Cabhd)-bp_bbh*bp_Cabhd;
% dY(14) = 0;
% dY(15) = 0;
%dY(16) = 0;
%dY(17) = 0;
% dY(18) = 0;

end


%%%%% Bipolar Kv %%%%%
function [am_Kv, bm_Kv, ah_Kv, bh_Kv] = bp_Kv(V)

am_Kv = 400/(exp(-(V-15)/36)+1.0);
bm_Kv = 1.0*exp(-V/13);
ah_Kv = 0.0003*exp(-V/7);
bh_Kv = 80/(exp(-(V+115)/15)+1.0)+0.02;

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
