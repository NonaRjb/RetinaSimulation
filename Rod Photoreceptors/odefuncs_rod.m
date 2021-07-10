function dY = odefuncs_rod(t, Y, jhvt, jhv)

% input current of light photons 
jhv = interp1(jhvt,jhv,t); % Interpolate the data set (gt,g) at time t

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
V = Y(1);
mKv = Y(2); 
hKv = Y(3);
mCa = Y(4);
mKCa = Y(5);
C1 = Y(6);
C2 = Y(7);
O1 = Y(8);
O2 = Y(9);
O3 = Y(10);
Cas = Y(11);
Caf = Y(12);
Cabls = Y(13);
Cabhs = Y(14);
Cablf = Y(15);
Cabhf = Y(16);
Rh = Y(17);
Rhi = Y(18);
Tr = Y(19);
PDE = Y(20);
Ca2 = Y(21);
Cab = Y(22);
cGMP = Y(23);


% explicit functions
%%% Kv %%%
[am_Kv, bm_Kv, ah_Kv, bh_Kv] = Kv(V);
iKv = gKv*mKv^3*hKv*(V-Ek);
%%% Ca %%%
Eca = -12.5*log(Cas/Ca0);
[am_Ca, bm_Ca] = Ca(V);
hCa = exp((40-V)/18)/(1+exp((40-V)/18));
iCa = gCa*mCa^4*hCa*(V-Eca);
%%% Cl_Ca %%%
mCl_Ca = 1/(1+exp((0.37-Cas)/0.09));
iCl = gClCa*mCl_Ca*(V-Eclca);
%%% K_Ca %%%
[am_KCa, bm_KCa] = K_Ca(V);
mKCas = Cas/(Cas+0.3);
iKCa = gKCa*mKCa^2*mKCas*(V-Ek);
%%% h %%%
[ah, bh] = h(V);
ih = gh*(O1+O2+O3)*(V-Eh);
%%% L %%%
iL =  gL*(V-EL);
%%% Cas %%%
iex = jex*exp(-(V+14)/70)*(Cas-Cae)/(Cas-Cae+2.3);
iex2 = jex2*(Cas-Cae)/(Cas-Cae+0.5);
%%% photo %%%
j = jmax*(cGMP)^3/(cGMP^3+10^3);
iPhoto = -j*(1.0-exp((V-8.5)/17.0));

% differential equations
dY = zeros(23,1);
dY(1) = 1/Cm*-(iKv+iCa+iCl+iKCa+ih+iL+iPhoto+iex+iex2);
dY(2) = am_Kv*(1-mKv)-bm_Kv*mKv;
dY(3) = ah_Kv*(1-hKv)-bh_Kv*hKv;
dY(4) = am_Ca*(1-mCa)-bm_Ca*mCa;
dY(5) = am_KCa*(1-mKCa)-bm_KCa*mKCa;
dY(6) = -4*ah*C1+bh*C2;
dY(7) = 4*ah*C1-3*(ah+bh)*C2+2*bh*O1;
dY(8) = 3*ah*C2-(2*ah+2*bh)*O1+3*bh*O2;
dY(9) = 2*ah*O1-(ah+3*bh)*O2+4*bh*O3;
dY(10) = ah*O2-4*bh*O3;
dY(11) = -(iCa+iex+iex2)/(2*F*V1)*10^(-6)-DCa*S1/(delta*V1)*(Cas-Caf)-Lb1*Cas*(BL-Cabls)+Lb2*Cabls-Hb1*Cas*(BH-Cabhs)+Hb2*Cabhs;
dY(12) = DCa*S1/(delta*V2)*(Cas-Caf)-Lb1*Caf*(BL-Cablf)+Lb2*Cablf-Hb1*Caf*(BH-Cabhf)+Hb2*Cabhf;
dY(13) = Lb1*Cas*(BL-Cabls)-Lb2*Cabls;
dY(14) = Hb1*Cas*(BH-Cabhs)-Hb2*Cabhs;
dY(15) = Lb1*Caf*(BL-Cablf)-Lb2*Cablf;
dY(16) = Hb1*Caf*(BH-Cabhf)-Hb2*Cabhf;
dY(17) = jhv-a1*Rh+a2*Rhi;
dY(18) = a1*Rh-(a2+a3)*Rhi;
dY(19) = e*Rh*(Ttot-Tr)-b1*Tr-tau1*Tr*(PDEtot-PDE)+tau2*PDE;
dY(20) = tau1*Tr*(PDEtot-PDE)-tau2*PDE;
dY(21) = b*j-gammaCa*(Ca2-Ca0p)-k1*(eT-Cab)*Ca2+k2*Cab;
dY(22) = k1*(eT-Cab)*Ca2-k2*Cab;
dY(23) = Amax/(1+(Ca2/Kc)^4)-cGMP*(Vbar+sigma*PDE);

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
