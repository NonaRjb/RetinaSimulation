function dY = odefuncs_bipolar(t, Y)

% constants
C = 0.02;   % [pF]
gKv = 2.0;  % [nS]
Ek = 58;    % [mV]
gA = 35;    % [nS]
gh = 0.975; % [nS]
Eh = -17.7; % [mV]
gCa = 1.1;  % [nS]
Ca0 = 2500; % [uM]
gKCa = 8.5; % [nS]
gL = 0.23;  % [nS]
El = -21;   % [mV]
F = 9.649e05;   % [cmol^(-1)] Faraday constant
DCa = 6e-08;    % [dm^2*sec^(-1)] Ca diffusion coefficient 
Vs = 1.692e-13; % [dm^(-3)] Volume of submembrane area 
Vd = 7.356e-13; % [dm^(-3)] Volume of the deep intracellular area
Ssd = 4e-08;    % [dm^(-2)] Surface area of the submembrane and the deep intracellular area spherical boundary
dsd = 5.9e05;   % [dm] Distance between submembrane area and the deep intracellular area
Cab1max = 400;  % [uM] Total low-affinity buffer concentration
Cabhmax = 200;  % [uM] Total high-affinity buffer concentration
abl = 0.4;      % [sec^(-1)uM^(-1)] On rate constants for the binding of Ca to low-affinity buffer
bbl = 0.2;      % [sec^(-1)uM^(-1)] Off rate constants for the binding of Ca to low-affinity buffer
abh = 100;      % [sec^(-1)uM^(-1)] On rate constants for the binding of Ca to high-affinity buffer
bbh = 90;       % [sec^(-1)uM^(-1)] Off rate constants for the binding of Ca to high-affinity buffer
Jex = 9;        % [pA] Maximum Na-Ca exchanger current
Jex2 = 9.5;     % [pA] Maximum Ca-ATPase exchanger current
Camin = 0.05;   % [uM] Minimum intracellular Ca concentration for Ca extrusion

% variables
V = Y(1);
mKv = Y(2);
hKv = Y(3);
mA = Y(4);
hA = Y(5);
C1 = Y(6);
C2 = Y(7);
O1 = Y(8);
O2 = Y(9);
O3 = Y(10);
mCa = Y(11);
mKCa = Y(12);
Cas = Y(13);
Cad = Y(14);
Cabls = Y(15);
Cabhs = Y(16);
Cabld = Y(17);
Cabhd = Y(18);

% explicit functions
%%% Kv %%%
[am_Kv, bm_Kv, ah_Kv, bh_Kv] = Kv(V);
iKv = gKv*mKv^3*hKv*(V-Ek);
%%% A %%%
[am_A, bm_A, ah_A, bh_A] = A(V);
iA = gA*mA^3*hA*(V-Ek);
%%% h %%%
[ah, bh] = h(V);
ih = gh*(O1+O2+O3)*(V-Eh);
%%% Ca %%%
Eca = 12.9*log(Ca0/Cas);
[am_Ca, bm_Ca] = Ca(V);
hCa = exp(-(V-50)/11)/(exp(-(V-50)/11)+1);
iCa = gCa*mCa^4*hCa*(V-Eca);
%%% K_Ca %%%
[am_KCa, bm_KCa] = K_Ca(V);
mKc1 = Cas/(Cas+0.2);
iKCa = gKCa*mKCa^2*mKc1*(V-Ek);
%%% L %%%
iL = gL*(V-El);
%%% Iex %%%
Iex = Jex*(Cas-Camin)/(Cas-Camin+2.3)*exp(-(V+14)/70);
Iex2 = Jex2*(Cas-Camin)/(Cas-Camin+0.5);

% differential equations
dY = zeros(18, 1);
dY(1) = 1/C*-(iKv+iA+ih+iCa+iKCa+iL);
dY(2) = am_Kv*(1-mKv)-bm_Kv*mKv;
dY(3) = ah_Kv*(1-hKv)-bh_Kv*hKv;
dY(4) = am_A*(1-mA)-bm_A*mA;
dY(5) = ah_A*(1-hA)-bh_A*hA;
dY(6) = -4*ah*C1+bh*C2;
dY(7) = 4*ah*C1-(3*ah+bh)*C2+2*bh*O1;
dY(8) = 3*ah*C2-(2*ah+2*bh)*O1+3*bh*O2;
dY(9) = 2*ah*O1-(ah+3*bh)*O2+4*bh*O3;
dY(10) = ah*O2-4*bh*O3;
dY(11) = am_Ca*(1-mCa)-bm_Ca*mCa;
dY(12) = am_KCa*(1-mKCa)-bm_KCa*mKCa;
dY(13) = -iCa/(2*F*Vs)-DCa*Ssd/(Vs*dsd)*(Cas-Cad)-(Iex+Iex2)/(2*F*Vs)+...
    bbl*Cabls-abl*Cas*(Cab1max-Cabls)+bbh*Cabhs-abh*Cas*(Cabhmax-Cabhs);
dY(14) = DCa*Ssd/(Vd*dsd)*(Cas-Cad)+bbl*Cabld-abl*Cad*(Cab1max-Cabld)+...
    bbh*Cabhd-abh*Cad*(Cabhmax-Cabhd);
dY(15) = abl*Cas*(Cab1max-Cabls)+bbl*Cabls;
dY(16) = abh*Cas*(Cabhmax-Cabhs)+bbh*Cabhs;
dY(17) = abl*Cad*(Cablmax-Cabld)+bbl*Cabld;
dY(18) = abh*Cad*(Cablmax-Cabhd)+bbh*Cabhd;

end


%%%%% Kv %%%%%
function [am_Kv, bm_Kv, ah_Kv, bh_Kv] = Kv(V)

am_Kv = 400/(exp(-(V-15)/36)+1);
bm_Kv = exp(-V/13);
ah_Kv = 0.003*exp(-V/7);
bh_Kv = 80/(exp(-(V+115)/15)+1)+0.02;

end
%%%%%%%%%%%%%%%%

%%%%% A %%%%%
function [am_A, bm_A, ah_A, bh_A] = A(V)

am_A = 1200/(exp((V-50)/28)+1);
bm_A = 6*exp(-V/10);
ah_A = 0.045*exp(-V/13);
bh_A = 75/(exp(-(V+50)/15)+1);

end

%%%%% h %%%%%
function [ah, bh] = h(V)

ah = 3/(exp((V+110)/15)+1);
bh = 1.5/(exp(-(V+115)/15)+1);

end


%%%%% Ca %%%%%
function [am_Ca, bm_Ca] = Ca(V)

am_Ca = 12000*(120-V)/(exp((V-120)/25)-1);
bm_Ca = 40000/(exp((V+68)/25)+1);

end


%%%%% K_Ca %%%%%
function [am_KCa, bm_KCa] = K_Ca(V)

am_KCa = 100*(230-V)/(exp((230-V)/52)-1);
bm_KCa = 120*exp(-V/95);

end