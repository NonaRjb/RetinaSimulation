function [Cas, Caf, Cals, Cahs, Calf, Cahf, Iex, Iex2] = ...
    Ca(V, Cas0, Caf0, Cals0, Cahs0, Cahf0, iex0, iex20)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% constants
F = 9.648 * 10 ^ 4; % Faraday constant
Cae = 0.05; % [uM] minimum intracellular Ca concentration
BH = 500; % [uM] maximum concentration of high-affinity buffer
BL = 300; % [uM] maximum concentration of low-affinity buffer
Lb1 = 0.4; Lb2 = 0.2; % [1/(sec*uM)] on and off rate constants for the binding of Ca to low affinity buffer
Hb1 = 100; Hb2 = 90; % [1/(sec*uM)] on and off rate constants for the binding of Ca to high affinity buffer
V1 = 3.664 * 10 ^ -13; % [dm^3] volume of submembrane area
V2 = 5.236 * 10 ^ -13; % [dm^3] volume of deep intracellular area
S1 = 3.14 * 10 ^ -6; % [cm^2] submembrane area
DCa = 6 * 10 ^ -8; % [dm^2/s] Ca diffusion coefficient
delta = 3 * 10 ^ -4; % [cm] diffusion distance
Imax = 20; % [pA] maximum Na+ - Ca2+ exchanger current
Imax2 = 30; % [pA] maximum Ca-ATPase pump current

Cas = iCa/(2*F*V1) - (iex0+iex20)/(2*F*V1) - DCa*S1/(V1*delta)*(Cas0-Caf0) + ...
    Lb2*Cals0 - Lb1*Cas0*(BL-Cals0) - Hb1*Cas0*(BH-Cahs0) + Hb2*Cahs0 + Cas0;
Caf = DCa*S1/(V2*delta)*(Cas0-Caf0) + Lb2*Calf0 - Lb1*Caf0*(BL-Calf0) - ...
    Hb1*Caf0*(BH-Cahf0) + Hb2*Cahf0 + Caf0;
Cals = Lb1*Cas0*(BL-Cals0) - Lb2*(Cals0) + Cals0;
Cahs = Hb1*Cas0*(BH-Cahs0) - Hb2*Cahs0 + Cahs0;
Calf = Lb1*Caf0*(BL-Calf0) - Lb2*Calf0 + Calf0;
Cahf = Hb1*Caf0*(BH-Cahf0) - Lb2*Cahf0 + Cahf0;

Iex = Imax*exp(-(V+14)/70)*(Cas0-Cae)/(Cas0-Cae+2.3);
Iex2 = Imax2*(Cas0-Cae)/(Cas0-Cae+0.5);

end

