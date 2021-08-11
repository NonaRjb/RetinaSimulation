clc
close all
clear variables


tspan = [0, 10];
tstep = 1e-06;
tspan_fix = tspan(1):tstep:tspan(2);

% [rod_V, rod_mKv, rod_hKv, rod_mCa, rod_mKCa, rod_C1, rod_C2, rod_O1,
% rod_O2, rod_O3, rod_Cas, rod_Caf, rod_Cabls, rod_Cabhs, rod_Cablf,
% rod_Cabhf, rod_Rh, rod_Rhi, rod_Tr, rod_PDE, rod_Ca2, rod_Cab, rod_cGMP,
% rbp_syn_S,
% bp_V, bp_mKv, bp_hKv, bp_mA, bp_hA, bp_C1, bp_C2, bp_O1, bp_O2, bp_O3,
% bp_mCa, bp_mKCa, bp_Cas, bp_Cad, bp_Cabls, bp_Cabhs, bp_Cabld, bp_Cabhd]
Y0 = [-36.185963, 0.43018469, 0.99927789, 0.43647161, 0.64228624,...
    0.64582664, 0.29823381, 0.051645089, 0.0039748312, 0.00011472013,...
    0.096557982, 0.096558392, 80.92926, 29.067444, 80.929703, 29.067556,...
    0.0, 0.0, 0.0, 0.0, 0.300000000000037, 34.883720929940061, 2.000000000017296,...
    0,...
    -36.424516776897130, 0.824374082910590, 0.109794106890060, ...
    0.186127778073054, 0.024443585120781, 0.928238175043767, 0.054905344261992,...
    0.001217870414180, 1.200618479085905e-05, 4.438540983730620e-08, ...
    0.290203193088928, 0.475464147362579, 0.011561930331728, 0.011563608687641,...
    6.780371247710756, 1.268364765067093, 11.302574980627850, 3.805639865311822];
jhvt = linspace(0,10,1000000);

%%% input sample #1
%jhv = (10)*ones(size(jhvt));

%%% input sample #2
%jhv = zeros(size(jhvt));
%jhv(10000:end) = 10;

%%% input sample #3
jhv = zeros(size(jhvt));
jhv(100000:102000) = 100;

tic
[t, y] = ode15s(@(t, Y) odefuncs_retina(t, Y, jhvt, jhv), tspan, Y0);
end_time = toc;

t_per_step = end_time/length(jhvt);
tot_t = end_time;

jmax = 5040;
j = jmax*(y(:,23)).^3./(y(:,23).^3+10^3);
iPhoto = -j.*(1.0-exp((y(:,1)-8.5)/17.0));
figure
subplot(2, 1, 1)
plot(jhvt, jhv); title('Input Stimulus'); xlabel('t [s]'); ylabel('J_{hv}(t) [Rh^*.s^{-1}] ')
grid on
grid minor
subplot(2, 1, 2)
plot(t, iPhoto); title('Photocurrent'); xlabel('t [s]'); ylabel('I_{Photo}(t) [pA]')
grid on
grid minor
% subplot(3, 1, 3)
% plot(t,y(:, 1)); title('Membrane Potential'); xlabel('t [s]'); ylabel('V_m [mV]')
% grid on
% grid minor

figure
subplot(3,1,1)
plot(t,y(:, 1)); title('Rod Membrane Potential'); xlabel('t [s]'); ylabel('V_m [mV]')
grid on
grid minor
subplot(3,1,2)
plot(t, y(:, 24)); title('Synaptic Value'); xlabel('t [s]'); ylabel('S(t) [mV]')
grid on
grid minor
subplot(3,1,3)
plot(t, y(:, 25)); title('RBP Membrane Potential'); xlabel('t [s]'); ylabel('V_m [mV]')
grid on
grid minor



% % this is the "finite difference" derivative. Note it is  one element shorter than y and x
% yd = diff(y(:, 1))./diff(t);
% % this is to assign yd an abscissa midway between two subsequent x
% xd = (t(2:end)+t(1:(end-1)))/2;
% % this should be a rough plot of your derivative
% figure
% plot(xd,yd)