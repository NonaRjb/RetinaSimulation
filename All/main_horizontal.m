clc
clear variables
close all


%%% Rod initial values (complete) %%%
rod0 = [-36.185963, 0.43018469, 0.99927789, 0.43647161, 0.64228624,...
    0.64582664, 0.29823381, 0.051645089, 0.0039748312, 0.00011472013,...
    0.096557982, 0.096558392, 80.92926, 29.067444, 80.929703, 29.067556,...
    0.0, 0.0, 0.0, 0.0, 0.300000000000037, 34.883720929940061, 2.000000000017296];

%%% Rod initial values (reduced)
% rod0 = [-36.185963, 0.99927789, 0.43647161, 0.64228624,...
%     0.29823381, 0.051645089, 0.0039748312, 0.00011472013,...
%     0.096557982, 0.096558392, 80.92926, 29.067444, 80.929703, 29.067556,...
%     0.0, 0.0, 0.0, 0.0, 0.300000000000037, 34.883720929940061, 2.000000000017296];

%%% Horizontal Cell initial values %%%
% [h_V, h_mA, h_hA, h_mKv, h_mCa, h_mNa, h_hNa, h_B, h_Ca, h_Gaba]
h0 = [-53.81, 0.03, 0.998, 0.139, 0.059, 0.026, 0.922, 0, 0.05, 0];

%%% Rod Bipolar intitial values %%%
% [bp_V, bp_mKv, bp_hKv, bp_mA, bp_hA, bp_C1, bp_C2, bp_O1, bp_O2, bp_O3,
% bp_mCa, bp_mKCa, bp_Cas, bp_Cad, bp_Cabls, bp_Cabhs, bp_Cabld, bp_Cabhd,
% rbp_syn_S]
bip0 = [-36.424516776897130, 0.824374082910590, 0.109794106890060, ...
    0.186127778073054, 0.024443585120781, 0.928238175043767, 0.054905344261992,...
    0.001217870414180, 1.200618479085905e-05, 4.438540983730620e-08, ...
    0.290203193088928, 0.475464147362579, 0.011561930331728, 0.011563608687641,...
    6.780371247710756, 1.268364765067093, 11.302574980627850, 3.805639865311822, 0];

%%% RGC initial values %%%
% [g_V, g_mNa, g_hNa, g_nCa, g_Cai, g_nK, g_mA, g_hA, bpg_S]
rgc0 = [-61.698524, 0.027940, 0.887161, 0.003019, 0.100007, 0.107809, ...
    0.070467, 0.300453, 0.0];

buffer_size = 1000000; % if set to 1000000 the result will be accurate but needs more time
t_start = 0;
t_end = 10;
dt = 2*1e-05;
eps = dt;
rod = RodPhotoReceptor_HC(rod0, buffer_size, dt);
h = Horizontal(h0, buffer_size, dt);
bip = Bipolar_complete(bip0, buffer_size, dt);
rgc = Ganglion(rgc0, buffer_size, dt);
jhvt = linspace(0,t_end,buffer_size);

%%% input sample #1
% jhv = (10)*ones(size(jhvt));

%%% input sample #2
% jhv = zeros(size(jhvt));
% jhv(500000:end) = 10;

%%% input sample #3
jhv = zeros(size(jhvt));
jhv(100000:102000) = 100;
%f0 = 0.1;
%jhv = sin(2*pi*f0*jhvt);

%iPhoto = zeros(size(jhvt))';
curr_t_rod = t_start;
%ab = zeros(10, buffer_size);
i = 1;
%maxD = [];
glu = zeros(buffer_size, 1);
iGlu = zeros(buffer_size, 1);
iCa = zeros(buffer_size, 1);
hc_feed = h0(10);

tic
while abs(curr_t_rod - t_end) > eps
%for i = 1 : buffer_size
    % input_j = interp1(jhvt, jhv, curr_t_rod);
    input_j = jhv(round(1+buffer_size/(t_end-t_start)*curr_t_rod));
    [y_rod, curr_t_rod, c_rod] = rod.solve(input_j, hc_feed);
    rod.update_time();
    [y_h, ~, c_h] = h.solve(y_rod(1));
    h.update_time();
    hc_feed = y_h(10);
    glu(i) = c_h(8);
    iGlu(i) = c_h(6);
    iCa(i) = c_rod(10);
    i = i + 1;
%     [y_bip, ~, ~] = bip.solve(y_rod(1));
%     bip.update_time();
%     [y_rgc, ~, ~] = rgc.solve(y_bip(1));
%     rgc.update_time();
    %maxD = [maxD, objdt];
    %iPhoto(i) = c_rod(end);
    %ab(:, i) = c(10:19);
    %i = i + 1;
end
end_time = toc;

t_per_step = end_time/length(jhvt);
tot_t = end_time;

t_vec = rod.get_tvec();
t_vec_end = find(t_vec == curr_t_rod);

v_rod = rod.get_V();
% v_bip = bip.get_V();
% s_bip = bip.get_S();
% v_rgc = rgc.get_V();
% s_rgc = rgc.get_S();
v_h = h.get_V();
gaba = h.get_gaba();


figure
subplot(3,2,1)
plot(jhvt, jhv); title('Input Photocurrent'); xlabel('t [s]'); ylabel('I [Rh/s]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
subplot(3,2,2)
plot(t_vec(1:t_vec_end), v_rod(1:t_vec_end)); title('Rod Membrane Potential'); xlabel('t [s]'); ylabel('V_m [mV]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
subplot(3,2,3)
plot(t_vec(1:t_vec_end), v_h(1:t_vec_end)); title('Horizontal Cell Membrane Potential'); xlabel('t [s]'); ylabel('GABA(t) [uM]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
subplot(3,2,4)
plot(t_vec(1:t_vec_end), gaba(1:t_vec_end)); title('GABA Concentration'); xlabel('t [s]'); ylabel('GABA(t) [uM]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
subplot(3,2,5)
plot(t_vec(1:t_vec_end), glu(1:t_vec_end)); title('Glutamate Concentration'); xlabel('t [s]'); ylabel('Glu(t) [uM]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
subplot(3,2,6)
plot(t_vec(1:t_vec_end), iGlu(1:t_vec_end)); title('Glutamate Current'); xlabel('t [s]'); ylabel('iGlu(t) [pA]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor

figure
subplot(2,1,1)
plot(t_vec(1:t_vec_end), v_rod(1:t_vec_end)); title('Rod Membrane Potential'); xlabel('t [s]'); ylabel('V_r [mV]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
hold on
plot(t_vec(1:t_vec_end), v_h(1:t_vec_end)); title('Horizontal Cell Membrane Potential'); xlabel('t [s]'); ylabel('V_h [mV]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
hold off
legend('V_Rod', 'V_HC')
subplot(2,1,2)
plot(t_vec(1:t_vec_end), iCa(1:t_vec_end)); title('Rod Calcium Current'); xlabel('t [s]'); ylabel('I_{Ca} [pA]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor

v1 = load('rod.mat');
v1 = v1.v_rod;
figure; 
plot(t_vec(1:t_vec_end), v1(1:t_vec_end));
hold on
plot(t_vec(1:t_vec_end), v_rod(1:t_vec_end));
hold off
legend('v_{rod}', 'v_{rodhc}')
grid on
grid minor
% subplot(3,2,3)
% plot(t_vec(1:t_vec_end), s_bip(1:t_vec_end)); title('Synaptic Value'); xlabel('t [s]'); ylabel('S(t) [mV]')
% xlim([t_vec(1), t_vec(t_vec_end)])
% grid on
% grid minor
% subplot(3,2,4)
% plot(t_vec(1:t_vec_end), v_bip(1:t_vec_end)); title('RBP Membrane Potential'); xlabel('t [s]'); ylabel('V_m [mV]')
% xlim([t_vec(1), t_vec(t_vec_end)])
% grid on
% grid minor
% subplot(3,2,5)
% plot(t_vec(1:t_vec_end), s_rgc(1:t_vec_end)); title('Synaptic Value'); xlabel('t [s]'); ylabel('S(t) [mV]')
% xlim([t_vec(1), t_vec(t_vec_end)])
% grid on
% grid minor
% subplot(3,2,6)
% plot(t_vec(1:t_vec_end), v_rgc(1:t_vec_end)); title('RGC Membrane Potential'); xlabel('t [s]'); ylabel('V_m [mV]')
% xlim([t_vec(1), t_vec(t_vec_end)])
% grid on
% grid minor
