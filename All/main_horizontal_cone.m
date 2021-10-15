clc
clear variables
close all


%%% Cone initial values %%%
cone0 = [-42.2321233629648844498660764656960964202880859375,...         % c_V
    0.004863366238889023189517768486211934941820800304412841796875,...  % c_Gf
    0.000000000000000000000000000000000000000000000000000000000000,...  % c_y1
    0.000000000000000000000000000000000000000000000000000000000000,...  % c_y2
    0.000000000000000000000000000000000000000000000000000000000000,...  % c_y3
    0.000000000000000000000000000000000000000000000000000000000000,...  % c_y4
    0.000000000000000000000000000000000000000000000000000000000000,...  % c_y5
    0.000000000000000000000000000000000000000000000000000000000000,...  % c_z1
    0.000000000000000000000000000000000000000000000000000000000000,...  % c_z2
    0.000000000000000000000000000000000000000000000000000000000000];    % c_z3

%%% Horizontal Cell initial values %%%
% [h_V, h_mA, h_hA, h_mKv, h_mCa, h_mNa, h_hNa, h_B, h_Ca, h_Gaba]
h0 = [-53.81, 0.03, 0.998, 0.139, 0.059, 0.026, 0.922, 0, 0.05, 0];

%%% Cone Bipolar intitial values %%%
% [bp_V, bp_mKv, bp_hKv, bp_mA, bp_hA, bp_C1, bp_C2, bp_O1, bp_O2, bp_O3,
% bp_mCa, bp_mKCa, bp_Cas, bp_Cad, bp_Cabls, bp_Cabhs, bp_Cabld, bp_Cabhd,
% rbp_syn_S]
bip0 = [-36.424516776897130,...
    0.824374082910590,...
    0.109794106890060,...
    0.186127778073054,...
    0.024443585120781,...
    0.928238175043767,...
    0.054905344261992,...
    0.001217870414180,...
    1.200618479085905e-05,...
    4.438540983730620e-08,...
    0.290203193088928,...
    0.475464147362579,...
    0.011561930331728,...
    0.011563608687641,...
    6.780371247710756,...
    1.268364765067093,...
    11.302574980627850,...
    3.805639865311822,...
    0.000000000000000];

%%% RGC initial values %%%
% [g_V, g_mNa, g_hNa, g_nCa, g_Cai, g_nK, g_mA, g_hA, bpg_S]
rgc0 = [-61.698524, 0.027940, 0.887161, 0.003019, 0.100007, 0.107809, ...
    0.070467, 0.300453, 0.0];

buffer_size = 1000000; % if set to 1000000 the result will be accurate but needs more time
t_start = 0;
t_end = 10;
dt = 2*1e-05;
eps = dt;
cone = ConePhotoReceptor_HC(cone0, buffer_size, dt);
h = Horizontal(h0, buffer_size, dt);
bip = Bipolar_complete(bip0, buffer_size, dt);
rgc = Ganglion(rgc0, buffer_size, dt);
jhvt = linspace(0,t_end,buffer_size);

%%% input sample #3
jhv = zeros(size(jhvt));
jhv(100000:102000) = 50000;
%f0 = 0.1;
%jhv = sin(2*pi*f0*jhvt);

%iPhoto = zeros(size(jhvt))';
curr_t = t_start;
%ab = zeros(10, buffer_size);
i = 1;
%maxD = [];
glu = zeros(buffer_size, 1);
iGlu = zeros(buffer_size, 1);
hc_feed = h0(10);

tic
while abs(curr_t - t_end) > eps
%for i = 1 : buffer_size
    % input_j = interp1(jhvt, jhv, curr_t_rod);
    input_j = jhv(round(1+buffer_size/(t_end-t_start)*curr_t));
    [y_cone, curr_t, c_cone] = cone.solve(input_j);
    cone.update_time();
    [y_h, ~, c_h] = h.solve(y_cone(1));
    h.update_time();
    hc_feed = y_h(10);
    glu(i) = c_h(8);
    iGlu(i) = c_h(6);
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

t_vec = cone.get_tvec();
t_vec_end = find(t_vec == curr_t);

v_cone = cone.get_V();
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
plot(t_vec(1:t_vec_end), v_cone(1:t_vec_end)); title('Cone Membrane Potential'); xlabel('t [s]'); ylabel('V_m [mV]')
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
plot(t_vec(1:t_vec_end), v_cone(1:t_vec_end)-cone0(1)); title('Cone Membrane Potential'); xlabel('t [s]'); ylabel('V_c [mV]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
hold on
plot(t_vec(1:t_vec_end), v_h(1:t_vec_end)-h0(1)); title('Horizontal Cell Membrane Potential'); xlabel('t [s]'); ylabel('V_h [mV]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
hold off
legend('V_Rod', 'V_HC')

% v1 = load('rod.mat');
% v1 = v1.v_rod;
% figure; 
% plot(t_vec(1:t_vec_end), v1(1:t_vec_end));
% hold on
% plot(t_vec(1:t_vec_end), v_cone(1:t_vec_end));
% hold off
% legend('v_{rod}', 'v_{rodhc}')
% grid on
% grid minor

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
