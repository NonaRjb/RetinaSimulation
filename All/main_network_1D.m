clc 
clear variables
close all

%%% Rod initial values (complete) %%%
rod0 = [-36.185963, 0.43018469, 0.99927789, 0.43647161, 0.64228624,...
    0.64582664, 0.29823381, 0.051645089, 0.0039748312, 0.00011472013,...
    0.096557982, 0.096558392, 80.92926, 29.067444, 80.929703, 29.067556,...
    0.0, 0.0, 0.0, 0.0, 0.300000000000037, 34.883720929940061, 2.000000000017296];

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
num = 5;

rods = RodPhotoReceptor_RK.empty(num, 0);
bips = Bipolar_complete.empty(num, 0);
rgcs = Ganglion.empty(num, 0);

for i = 1 : num
    rods(i) = RodPhotoReceptor_RK(rod0, buffer_size, dt);
    bips(i) = Bipolar_complete(bip0, buffer_size, dt);
    rgcs(i) = Ganglion(rgc0, buffer_size, dt);
end

jhvt = linspace(0,t_end,buffer_size);
f0 = [0.25, 0.5, 1, 2, 4];

%%% input sample #1
% jhv = (10)*ones(size(jhvt));

%%% input sample #2
% jhv = zeros(size(jhvt));
% jhv(500000:end) = 10;

%%% input sample #3
%jhv = zeros(size(jhvt));
%jhv(100000:102000) = 100;

curr_t_rod = t_start*ones(num, 1);
pulse_dur = [20, 50, 100, 500, 1000];

tic
for i = 1 : num
    %jhv = zeros(size(jhvt));
    %jhv(100000:100000+pulse_dur(i)*100) = 100;
    jhv = sin(2*pi*f0(i)*jhvt);
    while abs(curr_t_rod(i) - t_end) > eps
        input_j = jhv(round(1+buffer_size/(t_end-t_start)*curr_t_rod(i)));
        [y_rod, curr_t_rod(i), c_rod] = rods(i).solve(input_j);
        rods(i).update_time();
        [y_bip, ~, ~] = bips(i).solve(y_rod(1));
        bips(i).update_time();
        [y_rgc, ~, ~] = rgcs(i).solve(y_bip(1));
        rgcs(i).update_time();
    end
end
end_time = toc;

t_per_step = end_time/length(jhvt);
tot_t = end_time;

for i = 1 : num
    t_vec = rods(i).get_tvec();
    t_vec_end = find(t_vec == curr_t_rod(i));

    v_rod = rods(i).get_V();
    v_bip = bips(i).get_V();
    s_bip = bips(i).get_S();
    v_rgc = rgcs(i).get_V();
    s_rgc = rgcs(i).get_S();
    
    jhv = zeros(size(jhvt));
    jhv(100000:100000+pulse_dur(i)*100) = 100;

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
    plot(t_vec(1:t_vec_end), s_bip(1:t_vec_end)); title('Synaptic Value'); xlabel('t [s]'); ylabel('S(t) [mV]')
    xlim([t_vec(1), t_vec(t_vec_end)])
    grid on
    grid minor
    subplot(3,2,4)
    plot(t_vec(1:t_vec_end), v_bip(1:t_vec_end)); title('RBP Membrane Potential'); xlabel('t [s]'); ylabel('V_m [mV]')
    xlim([t_vec(1), t_vec(t_vec_end)])
    grid on
    grid minor
    subplot(3,2,5)
    plot(t_vec(1:t_vec_end), s_rgc(1:t_vec_end)); title('Synaptic Value'); xlabel('t [s]'); ylabel('S(t) [mV]')
    xlim([t_vec(1), t_vec(t_vec_end)])
    grid on
    grid minor
    subplot(3,2,6)
    plot(t_vec(1:t_vec_end), v_rgc(1:t_vec_end)); title('RGC Membrane Potential'); xlabel('t [s]'); ylabel('V_m [mV]')
    xlim([t_vec(1), t_vec(t_vec_end)])
    grid on
    grid minor
    
    spikes = double(v_rgc(1:t_vec_end) > 0);
    f = gausswin(10000);
    out = filter(f,1,spikes);
    figure; plot(t_vec(1:t_vec_end), out);
    
end
