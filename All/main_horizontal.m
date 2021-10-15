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


test_num = input('please enter the test number: ');
dt = 2*1e-05;
eps = dt;

if test_num == 1
    buffer_size = 1500000; % if set to 1000000 the result will be accurate but needs more time
    t_start = 0;
    t_end = 15;
    j_vals1 = [50, 50, 0, 0];
    j_vals2 = [50, 0, 50, 0];
    dur1 = 100000; dur2 = 102000;
    labels = {'center:50 surround:50', 'center:50 surround:0',...
        'center:0 surround:50', 'center:0 surround:0'};
end

jhvt = linspace(0,t_end,buffer_size);

figure(1);
tiledlayout(2, 2, 'TileSpacing', 'Compact');
ax1 = nexttile; ax2 = nexttile; ax3 = nexttile; ax4 = nexttile;
figure(2);
tiledlayout(4, 1, 'TileSpacing', 'Compact');
g1 = nexttile; g2 = nexttile; g3 = nexttile; g4 = nexttile;
figure(3);
rate_fig = axes;

for j = 1 : length(j_vals1)
    
if test_num == 1
    jhv = zeros(size(jhvt));
    jhv(dur1:dur2) = j_vals1(j);

    jhv2 = zeros(size(jhvt));
    jhv2(dur1:dur2) = j_vals2(j);
end

curr_t_rod = t_start;
rod = RodPhotoReceptor_HC(rod0, buffer_size, dt);
rod2 = RodPhotoReceptor_HC(rod0, buffer_size, dt);
rod3 = RodPhotoReceptor_HC(rod0, buffer_size, dt);
rod4 = RodPhotoReceptor_HC(rod0, buffer_size, dt);
rod5 = RodPhotoReceptor_HC(rod0, buffer_size, dt);
rod6 = RodPhotoReceptor_HC(rod0, buffer_size, dt);
h = Horizontal(h0, buffer_size, dt);
bip = Bipolar_complete(bip0, buffer_size, dt);
rgc = Ganglion(rgc0, buffer_size, dt);

i = 1;
glu = zeros(buffer_size, 1);
iGlu = zeros(buffer_size, 1);
hc_feed = h0(10);

tic
while abs(curr_t_rod - t_end) > eps

    input_j = jhv(round(1+buffer_size/(t_end-t_start)*curr_t_rod));
    input_j2 = jhv2(round(1+buffer_size/(t_end-t_start)*curr_t_rod));
    
    [y_rod, curr_t_rod, c_rod] = rod.solve(input_j, hc_feed);
    rod.update_time();
    [y_rod2, ~, ~] = rod2.solve(input_j2, hc_feed);
    rod2.update_time();
    [y_rod3, ~, ~] = rod3.solve(input_j2, hc_feed);
    rod3.update_time();
    [y_rod4, ~, ~] = rod4.solve(input_j2, hc_feed);
    rod4.update_time();
    [y_rod5, ~, ~] = rod5.solve(input_j2, hc_feed);
    rod5.update_time();
    [y_rod6, ~, ~] = rod6.solve(input_j2, hc_feed);
    rod6.update_time();
    
    h_in = (y_rod(1)+y_rod2(1)+y_rod3(1)+y_rod4(1)+y_rod5(1)+y_rod6(1));
    [y_h, ~, c_h] = h.solve(h_in);
    h.update_time();
    hc_feed = y_h(10);
    glu(i) = c_h(8);
    i = i + 1;
    [y_bip, ~, ~] = bip.solve(y_rod(1));
    bip.update_time();
    [y_rgc, ~, ~] = rgc.solve(y_bip(1));
    rgc.update_time();
end
end_time = toc;

t_per_step = end_time/length(jhvt);
tot_t = end_time;

t_vec = rod.get_tvec();
t_vec_end = find(t_vec == curr_t_rod);

v_rod = rod.get_V();
mCa = rod.get_mCa();
v_bip = bip.get_V();
v_rgc = rgc.get_V();
v_h = h.get_V();
gaba = h.get_gaba();

plot(ax1, t_vec(1:t_vec_end), v_rod(1:t_vec_end)); title(ax1, 'Rod Membrane Potential');
xlabel(ax1, 't [s]'); ylabel(ax1, 'V_m [mV]')
xlim(ax1, [t_vec(1), t_vec(t_vec_end)])
grid(ax1, 'on'); grid(ax1, 'minor')
hold(ax1, 'on')
plot(ax2, t_vec(1:t_vec_end), v_h(1:t_vec_end)); title(ax2, 'Horizontal Cell Membrane Potential');
xlabel(ax2, 't [s]'); ylabel(ax2, 'GABA(t) [uM]')
xlim(ax2, [t_vec(1), t_vec(t_vec_end)])
grid(ax2, 'on'); grid(ax2, 'minor')
hold(ax2, 'on')
plot(ax3, t_vec(1:t_vec_end), v_bip(1:t_vec_end)); title(ax3, 'Bipolar Cell Membrane Potential');
xlabel(ax3, 't [s]'); ylabel(ax3, 'V_{bip}(t) [mV]')
xlim(ax3, [t_vec(1), t_vec(t_vec_end)])
grid(ax3, 'on'); grid(ax3, 'minor')
hold(ax3, 'on')
plot(ax4, t_vec(1:t_vec_end), glu(1:t_vec_end)); title(ax4, 'Glutamate Concentration');
xlabel(ax4, 't [s]'); ylabel(ax4, 'Glu(t) [uM]')
xlim(ax4, [t_vec(1), t_vec(t_vec_end)])
grid(ax4, 'on'); grid(ax4, 'minor')
hold(ax4, 'on')

if j == 1
    plot(g1, t_vec(1:t_vec_end), v_rgc(1:t_vec_end))
    if test_num == 1
        title(g1, ['Ganglion Cell Membrane Potential ', char(labels(j))]);
    end
    xlabel(g1, 't [s]'); ylabel(g1, 'V_{rgc}(t) [mV]')
    xlim(g1, [t_vec(1), t_vec(t_vec_end)])
    grid(g1, 'on'); grid(g1, 'minor')
elseif j == 2
    plot(g2, t_vec(1:t_vec_end), v_rgc(1:t_vec_end))
    if test_num == 1
        title(g2, ['Ganglion Cell Membrane Potential ', char(labels(j))]);
    end
    xlabel(g2, 't [s]'); ylabel(g2, 'V_{rgc}(t) [mV]')
    xlim(g2, [t_vec(1), t_vec(t_vec_end)])
    grid(g2, 'on'); grid(g2, 'minor')
elseif j == 3
    plot(g3, t_vec(1:t_vec_end), v_rgc(1:t_vec_end))
    if test_num == 1
        title(g3, ['Ganglion Cell Membrane Potential ', char(labels(j))]);
    end
    xlabel(g3, 't [s]'); ylabel(g3, 'V_{rgc}(t) [mV]')
    xlim(g3, [t_vec(1), t_vec(t_vec_end)])
    grid(g3, 'on'); grid(g3, 'minor')
elseif j == 4
    plot(g4, t_vec(1:t_vec_end), v_rgc(1:t_vec_end))
    if test_num == 1
        title(g4, ['Ganglion Cell Membrane Potential ', char(labels(j))]);
    end
    xlabel(g4, 't [s]'); ylabel(g4, 'V_{rgc}(t) [mV]')
    xlim(g4, [t_vec(1), t_vec(t_vec_end)])
    grid(g4, 'on'); grid(g4, 'minor')
end
    

% v1 = load('rod.mat');
% v1 = v1.v_rod;
% figure; 
% plot(t_vec(1:t_vec_end), v1(1:t_vec_end));
% hold on
% plot(t_vec(1:t_vec_end), v_rod(1:t_vec_end));
% hold off
% title('Rod Membrane Potential With and Without HC'); xlabel('t [ms]'); ylabel('V [ms]')
% legend('v_{rod}', 'v_{rod+hc}')
% grid on
% grid minor

L = 10000;
out = firing_rate(v_rgc(1:t_vec_end), 0, L)/(L*dt);
plot(rate_fig, t_vec(1:t_vec_end), out); title(rate_fig, 'Firing Rate of RGC')
xlabel(rate_fig, 't [s]'); ylabel(rate_fig, 'f [Hz]')
hold(rate_fig, 'on')

end

legend(ax1, labels)
legend(ax2, labels)
legend(ax3, labels)
legend(ax4, labels)
legend(rate_fig, labels)
