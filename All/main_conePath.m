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
rgc0 = [-61.698524,...  % g_V
    0.027940,...        % g_mNa
    0.887161,...        % g_hNa
    0.003019,...        % g_nCa
    0.100007,...        % g_Cai
    0.107809,...        % g_nK
    0.070467,...        % g_mA
    0.300453,...        % g_hA
    0.000000];          % bpg_syn_S

% buffer_size = 1000000;
% t_start = 0;
% t_end = 10;
test_num = input('please enter the test number: ');
dt = 2*1e-05;
eps = dt;
method = 'euler';
% cone = ConePhotoReceptor(cone0, buffer_size, dt, method);
% bip = Bipolar_complete(bip0, buffer_size, dt, method);
% rgc = Ganglion(rgc0, buffer_size, dt, method);
% offset = -1*cone0(1);
% bip.set_Vth(cone0(1)+offset);
% rgc.set_Vth(bip0(1));

% jhvt = linspace(0,t_end,buffer_size);

%%% input sample #1
%jhv = (10)*ones(size(jhvt));

%%% input sample #2
% jhv = zeros(size(jhvt));
% jhv(500000:end) = 10;

%%% input sample #3
% jhv = zeros(size(jhvt));
% jhv(100000:102000) = 100;

% curr_t_cone = t_start;
%ab = zeros(10, buffer_size);
%i = 1;

if test_num == 1
    buffer_size = 1000000;
    t_start = 0;
    t_end = 10;
    j_vals = [0, 100, 5000, 50000, 890000];
    dur1 = 100000; dur2 = 102000;
    labels = {'0', '100', '5000', '50000', '890000'};
elseif test_num == 2
    buffer_size = 1000000;
    t_start = 0;
    t_end = 10;
    j_vals = [1000, 2000, 10000, 50000, 500000];
    labels = {'10^{ms}', '20^{ms}', '0.1^s', '0.5^s', '5^s'};
end

jhvt = linspace(0,t_end,buffer_size);

figure(1)
tiledlayout(3, 1, 'TileSpacing', 'Compact')
ax1 = nexttile; ax2 = nexttile; ax3 = nexttile;
figure(2)
tiledlayout(5, 1, 'TileSpacing', 'Compact')
g1 = nexttile; g2 = nexttile; g3 = nexttile; g4 = nexttile; g5 = nexttile;
figure(3)
rate_fig = axes;

for j = 1 : length(j_vals)
if test_num == 1
    jhv = zeros(size(jhvt));
    jhv(dur1:dur2) = j_vals(j);
elseif test_num == 2
    jhv = zeros(size(jhvt));
    jhv(100000:100000+j_vals(j)) = 50000;
end
curr_t_cone = t_start;
cone = ConePhotoReceptor(cone0, buffer_size, dt, method);
bip = Bipolar_complete(bip0, buffer_size, dt, method);
rgc = Ganglion(rgc0, buffer_size, dt, method);
offset = -1*cone0(1);
bip.set_Vth(cone0(1)+offset);
rgc.set_Vth(bip0(1));
tic
while abs(curr_t_cone - t_end) > eps
    input_j = jhv(round(1+buffer_size/(t_end-t_start)*curr_t_cone));
    [y_cone, curr_t_cone, c_cone] = cone.solve(input_j);
    cone.update_time();
    [y_bip, ~, ~] = bip.solve(y_cone(1)+offset);
    bip.update_time();
    [y_rgc, ~, ~] = rgc.solve(y_bip(1));
    rgc.update_time();
    %ab(:, i) = c(10:19);
    %i = i + 1;
end
end_time = toc;

t_per_step = end_time/length(jhvt);
tot_t = end_time;

t_vec = cone.get_tvec();
t_vec_end = find(t_vec == curr_t_cone);

v_cone = cone.get_V();
v_cone = v_cone+offset;
v_bip = bip.get_V();
s_bip = bip.get_S();
v_rgc = rgc.get_V();
s_rgc = rgc.get_S();

plot(ax1, jhvt, jhv); title(ax1, 'Input Stimulus');
xlabel(ax1, 't [s]'); ylabel(ax1, 'I_{in} [Photoisomerizations/s]')
xlim(ax1, [t_vec(1), t_vec(t_vec_end)])
grid(ax1, 'on'); grid(ax1, 'minor')
hold(ax1, 'on')

plot(ax2, t_vec(1:t_vec_end), v_cone(1:t_vec_end)); title(ax2, 'Cone Membrane Potential');
xlabel(ax2, 't [s]'); ylabel(ax2, 'V_c [mV]')
xlim(ax2, [t_vec(1), t_vec(t_vec_end)])
grid(ax2, 'on'); grid(ax2, 'minor')
hold(ax2, 'on')

plot(ax3, t_vec(1:t_vec_end), v_bip(1:t_vec_end)); title(ax3, 'Cone Bipolar Membrane Potential');
xlabel(ax3, 't [s]'); ylabel(ax3, 'V_{cbp} [mV]')
xlim(ax3, [t_vec(1), t_vec(t_vec_end)])
grid(ax3, 'on'); grid(ax3, 'minor')
hold(ax3, 'on')

if j == 1
    plot(g1, t_vec(1:t_vec_end), v_rgc(1:t_vec_end))
    if test_num == 1
        title(g1, ['Ganglion Membrane Potential ', '(', num2str(j_vals(j)), ' [photoisomerization/s])']);
    elseif test_num == 2
        title(g1, ['Ganglion Membrane Potential ', '(duration = ', num2str(j_vals(j)/100000), ' [s])']);
    end
    xlabel(g1, 't [s]'); ylabel(g1, 'V_{rgc} [mV]')
    xlim(g1, [t_vec(1), t_vec(t_vec_end)])
    grid(g1, 'on')
    grid(g1, 'minor')
elseif j == 2
    plot(g2, t_vec(1:t_vec_end), v_rgc(1:t_vec_end))
    if test_num == 1
        title(g2, ['Ganglion Membrane Potential ', '(', num2str(j_vals(j)), ' [photoisomerization/s])']);
    elseif test_num == 2
        title(g2, ['Ganglion Membrane Potential ', '(duration = ', num2str(j_vals(j)/100000), ' [s])']);
    end
    xlabel(g2, 't [s]'); ylabel(g2, 'V_{rgc} [mV]')
    xlim(g2, [t_vec(1), t_vec(t_vec_end)])
    grid(g2, 'on')
    grid(g2, 'minor')
elseif j == 3
    plot(g3, t_vec(1:t_vec_end), v_rgc(1:t_vec_end))
    if test_num == 1
        title(g3, ['Ganglion Membrane Potential ', '(', num2str(j_vals(j)), ' [photoisomerization/s])']);
    elseif test_num == 2
        title(g3, ['Ganglion Membrane Potential ', '(duration = ', num2str(j_vals(j)/100000), ' [s])']);
    end
    xlabel(g3, 't [s]'); ylabel(g3, 'V_{rgc} [mV]')
    xlim(g3, [t_vec(1), t_vec(t_vec_end)])
    grid(g3, 'on')
    grid(g3, 'minor')
elseif j == 4
    plot(g4, t_vec(1:t_vec_end), v_rgc(1:t_vec_end))
    if test_num == 1
        title(g4, ['Ganglion Membrane Potential ', '(', num2str(j_vals(j)), ' [photoisomerization/s])']);
    elseif test_num == 2
        title(g4, ['Ganglion Membrane Potential ', '(duration = ', num2str(j_vals(j)/100000), ' [s])']);
    end
    xlabel(g4, 't [s]'); ylabel(g4, 'V_{rgc} [mV]')
    xlim(g4, [t_vec(1), t_vec(t_vec_end)])
    grid(g4, 'on')
    grid(g4, 'minor')
elseif j == 5
    plot(g5, t_vec(1:t_vec_end), v_rgc(1:t_vec_end))
    if test_num == 1
        title(g5, ['Ganglion Membrane Potential ', '(', num2str(j_vals(j)), ' [photoisomerization/s])']);
    elseif test_num == 2
        title(g5, ['Ganglion Membrane Potential ', '(duration = ', num2str(j_vals(j)/100000), ' [s])']);
    end
    xlabel(g5, 't [s]'); ylabel(g5, 'V_{rgc} [mV]')
    xlim(g5, [t_vec(1), t_vec(t_vec_end)])
    grid(g5, 'on')
    grid(g5, 'minor')
end

L = 10000;
out = firing_rate(v_rgc(1:t_vec_end), 0, L)/(L*dt);
plot(rate_fig, t_vec(1:t_vec_end), out); title(rate_fig, 'Firing Rate of RGC')
xlabel(rate_fig, 't [s]'); ylabel(rate_fig, 'f [Hz]')
hold(rate_fig, 'on')

end

legend(ax1, labels)
legend(ax2, labels)
legend(ax3, labels)
legend(rate_fig, labels)
% figure 
% plot(t_vec(1:t_vec_end), 1./(ab(1,1:t_vec_end)+ab(2,1:t_vec_end)))
% hold on 
% plot(t_vec(1:t_vec_end), 1./(ab(3,1:t_vec_end)+ab(4,1:t_vec_end)))
% hold on 
% plot(t_vec(1:t_vec_end), 1./(ab(5,1:t_vec_end)+ab(6,1:t_vec_end)))
% hold on
% plot(t_vec(1:t_vec_end), 1./(ab(7,1:t_vec_end)+ab(8,1:t_vec_end)))
% hold on
% plot(t_vec(1:t_vec_end), 1./(ab(9,1:t_vec_end)+ab(10,1:t_vec_end)))
% hold off
% legend('\tau_{mKv}', '\tau_{hKv}', '\tau_{mCa}', '\tau_{mKCa}', '\tau_{h}')

