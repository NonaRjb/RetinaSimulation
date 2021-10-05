clc
clear variables
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for this test I_A is not included in calculation of dV/dt. Also I_syn has been
% replaced by Vpre to play the role of external stimualation and has a positive-sign
% contribution in calculating dV/dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% initial values for the bipolar cell
% [V, mKv, hKv, mA, hA, C1, C2, O1, O2, O3, mCa, mKCa, Cas, Cad, Cabls, 
%  Cabhs, Cabld, Cabhd]
bip0 = [-36.424516776897130, 0.824374082910590, 0.109794106890060, ...
    0.186127778073054, 0.024443585120781, 0.928238175043767, 0.054905344261992,...
    0.001217870414180, 1.200618479085905e-05, 4.438540983730620e-08, ...
    0.290203193088928, 0.475464147362579, 0.011561930331728, 0.011563608687641,...
    6.780371247710756, 1.268364765067093, 11.302574980627850, 3.805639865311822, 0];

test_num = input('please enter the test number: ');
dt = 2*1e-05;
eps = dt;
method = 'euler';

%curr_t_bip = t_start;
%ab = zeros(10, buffer_size);
%i = 1;
if test_num == 1
    buffer_size = 500000;
    t_start = 0;
    t_end = 5;
    j_vals = [10, 20, 30];
    jhvt = linspace(0,t_end,buffer_size);
    labels = {'10 pA', '20 pA', '30 pA'};
    dur1 = 50000; dur2 = 250000;
elseif test_num == 2
    buffer_size = 500000;
    t_start = 0;
    t_end = 5;
    j_vals = [-10, -20, -30];
    jhvt = linspace(0,t_end,buffer_size);
    labels = {'-10 pA', '-20 pA', '-30 pA'};
    dur1 = 50000; dur2 = 250000;
else
    error('The commanded test does not exist!')
end

tiledlayout(1, 6, 'TileSpacing', 'Compact')
ax1 = nexttile; ax2 = nexttile; ax3 = nexttile; ax4 = nexttile;
ax5 = nexttile; ax6 = nexttile;
for j = 1 : length(j_vals)
bip = Bipolar_complete(bip0, buffer_size, dt, method);
curr_t_bip = t_start;
jhv = zeros(size(jhvt));
jhv(dur1:dur2) = j_vals(j);
i = 1;
iKv = zeros(buffer_size, 1);
ih = zeros(buffer_size, 1);
iCa = zeros(buffer_size, 1);
iKCa = zeros(buffer_size, 1);
iL = zeros(buffer_size, 1);
%tic
while curr_t_bip < t_end && jhvt(end)>=curr_t_bip
    input_j = jhv(round(1+buffer_size/(t_end-t_start)*curr_t_bip));
    [y_bip, curr_t_bip, c_bip] = bip.solve(input_j);
    bip.update_time();
    iKv(i) = c_bip(5);
    ih(i) = c_bip(13);
    iCa(i) = c_bip(18);
    iKCa(i) = c_bip(22);
    iL(i) = c_bip(23);
    %ab(:, i) = c(10:19);
    i = i + 1;
end
%end_time = toc;

%t_per_step = end_time/length(jhvt);
%tot_t = end_time;

t_vec = bip.get_tvec();
t_vec_end = find(t_vec == curr_t_bip);

v = bip.get_V();

plot(ax1, t_vec(1:t_vec_end), v(1:t_vec_end)); title(ax1, 'Membrane Potential'); xlabel(ax1, 't [s]'); ylabel(ax1, 'V_m [mV]')
xlim(ax1, [t_vec(1), t_vec(t_vec_end)]); ylim(ax1, [-140, 20])
grid(ax1, 'on'); grid(ax1, 'minor');
hold(ax1, 'on')
plot(ax2, t_vec(1:t_vec_end), iKv(1:t_vec_end)); title(ax2, 'I_{Kv}'); xlabel(ax2, 't [s]'); ylabel(ax2, 'I_{Kv} [pA]')
xlim(ax2, [t_vec(1), t_vec(t_vec_end)]); ylim(ax2, [-25, 15])
grid(ax2, 'on'); grid(ax2, 'minor');
hold(ax2, 'on')
plot(ax3, t_vec(1:t_vec_end), iCa(1:t_vec_end)); title(ax3, 'I_{Ca}'); xlabel(ax3, 't [s]'); ylabel(ax3, 'I_{Ca} [pA]')
xlim(ax3, [t_vec(1), t_vec(t_vec_end)]); ylim(ax3, [-25, 15])
grid(ax3, 'on'); grid(ax3, 'minor'); 
hold(ax3, 'on')
plot(ax4, t_vec(1:t_vec_end), iKCa(1:t_vec_end)); title(ax4, 'I_{K(Ca)}'); xlabel(ax4, 't [s]'); ylabel(ax4, 'I_{K(Ca)} [pA]')
xlim(ax4, [t_vec(1), t_vec(t_vec_end)]); ylim(ax4, [-25, 15])
grid(ax4, 'on'); grid(ax4, 'minor'); 
hold(ax4, 'on')
plot(ax5, t_vec(1:t_vec_end), ih(1:t_vec_end)); title(ax5, 'I_{h}'); xlabel(ax5, 't [s]'); ylabel(ax5, 'I_{h} [pA]')
xlim(ax5, [t_vec(1), t_vec(t_vec_end)]); ylim(ax5, [-25, 15])
grid(ax5, 'on'); grid(ax5, 'minor'); 
hold(ax5, 'on')
plot(ax6, t_vec(1:t_vec_end), iL(1:t_vec_end)); title(ax6, 'I_{L}'); xlabel(ax6, 't [s]'); ylabel(ax6, 'I_{L} [pA]')
xlim(ax6, [t_vec(1), t_vec(t_vec_end)]); ylim(ax6, [-25, 15])
grid(ax6, 'on'); grid(ax6, 'minor'); 
hold(ax6, 'on')
%hold([ax1, ax2, ax3, ax4, ax5, ax6], 'on')
end