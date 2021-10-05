clc
clear variables
close all

%% Initial Values

%%% Rod initial values (complete) %%%
% [r_V, r_mKv, r_hKv, r_mCa, r_mKCa, r_C1, r_C2, r_O1, r_O2, r_O3, r_Cas,
% r_Caf, r_Cabls, r_Cabhs, r_Cablf, r_Cabhf, r_Rh, r_Rhi, r_Tr, r_PDE,
% r_Ca2, r_Cab, r_cGMP]
rod0 = [-35.15299217465973700000000000000000000000000000000000,... % r_V
    0.44095899101175912000000000000000000000000000000000,...       % r_mKv
    0.99912312462833097000000000000000000000000000000000,...       % r_hKv
    0.46586591802631688000000000000000000000000000000000,...       % r_mCa
    0.65322043938092211000000000000000000000000000000000,...       % r_mKCa
    0.64132528330066263000000000000000000000000000000000,...       % r_C1
    0.26412772593690498000000000000000000000000000000000,...       % r_C2
    0.04079255338037327700000000000000000000000000000000,...       % r_O1
    0.00280004667419985350000000000000000000000000000000,...       % r_O2
    0.00007207438056699847500000000000000000000000000000,...       % r_O3
    0.12381319321649724000000000000000000000000000000000,...       % r_Cas
    0.12381319648326869000000000000000000000000000000000,...       % r_Caf
    99.23899869893260500000000000000000000000000000000000,...      % r_Cabls
    36.28001495889955700000000000000000000000000000000000,...      % r_Cabhs
    99.23900194523179200000000000000000000000000000000000,...      % r_Cablf
    36.28001580109989300000000000000000000000000000000000,...      % r_Cabhf
    0.00000000000000000000000000000000000000000000000000,...       % r_Rh
    0.00000000000000000000000000000000000000000000000000,...       % r_Rhi
    0.00000000000000000000000000000000000000000000000000,...       % r_Tr
    0.00000000000000000000000000000000000000000000000000,...       % r_PDE
    0.30000000000000143000000000000000000000000000000000,...       % r_Ca2
    34.88372093023302500000000000000000000000000000000000,...      % r_Cab
    1.99999999999998670000000000000000000000000000000000];         % r_cGMP

%%% Cone initial values %%%
% [c_V, c_Gf, c_y1, c_y2, c_y3, c_y4, c_y5, c_z1, c_z2, c_z3] 
cone0 = [-45.76654665646767500000000000000000000000000000000000,...
    0.01137478132545304700000000000000000000000000000000,...
    0.00000000000000000000000000000000000000000000000000,...
    0.00000000000000000000000000000000000000000000000000,...
    0.00000000000000000000000000000000000000000000000000,...
    0.00000000000000000000000000000000000000000000000000,...
    0.00000000000000000000000000000000000000000000000000,...
    0.00000000000000000000000000000000000000000000000000,...
    0.00000000000000000000000000000000000000000000000000,...
    0.00000000000000000000000000000000000000000000000000];

%%% Rod Bipolar intitial values %%%
rbip0 = [-36.424516776897130,... % bp_V
    0.824374082910590,...        % bp_mKv
    0.109794106890060,...        % bp_hKv
    0.186127778073054,...        % bp_mA
    0.024443585120781,...        % bp_hA
    0.928238175043767,...        % bp_C1
    0.054905344261992,...        % bp_C2
    0.001217870414180,...        % bp_O1
    1.200618479085905e-05,...    % bp_O2
    4.438540983730620e-08,...    % bp_O3
    0.290203193088928,...        % bp_mCa
    0.475464147362579,...        % bp_mKCa
    0.011561930331728,...        % bp_Cas
    0.011563608687641,...        % bp_Cad
    6.780371247710756,...        % bp_Cabls
    1.268364765067093,...        % bp_Cabhs
    11.302574980627850,...       % bp_Cabld 
    3.805639865311822,...        % bp_Cabhd
    0.000000000000000];          % rbp_syn_S

%%% Cone Bipolar Initial Values %%%
cbip0 = [-40.85865173873610700638892012648284435272216796875000,... % bp_V
    0.75119613755960557011093214896391145884990692138672,...        % bp_mKv
    0.12833933300174665825821307407750282436609268188477,...        % bp_hKv
    0.11197954384331872124125339951206115074455738067627,...        % bp_mA
    0.04078196799520746040901286733060260303318500518799,...        % bp_hA
    0.91442826153121192778172598991659469902515411376953,...        % bp_C1
    0.06802394011168773346742710828038980253040790557861,...        % bp_C2
    0.00189760228726293194399799002525242030969820916653,...        % bp_O1
    0.00002352697486413260014254458796223445915529737249,...        % bp_O2
    0.00000010938512049644076438087637219778724961827265,...        % bp_O3
    0.23512534165859483303862020875385496765375137329102,...        % bp_mCa
    0.44667148688898017372395088386838324368000030517578,...        % bp_mKCa
    0.01069015667321450424598339168369420804083347320557,...        % bp_Cas
    0.01146115054172941634391680310045558144338428974152,...        % bp_Cad
    6.69114537743995629881510467384941875934600830078125,...        % bp_Cabls
    1.17390152355585408727733920386526733636856079101562,...        % bp_Cabhs
    11.29221858507642828328698669793084263801574707031250,...       % bp_Cabld 
    3.77261483646810091840961831621825695037841796875000,...        % bp_Cabhd
    0.00000000000000000000000000000000000000000000000000];          % cbp_syn_S

%%% Amacrine Cell initial values %%%
a0 = [-50.34951375105524107311794068664312362670898437500000,... % A_V
    0.24340676136333766477193307764537166804075241088867,...     % A_m
    0.15952849172211810979682411471003433689475059509277,...     % A_h
    0.54575595513993369056748861112282611429691314697266,...     % A_n
    0];                                                          % S

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

%% 

test_num = input('please enter the test number: ');
dt = 2*1e-05;
eps = dt;
method = 'euler';
gap_junc = 1;

if test_num == 1
    buffer_size = 1500000;
    t_start = 0;
    t_end = 15;
    rod_vals = [0, 10, 50, 100, 500];
    cone_vals = [0, 10, 50, 100, 500];
    dur1 = 100000; dur2 = 102000;
    labels = {'0', '10', '50', '100', '500'};
elseif test_num == 2
    buffer_size = 1500000;
    t_start = 0;
    t_end = 15;
    rod_vals = [1000, 2000, 10000, 50000, 500000];
    labels = {'10^{ms}', '20^{ms}', '0.1^s', '0.5^s', '5^s'};
end
jhvt = linspace(0,t_end,buffer_size);

figure(1)
tiledlayout(4, 2, 'TileSpacing', 'Compact')
ax1 = nexttile; ax2 = nexttile; ax3 = nexttile; ax4 = nexttile; ax5 = nexttile;
ax6 = nexttile; ax7 = nexttile;
figure(2)
tiledlayout(5, 1, 'TileSpacing', 'Compact')
g1 = nexttile; g2 = nexttile; g3 = nexttile; g4 = nexttile; g5 = nexttile;
figure(3)
rate_fig = axes;

for j = 1 : length(rod_vals)
if test_num == 1
    in_cone = zeros(size(jhvt));
    in_rod = zeros(size(jhvt));
    in_cone(dur1:dur2) = cone_vals(j);
    in_rod(dur1:dur2) = rod_vals(j);
elseif test_num == 2
    in_cone = zeros(size(jhvt));
    in_rod = zeros(size(jhvt));
    in_cone(100000:100000+rod_vals(j)) = 100;
    in_rod(100000:100000+rod_vals(j)) = 100;
end

curr_t = t_start;

cone = ConePhotoReceptor(cone0, buffer_size, dt, method, gap_junc);
rod = RodPhotoReceptor_RK(rod0, buffer_size, dt, method, gap_junc);
rod_bip = Bipolar_complete(rbip0, buffer_size, dt, method);
cone_bip = Bipolar_complete(cbip0, buffer_size, dt, method, gap_junc);
amacrine = Amacrine(a0, buffer_size, dt, method, gap_junc);
rgc = Ganglion(rgc0, buffer_size, dt, method);

offset = 42.2321233629648844498660764656960964202880859375;
cone_bip.set_Vth(cone0(1)+offset);
rod_bip.set_Vth(rod0(1));
amacrine.set_Vth(rbip0(1));
rgc.set_Vth(cbip0(1));

v_cone = cone0(1);
v_rod = rod0(1);
v_a = a0(1);
v_cbip = cbip0(1);

tic
while abs(curr_t - t_end) > eps
    
    in_rod_interp = in_rod(round(1+buffer_size/(t_end-t_start)*curr_t));
    in_cone_interp = in_cone(round(1+buffer_size/(t_end-t_start)*curr_t));
    
    [y_rod, curr_t, c_rod] = rod.solve(in_rod_interp, v_cone+offset);
    [y_cone, ~, ~] = cone.solve(in_cone_interp, v_rod-offset);
    [y_rbip, ~, ~] = rod_bip.solve(y_rod(1));
    [y_cbip, ~, ~] = cone_bip.solve(y_cone(1)+offset, v_a);
    [y_a, ~, ~] = amacrine.solve(y_rbip(1), v_cbip);
    [y_rgc, ~, ~] = rgc.solve(y_cbip(1));
    
    rod.update_time();
    cone.update_time();
    rod_bip.update_time();
    cone_bip.update_time();
    amacrine.update_time();
    rgc.update_time();
    
    v_cone = y_cone(1);
    v_rod = y_rod(1);
    v_a = y_a(1);
    v_cbip = y_cbip(1);
end
end_time = toc;

t_per_step = end_time/length(jhvt);
tot_t = end_time;

t_vec = rod.get_tvec();
t_vec_end = find(t_vec == curr_t);

v_rod = rod.get_V();
v_cone = cone.get_V();
v_rbip = rod_bip.get_V();
v_cbip = cone_bip.get_V();
v_a = amacrine.get_V();
v_rgc = rgc.get_V();

plot(ax1, jhvt, in_rod); title(ax1, 'Input to Rod Photoreceptor');
xlabel(ax1, 't [s]'); ylabel(ax1, 'I [PI/s]')
xlim(ax1, [t_vec(1), t_vec(t_vec_end)])
grid(ax1, 'on'); grid(ax1, 'minor')
hold(ax1, 'on')

plot(ax2, jhvt, in_cone); title(ax2, 'Input to Cone Photoreceptor');
xlabel(ax2, 't [s]'); ylabel(ax2, 'I [PI/s]')
xlim(ax2, [t_vec(1), t_vec(t_vec_end)])
grid(ax2, 'on'); grid(ax2, 'minor')
hold(ax2, 'on')

plot(ax3, t_vec(1:t_vec_end), v_rod(1:t_vec_end)); title(ax3, 'Rod Membrane Potential');
xlabel(ax3, 't [s]'); ylabel(ax3, 'V_{rod} [mV]')
xlim(ax3, [t_vec(1), t_vec(t_vec_end)])
grid(ax3, 'on'); grid(ax3, 'minor')
hold(ax3, 'on')

plot(ax4, t_vec(1:t_vec_end), v_cone(1:t_vec_end)); title(ax4, 'Cone Membrane Potential');
xlabel(ax4, 't [s]'); ylabel(ax4, 'V_{cone} [mV]')
xlim(ax4, [t_vec(1), t_vec(t_vec_end)])
grid(ax4, 'on'); grid(ax4, 'minor')
hold(ax4, 'on')

plot(ax5, t_vec(1:t_vec_end), v_rbip(1:t_vec_end)); title(ax5, 'RBP Membrane Potential');
xlabel(ax5, 't [s]'); ylabel(ax5, 'V_{rbp} [mV]')
xlim(ax5, [t_vec(1), t_vec(t_vec_end)])
grid(ax5, 'on'); grid(ax5, 'minor')
hold(ax5, 'on')

plot(ax6, t_vec(1:t_vec_end), v_cbip(1:t_vec_end)); title(ax6, 'CBP Membrane Potential');
xlabel(ax6, 't [s]'); ylabel(ax6, 'V_{cbp} [mV]')
xlim(ax6, [t_vec(1), t_vec(t_vec_end)])
grid(ax6, 'on'); grid(ax6, 'minor')
hold(ax6, 'on')

plot(ax7, t_vec(1:t_vec_end), v_a(1:t_vec_end)); title(ax7, 'Amacrine Membrane Potential');
xlabel(ax7, 't [s]'); ylabel(ax7, 'V_a [mV]')
xlim(ax7, [t_vec(1), t_vec(t_vec_end)])
grid(ax7, 'on'); grid(ax7, 'minor')
hold(ax7, 'on')

if j == 1
    plot(g1, t_vec(1:t_vec_end), v_rgc(1:t_vec_end))
    if test_num == 1
        title(g1, ['Ganglion Membrane Potential ', ...
            '(rod:', num2str(rod_vals(j)), ' cone: ', num2str(cone_vals(j)), ' [photoisomerization/s])']);
    elseif test_num == 2
        title(g1, ['Ganglion Membrane Potential ', '(duration = ', num2str(rod_vals(j)/100000), ' [s])']);
    end
    xlabel(g1, 't [s]'); ylabel(g1, 'V_{rgc} [mV]')
    xlim(g1, [t_vec(1), t_vec(t_vec_end)])
    grid(g1, 'on')
    grid(g1, 'minor')
elseif j == 2
    plot(g2, t_vec(1:t_vec_end), v_rgc(1:t_vec_end))
    if test_num == 1
        title(g2, ['Ganglion Membrane Potential ', ...
            '(rod:', num2str(rod_vals(j)), ' cone: ', num2str(cone_vals(j)), ' [photoisomerization/s])']);
    elseif test_num == 2
        title(g2, ['Ganglion Membrane Potential ', '(duration = ', num2str(rod_vals(j)/100000), ' [s])']);
    end
    xlabel(g2, 't [s]'); ylabel(g2, 'V_{rgc} [mV]')
    xlim(g2, [t_vec(1), t_vec(t_vec_end)])
    grid(g2, 'on')
    grid(g2, 'minor')
elseif j == 3
    plot(g3, t_vec(1:t_vec_end), v_rgc(1:t_vec_end))
    if test_num == 1
        title(g3, ['Ganglion Membrane Potential ', ...
            '(rod:', num2str(rod_vals(j)), ' cone: ', num2str(cone_vals(j)), ' [photoisomerization/s])']);
    elseif test_num == 2
        title(g3, ['Ganglion Membrane Potential ', '(duration = ', num2str(rod_vals(j)/100000), ' [s])']);
    end
    xlabel(g3, 't [s]'); ylabel(g3, 'V_{rgc} [mV]')
    xlim(g3, [t_vec(1), t_vec(t_vec_end)])
    grid(g3, 'on')
    grid(g3, 'minor')
elseif j == 4
    plot(g4, t_vec(1:t_vec_end), v_rgc(1:t_vec_end))
    if test_num == 1
        title(g4, ['Ganglion Membrane Potential ', ...
            '(rod:', num2str(rod_vals(j)), ' cone: ', num2str(cone_vals(j)), ' [photoisomerization/s])']);
    elseif test_num == 2
        title(g4, ['Ganglion Membrane Potential ', '(duration = ', num2str(rod_vals(j)/100000), ' [s])']);
    end
    xlabel(g4, 't [s]'); ylabel(g4, 'V_{rgc} [mV]')
    xlim(g4, [t_vec(1), t_vec(t_vec_end)])
    grid(g4, 'on')
    grid(g4, 'minor')
elseif j == 5
    plot(g5, t_vec(1:t_vec_end), v_rgc(1:t_vec_end))
    if test_num == 1
        title(g5, ['Ganglion Membrane Potential ', ...
            '(rod:', num2str(rod_vals(j)), ' cone: ', num2str(cone_vals(j)), ' [photoisomerization/s])']);
    elseif test_num == 2
        title(g5, ['Ganglion Membrane Potential ', '(duration = ', num2str(rod_vals(j)/100000), ' [s])']);
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

legend(ax1, labels); legend(ax2, labels); legend(ax3, labels)
legend(ax4, labels); legend(ax5, labels); legend(ax6, labels)
legend(ax7, labels)
legend(rate_fig, labels)