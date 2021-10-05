clc 
clear variables 
close all

%%% RGC initial values %%%
% [g_V, g_mNa, g_hNa, g_nCa, g_Cai, g_nK, g_mA, g_hA, bpg_S]
rgc0 = [-65.698524, 0.027940, 0.887161, 0.003019, 0.100007, 0.107809, ...
    0.070467, 0.300453, 0.0];

test_num = input('please enter the test number: ');
dt = 2*1e-05;
eps = dt;

%%% input sample #1
%jhv = (20)*ones(size(jhvt));

%%% input sample #2
%jhv = zeros(size(jhvt));
%jhv(10000:end) = 10;

%%% input sample #3
%jhv = zeros(size(jhvt));
%jhv(100000:102000) = 100;

%ab = zeros(10, buffer_size);
%i = 1;

if test_num == 1
    buffer_size = 350000;
    t_start = 0;
    t_end = 0.35;
    j_vals = [10, 20, 30, 40, 50];
else
    error('No such a test number')
end
jhvt = linspace(0,t_end,buffer_size);

figure
for j = 1 : length(j_vals)
rgc = Ganglion(rgc0, buffer_size, dt);
curr_t_rgc = t_start;
jhv = (j_vals(j))*ones(size(jhvt));
% tic
while abs(curr_t_rgc - t_end) > eps
    input_j = jhv(round(1+buffer_size/(t_end-t_start)*curr_t_rgc));
    [y_rgc, curr_t_rgc, c_rgc] = rgc.solve(input_j);
    rgc.update_time();
end
% end_time = toc;

% t_per_step = end_time/length(jhvt);
% tot_t = end_time;

t_vec = rgc.get_tvec();
t_vec_end = find(t_vec == curr_t_rgc);

v_rgc = rgc.get_V();
s_rgc = rgc.get_S();

subplot(length(j_vals), 1, j)
plot(t_vec(1:t_vec_end), v_rgc(1:t_vec_end))
title(['RGC Membrane Potential I = ', num2str(j_vals(j)), ' pA']);
xlabel('t [s]'); ylabel('V_m [mV]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
end
