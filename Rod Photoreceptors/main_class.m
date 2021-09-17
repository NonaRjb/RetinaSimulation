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

test_num = input('please enter the test number: ');
dt = 2*1e-05;
eps = dt;

%iPhoto = zeros(size(jhvt))';
%curr_t_rod = t_start;
%ab = zeros(10, buffer_size);
%i = 1;

figure
if test_num == 1
    buffer_size = 60000;
    t_start = 0;
    t_end = 0.6;
    j_vals = [10^(0), 10^(0.6), 10^(1.2), 10^(1.8), 10^(2.4), 10^(3.0), 10^(3.6)];
    t_dur = 1300;
    labels = {'10^0', '10^{0.6}', '10^{1.2}', '10^{1.8}', '10^{2.4}', '10^{3.0}',...
    '10^{3.6}'};
elseif test_num == 2
    buffer_size = 1200000;
    t_start = 0;
    t_end = 12;
    j_vals = [10^(2.4), 10^(3.0), 10^(3.6)];
    t_dur = 1300;
    labels = {'10^{2.4}', '10^{3.0}', '10^{3.6}'};
elseif test_num == 3
    buffer_size = 60000;
    t_start = 0;
    t_end = 0.6;
    j_vals = [10^(0.6), 10^(1.2), 10^(1.8), 10^(2.4), 10^(3.0), 10^(3.6), ...
        10^(4.2), 10^(4.8), 10^(5.4)];
    t_dur = buffer_size;
    labels = {'10^{0.6}', '10^{1.2}', '10^{1.8}', '10^{2.4}', '10^{3.0}',...
    '10^{3.6}', '10^{4.2}', '10^{4.8}', '10^{5.4}'};
elseif test_num == 4
    buffer_size = 1200000;
    t_start = 0;
    t_end = 12;
    j_vals = [10^(2.4), 10^(3.0), 10^(3.6), 10^(4.2), 10^(4.8), 10^(5.4)];
    t_dur = 100000;
    abels = {'10^{2.4}', '10^{3.0}', '10^{3.6}', '10^{4.2}', '10^{4.8}', '10^{5.4}'};
end
jhvt = linspace(0,t_end,buffer_size);

for j = 1 : length(j_vals)
jhv = zeros(size(jhvt));
jhv(1:t_dur) = j_vals(j);
curr_t_rod = t_start;
rod = RodPhotoReceptor_RK(rod0, buffer_size, dt);
% tic
while abs(curr_t_rod - t_end) > eps
    input_j = jhv(round(1+buffer_size/(t_end-t_start)*curr_t_rod));
    [y_rod, curr_t_rod, ~] = rod.solve(input_j);
    rod.update_time();
    %iPhoto(i) = c(end);
    %ab(:, i) = c(10:19);
    %i = i + 1;
end
%end_time = toc;

%t_per_step = end_time/length(jhvt);
%tot_t = end_time;

t_vec = rod.get_tvec();
t_vec_end = find(t_vec == curr_t_rod);

v = rod.get_V();
%figure
plot(t_vec(1:t_vec_end), v(1:t_vec_end)-rod0(1)); title('Membrane Potential'); xlabel('t [s]'); ylabel('V_m [mV]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor

hold on
end
legend(labels, 'Location', 'southeast');
hold off
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

% cas = rod.get_Cas();
% figure 
% plot(t_vec(1:t_vec_end), cas(1:t_vec_end));  title('[Ca_s]'); xlabel('t [s]'); ylabel('[Ca_s] [uM]')
% xlim([t_vec(1), t_vec(t_vec_end)])
% grid on
% grid minor

