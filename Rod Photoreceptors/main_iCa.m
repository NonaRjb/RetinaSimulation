clc
clear variables
% close all


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

buffer_size = 1000000;
t_start = 0;
t_end = 10;
dt = 2*1e-05;
eps = dt;
rod = RodPhotoReceptor_RK(rod0, buffer_size, dt);
jhvt = linspace(0,t_end,buffer_size);

jhv = zeros(size(jhvt));
jhv(100000:102000) = 100;

%iPhoto = zeros(size(jhvt))';
iCa = zeros(buffer_size, 1);
curr_t_rod = t_start;
%ab = zeros(10, buffer_size);
i = 1;

tic
while abs(curr_t_rod - t_end) > eps
    input_j = jhv(round(1+buffer_size/(t_end-t_start)*curr_t_rod));
    [y_rod, curr_t_rod, c] = rod.solve(input_j);
    rod.update_time();
    %iPhoto(i) = c(end);
    iCa(i) = c(10);
    %ab(:, i) = c(10:19);
    i = i + 1;
end
end_time = toc;

t_per_step = end_time/length(jhvt);
tot_t = end_time;

t_vec = rod.get_tvec();
t_vec_end = find(t_vec == curr_t_rod);

v = rod.get_V();
cas = rod.get_Cas();

figure

subplot(3, 1, 1)
plot(t_vec(1:t_vec_end), v(1:t_vec_end)-rod0(1)); title('Membrane Potential'); xlabel('t [s]'); ylabel('V_m [mV]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor

subplot(3, 1, 2)
plot(t_vec(1:t_vec_end), cas(1:t_vec_end));  title('[Ca_s]'); xlabel('t [s]'); ylabel('[Ca_s] [uM]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor

subplot(3, 1, 3)
plot(t_vec(1:t_vec_end), iCa(1:t_vec_end));  title('I_{Ca}'); xlabel('t [s]'); ylabel('I_{Ca} [pA]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor

