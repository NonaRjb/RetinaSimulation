clc
clear variables
close all


%%% Cone initial values %%%
cone0 = [-42.2321233629648844498660764656960964202880859375,...
    0.004863366238889023189517768486211934941820800304412841796875,...
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

buffer_size = 1000000; % if set to 1000000 the result will be accurate but needs more time
t_start = 0;
t_end = 1;
dt = 2*1e-05;
eps = dt;
cone = ConePhotoReceptor(cone0, buffer_size, dt);
jhvt = linspace(0,t_end,buffer_size);

%%% input sample #1
% jhv = (10)*ones(size(jhvt));

%%% input sample #2
% jhv = zeros(size(jhvt));
% jhv(500000:end) = 10;

%%% input sample #3
jhv = zeros(size(jhvt));
jhv(10000:20000) = 890000;

%iPhoto = zeros(size(jhvt))';
curr_t_cone = t_start;
%ab = zeros(10, buffer_size);
%i = 1;
%maxD = [];

tic
while abs(curr_t_cone - t_end) > eps
%for i = 1 : buffer_size
    % input_j = interp1(jhvt, jhv, curr_t_rod);
    input_j = jhv(round(1+buffer_size/(t_end-t_start)*curr_t_cone));
    [y_rod, curr_t_cone, c_rod] = cone.solve(input_j);
    cone.update_time();
    %maxD = [maxD, objdt];
    %iPhoto(i) = c_rod(end);
    %ab(:, i) = c(10:19);
    %i = i + 1;
end
end_time = toc;

t_per_step = end_time/length(jhvt);
tot_t = end_time;

t_vec = cone.get_tvec();
t_vec_end = find(t_vec == curr_t_cone);

v_cone = cone.get_V();

figure
subplot(2,1,1)
plot(jhvt, jhv); title('Input Photocurrent'); xlabel('t [s]'); ylabel('I [Rh/s]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
subplot(2,1,2)
plot(t_vec(1:t_vec_end), v_cone(1:t_vec_end)); title('Cone Membrane Potential'); xlabel('t [s]'); ylabel('V_m [mV]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor

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

