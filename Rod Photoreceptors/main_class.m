clc
clear variables
%close all


rod0 = [-36.186, 0.430, 0.999, 0.436, 0.642, 0.646, 0.298, 0.0517, 0.00398, ...
    0.000115, 0.0966, 0.0966, 80.929, 29.068, 80.929, 29.068, 0, 0, 0, 0, ...
    0.3, 34.88, 2.0];

% rod0 = [-36.186, 0.999, 0.436, 0.646, 0.298, 0.0517, 0.00398, ...
%     0.000115, 0.0966, 0.0966, 80.929, 29.068, 80.929, 29.068, 0, 0, 0, 0, ...
%     0.3, 34.88, 2.0];


%%%
% bip0 = [-36.186, 0.430, 0.999, 0, 0, 0.646, 0.298, 0.0517, 0.00398,...
%     0.000115, 0.436, 0.642, 0.0966, 0, 80.929, 29.068, 0, 0, 0];

buffer_size = 2; % if set to 1000000 the result will be accurate but needs more time
t_start = 0;
t_end = 10;
rod = RodPhotoReceptor_RK(rod0, buffer_size);
%bip = Bipolar_complete(bip0, buffer_size);
jhvt = linspace(0,t_end,1000000);

%%% input sample #1
%jhv = (10)*ones(size(jhvt));

%%% input sample #2
%jhv = zeros(size(jhvt));
%jhv(10000:end) = 10;

%%% input sample #3
jhv = zeros(size(jhvt));
jhv(100000:102000) = 100;

%iPhoto = zeros(size(jhvt))';
curr_t_rod = t_start;
%ab = zeros(10, buffer_size);
%i = 1;
%maxD = [];

tic
figure
h = animatedline;
axis([0 t_end -60 -35])
addpoints(h,curr_t_rod,rod0(1));
drawnow
while curr_t_rod < t_end && jhvt(end)>=curr_t_rod
    input_j = interp1(jhvt, jhv, curr_t_rod);
    [y_rod, curr_t_rod, ~] = rod.solve(input_j);
    rod.update_time();
    addpoints(h,curr_t_rod,y_rod(1));
    drawnow limitrate
%     [y_bip, curr_t_bip, c_bip] = bip.solve(y_rod(1));
%     bip.update_time();
    %maxD = [maxD, objdt];
    %iPhoto(i) = c(end);
    %ab(:, i) = c(10:19);
    %i = i + 1;
end
drawnow
end_time = toc;

t_per_step = end_time/length(jhvt);
tot_t = end_time;

t_vec = rod.get_tvec();
t_vec_end = find(t_vec == curr_t_rod);

v = rod.get_V();
figure
plot(t_vec(1:t_vec_end), v(1:t_vec_end)); title('Membrane Potential'); xlabel('t [s]'); ylabel('V_m [mV]')
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

