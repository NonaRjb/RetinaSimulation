clc
clear variable
close all


h0 = [-60, 0.03, 0.998, 0.139, 0.059, 0.026, 0.922, 0, 0, 0];

buffer_size = 1000000; % if set to 1000000 the result will be accurate but needs more time
t_start = 0;
t_end = 10;
dt = 2*1e-05;
eps = dt;
%rod = RodPhotoReceptor_RK(rod0, buffer_size, dt);
h = Horizontal(h0, buffer_size, dt);
jhvt = linspace(0,t_end,buffer_size);

%%% input sample #1
%jhv = 10*ones(size(jhvt));

%%% input sample #2
jhv = zeros(size(jhvt));
jhv(100000:end) = 100;

%%% input sample #3
%jhv = zeros(size(jhvt));
%jhv(100000:102000) = 100;
%f0 = 0.5;
%jhv = sin(2*pi*f0*jhvt);

%iPhoto = zeros(size(jhvt))';
curr_t_h = t_start;
%ab = zeros(10, buffer_size);
i = 1;
%maxD = [];
glu = zeros(buffer_size, 1);
iGlu = zeros(buffer_size, 1);

tic
while abs(curr_t_h - t_end) > eps
%for i = 1 : buffer_size
    % input_j = interp1(jhvt, jhv, curr_t_rod);
    input_j = jhv(round(1+buffer_size/(t_end-t_start)*curr_t_h));
    [~, curr_t_h, c] = h.solve(input_j);
    h.update_time;
    glu(i) = c(8);
    iGlu(i) = c(6);
    i = i + 1;
end
end_time = toc;

t_per_step = end_time/length(jhvt);
tot_t = end_time;

t_vec = h.get_tvec();
t_vec_end = find(t_vec == curr_t_h);

%v_rod = rod.get_V();
v_h = h.get_V();
gaba = h.get_gaba();


figure
subplot(4,1,1)
plot(jhvt, jhv); title('Input Photocurrent'); xlabel('t [s]'); ylabel('I [Rh/s]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
subplot(4,1,2)
plot(t_vec(1:t_vec_end), glu(1:t_vec_end)); title('Glutamate concentration'); xlabel('t [s]'); ylabel('V_m [mV]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
subplot(4,1,3)
plot(t_vec(1:t_vec_end), iGlu(1:t_vec_end)); title('Glutamte-induced current'); xlabel('t [s]'); ylabel('V_m [mV]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
subplot(4,1,4)
plot(t_vec(1:t_vec_end), v_h(1:t_vec_end)); title('Horizontal Membrane Potential'); xlabel('t [s]'); ylabel('S(t) [mV]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
% subplot(3,2,4)
% plot(t_vec(1:t_vec_end), v_bip(1:t_vec_end)); title('RBP Membrane Potential'); xlabel('t [s]'); ylabel('V_m [mV]')
% xlim([t_vec(1), t_vec(t_vec_end)])
% grid on
% grid minor
% subplot(3,2,5)
% plot(t_vec(1:t_vec_end), s_rgc(1:t_vec_end)); title('Synaptic Value'); xlabel('t [s]'); ylabel('S(t) [mV]')
% xlim([t_vec(1), t_vec(t_vec_end)])
% grid on
% grid minor
% subplot(3,2,6)
% plot(t_vec(1:t_vec_end), v_rgc(1:t_vec_end)); title('RGC Membrane Potential'); xlabel('t [s]'); ylabel('V_m [mV]')
% xlim([t_vec(1), t_vec(t_vec_end)])
% grid on
% grid minor