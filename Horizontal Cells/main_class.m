 clc
clear variable
close all


%%% Horizontal Cell initial values %%%
% [h_V, h_mA, h_hA, h_mKv, h_mCa, h_mNa, h_hNa, h_B, h_Ca, h_Gaba]
h0 = [-70, 0.03, 0.998, 0.139, 0.059, 0.026, 0.922, 0, 0.05, 0.1];

%%% Horizontal Cell 2 initial values %%%
% h0 = [-80, 0.026, 0.922, 0.059, 0.139, 0.139, 0.932, 0.03, 0.998]; 

buffer_size = 10000000; % if set to 1000000 the result will be accurate but needs more time
t_start = 0;
t_end = 150;
dt = 2*1e-05;
eps = dt;
%rod = RodPhotoReceptor_RK(rod0, buffer_size, dt);
h = Horizontal(h0, buffer_size, dt);
jhvt = linspace(0,t_end,buffer_size);

%%% input sample #1
%jhv = -36*ones(size(jhvt));
%jhv = 10+(10-36)*exp(-jhvt/10);

%%% input sample #2
jhv = zeros(size(jhvt));
%jhv(100000:end) = 1;

%%% input sample #3
%jhv = zeros(size(jhvt));
%jhv(100000:102000) = 100;
%f0 = 0.5;
%jhv = sin(2*pi*f0*jhvt);

%jhv = importdata('v_rod.mat');

%iPhoto = zeros(size(jhvt))';
curr_t_h = t_start;
%ab = zeros(10, buffer_size);
i = 1;
%maxD = [];
glu = zeros(buffer_size, 1);
iGlu = zeros(buffer_size, 1);
iCa = zeros(buffer_size, 1);
D8 = zeros(buffer_size, 1);
dCain = zeros(buffer_size, 1);
dCaef = zeros(buffer_size, 1);

tic
while abs(curr_t_h - t_end) > eps
%for i = 1 : buffer_size
    % input_j = interp1(jhvt, jhv, curr_t_rod);
    input_j = jhv(round(1+buffer_size/(t_end-t_start)*curr_t_h));
    [~, curr_t_h, c] = h.solve(input_j);
    h.update_time();
    % glu(i) = c(8);
    iGlu(i) = c(6);
    iCa (i) = c(1);
    D8(i) = c(9);
    dCain(i) = c(10);
    dCaef(i) = c(11);
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
ca = h.get_Ca();
sb = h.get_CaStore();


figure
subplot(3,2,1)
plot(t_vec(1:t_vec_end), iCa(1:t_vec_end)); title('Calcium current'); xlabel('t [s]'); ylabel('I_{Ca} [pA]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
subplot(3,2,2)
plot(t_vec(1:t_vec_end), sb(1:t_vec_end)); title('Calcium Store concentration'); xlabel('t [s]'); ylabel('[SB] [\muM]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
subplot(3,2,3)
plot(t_vec(1:t_vec_end), iGlu(1:t_vec_end)); title('Glutamte-induced current'); xlabel('t [s]'); ylabel('V_m [mV]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
subplot(3,2,4)
plot(t_vec(1:t_vec_end), v_h(1:t_vec_end)); title('Horizontal Membrane Potential'); xlabel('t [s]'); ylabel('S(t) [mV]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
subplot(3,2,5)
plot(t_vec(1:t_vec_end), gaba(1:t_vec_end)); title('GABA concentration'); xlabel('t [s]'); ylabel('GABA(t)')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
subplot(3,2,6)
plot(t_vec(1:t_vec_end), ca(1:t_vec_end)); title('Ca concentration'); xlabel('t [s]'); ylabel('Ca_i(t) [uM]')
xlim([t_vec(1), t_vec(t_vec_end)])
%ylim([0, 2])
grid on
grid minor

figure
subplot(3, 1, 1)
plot(t_vec(1:t_vec_end), D8(1:t_vec_end)); title('D(8)')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
subplot(3, 1, 2)
plot(t_vec(1:t_vec_end), dCain(1:t_vec_end)); title('dCain')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor
subplot(3, 1, 3)
plot(t_vec(1:t_vec_end), dCaef(1:t_vec_end)); title('dCaef')
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