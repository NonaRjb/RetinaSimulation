clc
clear variables
close all

%%% initial values for the bipolar cell
% [V, mKv, hKv, mA, hA, C1, C2, O1, O2, O3, mCa, mKCa, Cas, Cad, Cabls, 
%  Cabhs, Cabld, Cabhd]
bip0 = [-36.424516776897130, 0.824374082910590, 0.109794106890060, ...
    0.186127778073054, 0.024443585120781, 0.928238175043767, 0.054905344261992,...
    0.001217870414180, 1.200618479085905e-05, 4.438540983730620e-08, ...
    0.290203193088928, 0.475464147362579, 0.011561930331728, 0.011563608687641,...
    6.780371247710756, 1.268364765067093, 11.302574980627850, 3.805639865311822, 0];

buffer_size = 100000; % if set to 1000000 the result will be accurate but needs more time
t_start = 0;
t_end = 1;
bip = Bipolar_complete(bip0, buffer_size);
jhvt = linspace(0,t_end,buffer_size);

%%% input sample #1
jhv = (10)*ones(size(jhvt));

%%% input sample #2
%jhv = zeros(size(jhvt));
%jhv(10000:end) = 10;

%%% input sample #3
% jhv = zeros(size(jhvt));
% jhv(100000:102000) = 100;

%iPhoto = zeros(size(jhvt))';
curr_t_bip = t_start;
%ab = zeros(10, buffer_size);
%i = 1;
%maxD = [];

tic
while curr_t_bip < t_end && jhvt(end)>=curr_t_bip
    input_j = interp1(jhvt, jhv, curr_t_bip);
    %[y_rod, curr_t_rod, c_rod] = rod.solve(input_j);
    %rod.update_time();
    [y_bip, curr_t_bip, c_bip] = bip.solve(input_j);
    bip.update_time();
    %maxD = [maxD, objdt];
    %iPhoto(i) = c(end);
    %ab(:, i) = c(10:19);
    %i = i + 1;
end
end_time = toc;

t_per_step = end_time/length(jhvt);
tot_t = end_time;

t_vec = bip.get_tvec();
t_vec_end = find(t_vec == curr_t_bip);

v = bip.get_V();
figure
plot(t_vec(1:t_vec_end), v(1:t_vec_end)); title('Membrane Potential'); xlabel('t [s]'); ylabel('V_m [mV]')
xlim([t_vec(1), t_vec(t_vec_end)])
grid on
grid minor