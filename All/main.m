clc
clear variables
close all


% rod0 = [-36.186, 0.430, 0.999, 0.436, 0.642, 0.646, 0.298, 0.0517, 0.00398, ...
%     0.000115, 0.0966, 0.0966, 80.929, 29.068, 80.929, 29.068, 0, 0, 0, 0, ...
%     0.3, 34.88, 2.0];

rod0 = [-36.185963, 0.43018469, 0.99927789, 0.43647161, 0.64228624,...
    0.64582664, 0.29823381, 0.051645089, 0.0039748312, 0.00011472013,...
    0.096557982, 0.096558392, 80.92926, 29.067444, 80.929703, 29.067556,...
    0.0, 0.0, 0.0, 0.0, 0.300000000000037, 34.883720929940061, 2.000000000017296];

% rod0 = [-36.186, 0.999, 0.436, 0.646, 0.298, 0.0517, 0.00398, ...
%     0.000115, 0.0966, 0.0966, 80.929, 29.068, 80.929, 29.068, 0, 0, 0, 0, ...
%     0.3, 34.88, 2.0];


%%%
bip0 = [-36.424516776897130, 0.824374082910590, 0.109794106890060, ...
    0.186127778073054, 0.024443585120781, 0.928238175043767, 0.054905344261992,...
    0.001217870414180, 1.200618479085905e-05, 4.438540983730620e-08, ...
    0.290203193088928, 0.475464147362579, 0.011561930331728, 0.011563608687641,...
    6.780371247710756, 1.268364765067093, 11.302574980627850, 3.805639865311822, 0];

buffer_size = 100000; % if set to 1000000 the result will be accurate but needs more time
t_start = 0;
t_end = 1;
rod = RodPhotoReceptor_RK(rod0, buffer_size);
bip = Bipolar_complete(bip0, buffer_size);
jhvt = linspace(0,t_end,buffer_size);