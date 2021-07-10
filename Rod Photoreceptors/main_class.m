clc
%clear variables
%close all


Y0 = [-36.186, 0.430, 0.999, 0.436, 0.642, 0.646, 0.298, 0.0517, 0.00398, ...
    0.000115, 0.0966, 0.0966, 80.929, 29.068, 80.929, 29.068, 0, 0, 0, 0, ...
    0.3, 34.88, 2.0];

buffer_size = 100000;
rod = RodPhotoReceptor(Y0, buffer_size);
jhvt = linspace(0,1,100000);
jhv = (10)*ones(size(jhvt));
%jhv(10000:end) = 0;
%jhv = zeros(size(jhvt));
%jhv(10000:10400) = 100;

iPhoto = zeros(size(jhvt))';

tic
for i = 1 : length(jhvt)
    [y, c] = rod.solve(jhv(i));
    rod.update_time();
    iPhoto(i) = c(end);
end
end_time = toc;

t_per_step = end_time/length(jhvt);
tot_t = end_time;

v = rod.get_V();
figure
plot(jhvt, v)
grid on

cas = rod.get_Cas();
figure 
plot(jhvt, cas)
grid on

d = rod.get_diff();
