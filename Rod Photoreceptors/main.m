clc
%close all
%clear variables


tspan = [0, 100];
Y0 = [-36.186, 0.430, 0.999, 0.436, 0.642, 0.646, 0.298, 0.0517,...
    0.00398, 0.000115, 0.0966, 0.0966, 80.929, 29.068, 80.929, 29.068, 0, 0, 0, 0, 0.3, 34.88, 2.0];
jhvt = linspace(0,100,10000000);
jhv = (10)*ones(size(jhvt));
%jhv(60000:end) = 0;
%jhv = zeros(size(jhvt));
%jhv(10000:10400) = 100;
tic
[t, y] = ode15s(@(t, Y) odefuncs_rod(t, Y, jhvt, jhv), tspan, Y0);
end_time = toc;

t_per_step = end_time/length(jhvt);
tot_t = end_time;

jmax = 5040;
j = jmax*(y(:,23)).^3./(y(:,23).^3+10^3);
iPhoto = -j.*(1.0-exp((y(:,1)-8.5)/17.0));
figure
subplot(3, 1, 1)
plot(jhvt, jhv); title('Input Stimulus'); xlabel('t [s]'); ylabel('J_{hv}(t) [Rh^*.s^{-1}] ')
grid on
grid minor
subplot(3, 1, 2)
plot(t, iPhoto); title('Photocurrent'); xlabel('t [s]'); ylabel('I_{Photo}(t) [pA]')
grid on
grid minor
subplot(3, 1, 3)
plot(t,y(:, 1)); title('Membrane Potential'); xlabel('t [s]'); ylabel('V_m [mV]')
grid on
grid minor

figure
plot(t, y(:, 11))
grid on 
grid minor