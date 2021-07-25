clc
close all
clear variables


v = -80 : 0.01 : 80;
% v = -80;
ab = zeros(10, length(v));

for i = 1 : length(v)
    ab(1:4, i) = Kv(v(i));
    ab(5:6, i) = Ca(v(i));
    ab(7:8, i) = K_Ca(v(i));
    ab(9:10, i) = h(v(i));
end

a = zeros(4, 1);
a(:) = Kv(-80);

figure 
subplot(2, 3, 1)
plot(v, 1./(ab(1,:)+ab(2,:)))
title('\tau_{mKv}')
subplot(2, 3, 2)
plot(v, 1./(ab(3,:)+ab(4,:)))
title('\tau_{hKv}')
subplot(2, 3, 3)
plot(v, 1./(ab(5,:)+ab(6,:)))
title('\tau_{mCa}')
subplot(2, 3, 4)
plot(v, 1./(ab(7,:)+ab(8,:)))
title('\tau_{mKCa}')
subplot(2, 3, 5)
plot(v, 1./(ab(9,:)+ab(10,:)))
title('\tau_{h}')

figure 
subplot(2, 3, 1)
plot(v, ab(1,:)./(ab(1,:)+ab(2,:)))
title('mKv_{\infty}')
subplot(2, 3, 2)
plot(v, ab(3,:)./(ab(3,:)+ab(4,:)))
title('hKv_{\infty}')
subplot(2, 3, 3)
plot(v, ab(5,:)./(ab(5,:)+ab(6,:)))
title('mCa_{\infty}')
subplot(2, 3, 4)
plot(v, ab(7,:)./(ab(7,:)+ab(8,:)))
title('mKCa_{\infty}')
subplot(2, 3, 5)
plot(v, ab(9,:)./(ab(9,:)+ab(10,:)))
title('h_{\infty}')

figure 
plot(v, 1./(ab(1,:)+ab(2,:)))
hold on
plot(v, 1./(ab(3,:)+ab(4,:)))
hold on
plot(v, 1./(ab(5,:)+ab(6,:)))
hold on
plot(v, 1./(ab(7,:)+ab(8,:)))
hold on
plot(v, 1./(ab(9,:)+ab(10,:)))
hold off
legend('\tau_{mKv}', '\tau_{hKv}', '\tau_{mCa}', '\tau_{mKCa}', '\tau_{h}')


%%%%% Kv %%%%%
function ab = Kv(V)

am_Kv = 5*(100-V)/(exp((100-V)/42)-1);
bm_Kv = 9*exp(-(V-20)/40);
ah_Kv = 0.15*exp(-V/22);
bh_Kv = 0.4125/(exp((10-V)/7)+1);
ab = [am_Kv; bm_Kv; ah_Kv; bh_Kv];

end
%%%%%%%%%%%%%%%%

%%%%% Ca %%%%%
function ab = Ca(V)

am_Ca = 3*(80-V)/(exp((80-V)/25)-1);
bm_Ca = 10/(1+exp((V+38)/7));
ab = [am_Ca; bm_Ca];

end
%%%%%%%%%%%%%%%%

%%%%% K_Ca %%%%%
function ab = K_Ca(V)

am_KCa = 15*(80-V)/(exp((80-V)/40)-1);
bm_KCa = 20*exp(-V/35);
ab = [am_KCa; bm_KCa];

end
%%%%%%%%%%%%%%%%

%%%%% h %%%%%
function ab = h(V)

ah = 8/(exp((V+78)/14)+1);
bh = 18/exp(-(V+8)/19+1);
ab = [ah, bh];

end
%%%%%%%%%%%%%%%%