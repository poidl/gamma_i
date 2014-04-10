% figure(5), subplot(2,1,1), plot(log_c(:),z1'), grid on
% xlabel('log(c)'), ylabel('mean log(D_V)')
% figure(5), subplot(2,1,2), plot(log_c(:),z2'), grid on
% xlabel('log(c)'), ylabel('%')

figure(5), subplot(2,1,1), plot(handles_maxitss,z1), grid on
xlabel('# iterations'), ylabel('mean log(D_V)')
figure(5), subplot(2,1,2), plot(handles_maxitss,z2), grid on
xlabel('# iterations'), ylabel('%')