close all;
clear all;
%% 1.2
filt = [ones(1,101),zeros(1,400),ones(1,100)];

impres=ifft(filt);

impres=[impres(301:end),impres(1:300)]; %flip impulse response
n=290;
n2=100;

window = [zeros(1,n),ones(1,601-2*n),zeros(1,n)];   %windowing functions
window2 = [zeros(1,n2),ones(1,601-2*n2),zeros(1,n2)];
windowed=impres.*window;
windowed2=impres.*window2;

% figure 1 

figure('DefaultAxesFontSize',16)
subplot(3,1,1)
plot(impres,'LineWidth',1.5,'color','#0072BD')

%beautify
set(gca,'XMinorTick','on','YMinorTick','on')
axis([301 601 -0.1 0.35])
yticks([ 0 0.2])
ylabel('Amplitude ', 'FontSize', 16)
title('n=0','FontSize', 16)

subplot(3,1,2)
plot(windowed2,'LineWidth',1.5,'color','#D95319')

%beautify
set(gca,'XMinorTick','on','YMinorTick','on')
axis([301 601 -0.1 0.35])
yticks([ 0 0.2])
ylabel('Amplitude ', 'FontSize', 16)
title('n=100','FontSize', 16)


subplot(3,1,3)
plot(windowed,'LineWidth',1.5,'color','#77AC30')

%beautify
set(gca,'XMinorTick','on','YMinorTick','on')
axis([301 601 -0.1 0.35])
yticks([ 0 0.2])
ylabel('Amplitude ', 'FontSize', 16)
title('n=290','FontSize', 16)


sgtitle('Impulse response','FontSize', 16)
xlabel('Samples','FontSize', 16) 


%figure2
figure('DefaultAxesFontSize',16)

plot((abs(fft(windowed2))),'LineWidth',1.5,'color','#D95319')

%beautify
set(gca,'XMinorTick','on','YMinorTick','on')
axis([0 301 0 1.15])

hold on 

plot((abs(fft(windowed))),'LineWidth',1.5,'color','#77AC30')

%beautify
set(gca,'XMinorTick','on','YMinorTick','on')
axis([0 301 0 1.15])



plot(filt,'LineWidth',1.5,'color','#0072BD')
set(gca,'XMinorTick','on','YMinorTick','on')
xticks([0 101, 201, 301])
xticklabels({'0','5','10','15'})
legend('Filter after windowing n=100','Filter after windowing n=290','Original filter')

title('Frequency reponse','FontSize', 16)
xlabel('Frequency [kHz]','FontSize', 16) %positive and negative time
ylabel('Amplitude ', 'FontSize', 16)