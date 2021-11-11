%% Assignment 2
%% 1.1

clc; close all; clear;

%% Butterworth filter design

fs = 44100; %Sampling frequency. Normally 44100 in sound signals
Wp = [0.2 0.3]; %Passband from 0.2 to 0.3
Ws = [0.1 0.4]; %Passband from 0.1 to 0.4
Rp = 2; % 2 dB ripple
Rs = 100; %Attenuation in stop-band at least -100 db

[n,Wn] = buttord(Wp,Ws,Rp,Rs); % Calculations the filter order and cutt-off frequencies
[b,a] = butter(n,Wn); %Calculations of the transer function coefficients


figure;
freqz(b,a);

figure;
[h,w] = freqz(b,a);
plot(w/pi,20*log10(abs(h)),'Linewidth', 3)
xlabel('Normalized Frequency (\times \pi rad/sample)')
ylabel('Magnitude (dB)')


legend('n=12,wc1=0.20, wc2=0.30');

title('Butterworth')

currentHandle = gca; % axis handle
set(currentHandle, 'FontSize', 20)
grid on


figure;
impz(b,a) %Impulse response

figure;
phasez(b,a,1000) %Phase response

%fvtool(b,a,'polezero')

[z,p,k] = tf2zp(b,a); % Find poles from transfer function coefficients

figure;
zplane(b,a) % Poles and zeros in the z-plane


% Approximations orden: 12


%% Chebyshev Type I filter design

Wp = [0.2 0.3];
Ws = [0.1 0.4];
Rp = 2;
Rs = 100;

[n1,Wp1] = cheb1ord(Wp,Ws,Rp,Rs);

[b1,a1] = cheby1(n1,1,Wp);

figure;
freqz(b1,a1)

figure;
impz(b1,a1)

figure;
phasez(b1,a1,n1)


[z1,p1,k1] = tf2zp(b1,a1);

figure;
zplane(z1,p1)


%% Plots Butterworth and Chev I


[h,w] = freqz(b,a);
[h1,w1] = freqz(b1,a1);

figure;
subplot(2,1,1);
plot(w/pi,20*log10(abs(h)),'Linewidth', 1.5)
xlabel('Normalized Frequency (\times \pi rad/sample)')
ylabel('Magnitude (dB)')
legend('n=12,wc1=0.199, wc2=0.302');
title('Butterworth')
currentHandle = gca; % axis handle
set(currentHandle, 'FontSize', 17)
grid on

subplot(2,1,2)
plot(w/pi,20*log10(abs(h1)),'Linewidth', 1.5)
xlabel('Normalized Frequency (\times \pi rad/sample)')
ylabel('Magnitude (dB)')
legend('n=8,wc1=0.200, wc2=0.300');
title('Chebyshev I')
currentHandle = gca; % axis handle
set(currentHandle, 'FontSize', 17)
grid on

figure;
subplot(2,1,1);
impz(b,a)
currentHandle = gca; % axis handle
set(currentHandle, 'FontSize', 17)
xlim([0 500])
title('Butterworth')
grid on


subplot(2,1,2)
impz(b1,a1)
currentHandle = gca; % axis handle
set(currentHandle, 'FontSize', 17)
xlim([0 500])
title('Chebyshev I')
grid on

figure;
subplot(2,1,1);
zplane(b,a)
currentHandle = gca; % axis handle
set(currentHandle, 'FontSize', 17)
title('Butterworth')
grid on


subplot(2,1,2)
zplane(b1,a1)
currentHandle = gca; % axis handle
set(currentHandle, 'FontSize', 17)
title('Chebyshev I')
grid on






%% Chebyshev Type II filter design

Wp = [0.2 0.3];
Ws = [0.1 0.4];
Rp = 2;
Rs = 100;

[n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);

[b,a] = cheby2(n,Rs,Ws);

figure;
freqz(b,a)

figure;
impz(b,a)

figure;
phasez(b,a,n)


[z,p,k] = tf2zp(b,a);

figure;
zplane(z,p)

%% Elliptic filter design

Wp = [0.2 0.3];
Ws = [0.1 0.4];
Rp = 2;
Rs = 100;

[n,Wn] = ellipord(Wp,Ws,Rp,Rs);

[b,a] = ellip(n,Rp,Rs,Wp);

figure;
freqz(b,a)

figure;
impz(b,a)

figure;
phasez(b,a,n)

%figure;
%fvtool(b,a)

[z,p,k] = tf2zp(b,a);

figure;
zplane(z,p)



%% Windowing with Rect. Effective length
% Run butterworth code first before continuing

figure;
impz(b,a) % Check the impulse respinse

[h,t]=impz(b,a); % Save the impulse response in t and h

maximp=max(abs(h)); % Find maximum amplitude value

cut=maximp*0.1; %Find the cut-value. 10% of max

index2=142; %The index which makes the impulse response less than 10% is read from the impulse response

w = rectwin(index2); %Rect-signal is made

sig_w = [w.*h(1:index2); zeros((length(t)-index2),1)]; %The signal is multiplied with the rect-signal. The rest is zero-padded

figure;
stem(t,sig_w,'filled') % The cutted impulse response

delta_f= (fs)/(length(sig_w)); %sampling frequency / number of samples

sig_wfreq=fft(sig_w); % Fourier transform of the cutted signal

Newsig_wfreq=fftshift(sig_wfreq); %shifts the frequencies so they are alligned around w=0

n=length(sig_wfreq); 

fshift = (-n/2:n/2-1)*(fs/n); %Shifts the frequencies

figure;
plot(fshift,mag2db(abs(Newsig_wfreq))); % Plots the frequency response for the cutted impulse response
grid on




%% Windowing with Rect. 75%,60%,40% and 10% of the effective length


ha1=round(index2*0.75);
ha2=round(index2*0.60);
ha3=round(index2*0.40);
ha4=round(index2*0.10);

w1= rectwin(ha1);
w2= rectwin(ha2);
w3= rectwin(ha3);
w4= rectwin(ha4);

sig_w1 = [w1.*h(1:ha1); zeros((length(t)-ha1),1)];
sig_w2 = [w2.*h(1:ha2); zeros((length(t)-ha2),1)];
sig_w3 = [w3.*h(1:ha3); zeros((length(t)-ha3),1)];
sig_w4 = [w4.*h(1:ha4); zeros((length(t)-ha4),1)];

figure;
stem(t,sig_w1,'filled')
figure;
stem(t,sig_w2,'filled')
figure;
stem(t,sig_w3,'filled')
figure;
stem(t,sig_w4,'filled')
figure;
stem(t,sig_w,'filled')


sig_wfreq1=fft(sig_w1);
sig_wfreq2=fft(sig_w2);
sig_wfreq3=fft(sig_w3);
sig_wfreq4=fft(sig_w4);
sig_wfreq=fft(sig_w);
sig_wfreqTRUE=fft(h);


Newsig_wfreq1=fftshift(sig_wfreq1);
Newsig_wfreq2=fftshift(sig_wfreq2);
Newsig_wfreq3=fftshift(sig_wfreq3);
Newsig_wfreq4=fftshift(sig_wfreq4);
Newsig_wfreq=fftshift(sig_wfreq);
Newsig_wfreqTRUE=fftshift(sig_wfreqTRUE);

n=length(sig_wfreq);

fshift = (-n/2:n/2-1)*(fs/n);

figure;
hold on
plot(fshift,mag2db(abs(Newsig_wfreqTRUE)),'Linewidth',1.5);
plot(fshift,mag2db(abs(Newsig_wfreq)),'Linewidth',1.5);
plot(fshift,mag2db(abs(Newsig_wfreq1)),'Linewidth',1.5);
plot(fshift,mag2db(abs(Newsig_wfreq2)),'Linewidth',1.5);
plot(fshift,mag2db(abs(Newsig_wfreq3)),'Linewidth',1.5);
plot(fshift,mag2db(abs(Newsig_wfreq4)),'Linewidth',1.5);
grid on
currentHandle = gca; % axis handle
set(currentHandle, 'FontSize', 17)
ylabel('Magnitude (dB)')
xlabel('Frequency (Hz)')
title('Windowing with rect function')

legend('True','Effective','75%','60%','40%','10%')

% 'True' shows the frequency response of the un-cutted filter, so just the
% butterworth filter. 'Effective' shows the frequency response of the
% filter cutted in the impulse response at the effective length. '75%'
% shows the frequency response of the filter cutted at 75% of the effective
% length ect.

%% Windowing with Blackman. Effective, 75%,60%,40% and 10% of the effective length


ha1=round(index2*0.75);
ha2=round(index2*0.60);
ha3=round(index2*0.40);
ha4=round(index2*0.10);

w=blackman(index2);
w1= blackman(ha1);
w2= blackman(ha2);
w3= blackman(ha3);
w4= blackman(ha4);

sig_w = [w.*h(1:index2); zeros((length(t)-index2),1)];
sig_w1 = [w1.*h(1:ha1); zeros((length(t)-ha1),1)];
sig_w2 = [w2.*h(1:ha2); zeros((length(t)-ha2),1)];
sig_w3 = [w3.*h(1:ha3); zeros((length(t)-ha3),1)];
sig_w4 = [w4.*h(1:ha4); zeros((length(t)-ha4),1)];

figure;
stem(t,sig_w1,'filled')
figure;
stem(t,sig_w2,'filled')
figure;
stem(t,sig_w3,'filled')
figure;
stem(t,sig_w4,'filled')
figure;
stem(t,sig_w,'filled')


sig_wfreq1=fft(sig_w1);
sig_wfreq2=fft(sig_w2);
sig_wfreq3=fft(sig_w3);
sig_wfreq4=fft(sig_w4);
sig_wfreq=fft(sig_w);
sig_wfreqTRUE=fft(h);


Newsig_wfreq1=fftshift(sig_wfreq1);
Newsig_wfreq2=fftshift(sig_wfreq2);
Newsig_wfreq3=fftshift(sig_wfreq3);
Newsig_wfreq4=fftshift(sig_wfreq4);
Newsig_wfreq=fftshift(sig_wfreq);
Newsig_wfreqTRUE=fftshift(sig_wfreqTRUE);

n=length(sig_wfreq);

fshift = (-n/2:n/2-1)*(fs/n);

figure;
hold on
plot(fshift,mag2db(abs(Newsig_wfreqTRUE)),'Linewidth',1.5);
plot(fshift,mag2db(abs(Newsig_wfreq)),'Linewidth',1.5);
plot(fshift,mag2db(abs(Newsig_wfreq1)),'Linewidth',1.5);
plot(fshift,mag2db(abs(Newsig_wfreq2)),'Linewidth',1.5);
plot(fshift,mag2db(abs(Newsig_wfreq3)),'Linewidth',1.5);
plot(fshift,mag2db(abs(Newsig_wfreq4)),'Linewidth',1.5);
grid on
currentHandle = gca; % axis handle
set(currentHandle, 'FontSize', 17)
ylabel('Magnitude (dB)')
xlabel('Frequency (Hz)')
title('Windowing with Blackmann')
currentHandle = gca; % axis handle
%set(currentHandle, 'FontSize', 17)

legend('True','Effective','75%','60%','40%','10%')















