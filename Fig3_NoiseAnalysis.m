% Feng et al, Ultra-wideband optical coherence elastography from acoustic to ultrasonic frequencies
% This script is used to show the phase noise correction in Fig 3 and Fig S4
% @ Guoyang Li & Xu Feng

clear; clc; close all;

%% READ DATA 
k0 = 2*pi/1280e-9; % central wavenumber
M = 108; % number of Alines per M-scan
N = 25; % cycles of aliasing 
load('Data_679.6kHz.mat');

% HARD CODE %
shift = 12; % pixel shift
center = 107; % pixel index of mainlobe
% peak1 = center + shift; % positive sideband
% peak2 = center - shift; % negative sideband

fA = 43.2e3; % Aline rate
apparent_freq = 10e3; % Hz
Na = 31;
fs = Na*(fA/2) + apparent_freq; % modulation frequency

%% EXPERIMENT 
meas = 900; % number of M-scan traces
istart = 26; 

if mod(Na,2) == 1
    Pixel1 = N + 1;
    Pixel2 = M - N + 1;
    pinc = (Na+1)*pi - fs/fA*2*pi; % theoretical phase
elseif mod(Na,2) == 0
    Pixel1 = M/2 + N + 1;
    Pixel2 = M/2 - N + 1;
    pinc = fs/fA*2*pi - Na*pi; % theoretical phase
end

D1 = D1(:,1:meas);
D1_uncorrected = D1;

%% Phase noice correction
s1 = s1(:,1:meas);
phas1 = s1;
for ii = 1:meas
%     s1(:,ii) = hilbert(s1(:,ii));
    phas1(:,ii) = unwrap(angle(s1(:,ii)));
end

phas0 = phas1(istart,:); % initial phase
phas1 = phas1(istart:M+istart-1,:);
phas1 = phas1 - phas1(1,:);
dphas1 = diff(phas1);
dphiLin = ones(M-1,meas)*pinc;
dphac = dphas1 - dphiLin; % phase jitter

Amp1 = abs(D1);
Pha1 = unwrap(angle(D1));
dPha1 = diff(Pha1) + dphac; % with correction
dPha1_uncorrected = diff(Pha1); % without correction
Pha1 = cumsum([Pha1(1,:) + phas0; dPha1]);
D1 = Amp1.*exp(i*Pha1);

FD1 = fftshift((fft(D1,[],1)))/M;
Amp1 = abs(FD1(Pixel1,:))*1e9;
Phi1 = angle(FD1(Pixel1,:));

% uncorrected
FD1_uncorrected = fftshift((fft(D1_uncorrected,[],1)))/M;
Amp1_uncorrected = abs(FD1_uncorrected(Pixel1,:))*1e9;
Phi1_uncorrected = angle(D1_uncorrected(Pixel1,:));


%% Plot
% Fig. S5, laser time jitter
figure(1); hold on;
dt = dphac/(2*pi*fs); % time jitter
dt = reshape(dt,[(M-1)*meas 1]);
dt0 = std(dt);
histfit(dt*1e9,50);
title('Laser time jitter');
xlabel('Time jitter (ns)');
fprintf('Standard deviation of laser time jitter is %4.1f ns\n',dt0*1e9);

% Reference waveform applied to the PZT
figure(2); hold on;
subplot(2,1,1)
plot(real(0.5*s1(:,1:30:900)));
xlabel('A-lines (#)'); ylabel('Voltage (V)');
title('Reference waveform from function generator')

subplot(2,1,2)
plot(phas1(:,1:30:900));
xlabel('A-lines (#)'); ylabel('Phase (rad)')

% Fig. 3c, M-scan traces before and after correction
figure(3);
subplot(2,1,1)
plot(real(D1_uncorrected(:,:)),'-')
title('Uncorrected');
xlabel('A-lines (#)');ylabel('Disp (nm)');
xlim([1 M]);

subplot(2,1,2)
plot(real(D1(:,:)),'-')
title('Corrected');
xlabel('A-lines (#)');ylabel('Dips (nm)');
xlim([1 M]);

% Fig. 3d, Vibration phase noise after correction
figure(4); 
subplot(2,1,1),hold on,
plot(Phi1_uncorrected,'r.'); plot(Phi1_uncorrected,'b-') 
xlabel('M-scan traces (#)');ylabel('Vibration phase (rad)');
title('Uncorrected')
xlim([0,900])
% ylim([0 0.5]);

subplot(2,1,2),hold on,
plot(Phi1,'r.'); plot(Phi1,'b-') 
xlabel('M-scan traces (#)');ylabel('Vibration phase (rad)');
title('Corrected')
xlim([0,900])
ylim([0 0.5]);

fprintf('Before correction: mean amplitude is %4.2f nm, phase noise is %4.4f rad\n',mean(Amp1_uncorrected), std(Phi1_uncorrected));
fprintf('After correction: mean amplitude is %4.2f nm, phase noise is %4.4f rad\n',mean(Amp1), std(Phi1));


