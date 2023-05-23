% Feng et al, Ultra-wideband optical coherence elastography from acoustic to ultrasonic frequencies
% This script is used to generate the simulation data for Fig. 2
% @ Guoyang Li & Xu Feng

clear all; clc; close all; fclose all;
% rng('defmult');

filnam = 'PLT_threescatters';

%% PARAMETER DEFINITION
% parameters for OCT signal
SNR = 40; % SNR
[r, fA, z0, dz] = deal(0.8, 43.2e3, 3e-3, 2e-3); % reflectivity, Aline rate, and coordinates of the scatter
npts = 2; % add random scatters

% sweep-source laser parameters
[P, t] = laser(fA, 0.5); % laser power, alpha
[lamda0, dl] = deal(1280e-9, 80e-9); % wavelength tuning

% parameters  for vibration signal
apparent_freq = 2.4e3; % define an apparent freq between 0 ~ 0.5*Aline rate, Hz
fm = 16*fA - apparent_freq; % modulation frequency
M = 216; % number of Alines per M-scan trace
N = round(M/fA*apparent_freq); % cycles of aliasing

% define peak index
if mod(N,2) == 0
    Pixel1 = round(M/2 + N + 1);
    Pixel2 = round(M/2 - N + 1);
elseif mod(N,2) == 1
    Pixel1 = N + 1;
    Pixel2 = M - N + 1;
end

% Main Lobe signal
ind = 277;

% Sidelobe signal
%%% for fm = 688.8 kHz (used in this paper)
dind = 16 % peak shift, determined by fm
zpts = 0.088*[-1 1]*dz + z0; % 0.088 to make sure 3 scatterers coincide

%%% for fm = 780 kHz (another example)
% % fm = 18*fA - apparent_freq;
% % dind = 18 % peak shift, determined by fm
% % zpts = 0.099*[-1 1]*dz + z0; % 0.099 to make sure 3 scatterers coincide

T = 1/fA;
k0 = 2*pi/lamda0; k1 = (2*pi/lamda0 - 2*pi/(lamda0+dl))/T;
k = k0 + k1*t;
A = length(t);

rpts = exp(-(zpts - z0).^2/dz^2*128); % reflectivity

%% Three scatterers
% scatterer 2
z = z0;

deltaz = 10e-9; % amplitude, unit: m
phi0 = pi/6; % vibration phase

for m = 0:M - 1
    
    I(:,m + 1) = r*P.*cos(2*k0*z + 2*k1*t*z + 2*k0*deltaz*sin(2*pi*fm*t + 2*pi*fm*m*T + phi0)...
        + 2*k1*deltaz*t.*sin(2*pi*fm*t + 2*pi*fm*m*T + phi0));
    
    I(:,m + 1) = awgn(I(:,m + 1),SNR);
    F0(:,m+1) = fft(I(:,m+1));
    
end

% scatterer 1
z = zpts(1); r = rpts(1);
deltaz = 10e-9; % amplitude, unit: m
for m = 0:M - 1
    
    I(:,m + 1) = r*P.*cos(2*k0*z + 2*k1*t*z + 2*k0*deltaz*sin(2*pi*fm*t + 2*pi*fm*m*T + phi0)...
        + 2*k1*deltaz*t.*sin(2*pi*fm*t + 2*pi*fm*m*T + phi0));
    
    I(:,m + 1) = awgn(I(:,m + 1),SNR);
    F1(:,m+1) = fft(I(:,m+1));
    
end

% scatterer 3
z = zpts(2); r = rpts(2);
deltaz = 30e-9; % amplitude, unit: m
for m = 0:M - 1
    
    I(:,m + 1) = r*P.*cos(2*k0*z + 2*k1*t*z + 2*k0*deltaz*sin(2*pi*fm*t + 2*pi*fm*m*T + phi0)...
        + 2*k1*deltaz*t.*sin(2*pi*fm*t + 2*pi*fm*m*T + phi0));
    
    I(:,m + 1) = awgn(I(:,m + 1),SNR);
    F2(:,m+1) = fft(I(:,m+1));
    
end

% combined signal
F = F0 + F1 + F2;

%% Fig. 2a
figure(1)
subplot(4,1,1) % carrier signal
plot(10*log10(F0(1:1024,1).*conj(F0(1:1024,1))));xlim([100 500]); ylim([-20 55])
legend('Scatterer 1')

subplot(4,1,2) % negative sideband
plot(10*log10(F2(1:1024,1).*conj(F2(1:1024,1))));xlim([100 500]); ylim([-20 55])
legend('Scatterer 2')

subplot(4,1,3) % positive sideband
plot(10*log10(F1(1:1024,1).*conj(F1(1:1024,1))));xlim([100 500]); ylim([-20 55])
legend('Scatterer 3')

subplot(4,1,4) % all combined
plot(10*log10(F(1:1024,1).*conj(F(1:1024,1))));xlim([100 500]); ylim([-20 55])
legend('Combined')

F0 = F(ind - dind,:)';
F1 = F(ind,:)';
F0_bgd = mean(F0); % static signal

%% Fig. 2b
figure(2);
subplot(5,1,1); hold on;
plot(real(F1),'k-'); plot(imag(F1),'k--');
title('Original data');
xlabel('M-scan'); xlim([1 216]);
F1_bgd = mean(F1);ylim([-500 300]);

subplot(5,1,2); hold on;
plot(repmat(real(F0_bgd),1,M),'k-'); plot(repmat(imag(F0_bgd),1,M),'k--');
title('Background static signal');
xlabel('M-scan'); xlim([1 216]); ylim([-200 100]);

F1 = F1 - mean(F1);  % subtract the static signal
subplot(5,1,3); hold on;
plot(real(F1),'k-'); plot(imag(F1),'k--');
title('After the static signal is subtracted');
xlabel('M-scan'); xlim([1 216]); ylim([-50 50]);

F1 = F1/mean(F0)/k0*1e9;  % subtract static phase
subplot(5,1,4); hold on;
plot(real(F1),'k-'); plot(imag(F1),'k--');
% title('Substract optical phase at z_0');
title('After normalization with the background signal')
xlabel('M-scan'); xlim([1 216]);
ylabel('Amp. (nm)'); ylim([-50 50]);

FF1 = fftshift(fft(F1))/M; % FFT
hz = linspace(-fA/2,fA/2,M)/1e3; % kHz
subplot(5,1,5); hold on;
plot(hz,abs(FF1),'k-');
plot(hz(Pixel1),abs(FF1(Pixel1)),'ro');
plot(hz(Pixel2),abs(FF1(Pixel2)),'bo')
title('Fourier domain');
xlabel('Frequency');
ylabel('Amp. (nm)');
ylim([0 40]);

% save the entire workspace
% save([filnam, '.mat']);

%%
% ---/--- functions ---/--- %
function [P, t] = laser(fA, sigma)
% This function is used to generate sweep-source laser, the output is the
% laser power and time of each scan.
% Hardcode: sampling frequency of the photodector.

fe = 2048*fA; % Sampling frequency of the photodector

T = 1/fA;
dt = 1/fe;
t = [-T/2:dt:T/2]';

P = exp(-4*log(2)*t.^2/(sigma*T)^2); % Laser power

end

