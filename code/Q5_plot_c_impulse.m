% Question 5
% Plot impulse.wav against the C DFT and IDFT results.

clear; close all; clc;

N = 10;
[Xwav, Fs] = audioread('impulse.wav');
h = Xwav(:, 1);
h = h(:);
L = length(h);
t = (0:L-1).' / Fs;

dft = readmatrix('impulse_dft.txt', 'FileType', 'text', 'CommentStyle', '#');
idft = readmatrix('impulse_idft.txt', 'FileType', 'text', 'CommentStyle', '#');

seg = dft(:, 1);
mag = dft(:, 5);
nSeg = max(seg) + 1;
xr = idft(:, 3);
xi = idft(:, 4);
hhat = xr(1:L);
pad = numel(xr) - L;

fprintf('Fs = %g Hz, samples = %d, segments = %d, pad = %d\n', Fs, L, nSeg, pad);
fprintf('Max |C IDFT - impulse| = %.3g\n', max(abs(hhat - h)));
fprintf('Max |imag(C IDFT)| = %.3g\n', max(abs(xi)));
fprintf('MATLAB IDFT error was 1.44e-15.\n');

figure(1);
plot(t, h, 'LineWidth', 0.8); hold on;
plot(t, hhat, '--', 'LineWidth', 0.8);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('impulse.wav and C 10-point IDFT');
legend('impulse.wav', 'C IDFT (real part)', 'Location', 'best');
saveas(gcf, 'Fig17_impulse_vs_c_idft.png');

figure(2);
plot(t, hhat - h, 'LineWidth', 0.6);
grid on;
xlabel('Time (s)');
ylabel('Error');
title('C IDFT minus impulse.wav');
saveas(gcf, 'Fig18_impulse_c_idft_error.png');

figure(3);
plot(mag, 'LineWidth', 0.5);
grid on;
xlabel('Concatenated bin index (10 bins per segment)');
ylabel('|H[k]| from C DFT');
title('C 10-point DFT magnitude of impulse.wav');
saveas(gcf, 'Fig19_impulse_c_dft_magnitude.png');

fprintf('Saved Fig17, Fig18, Fig19.\n');