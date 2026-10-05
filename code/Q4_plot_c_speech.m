% Question 4
% Plot speech.wav against the C DFT and IDFT results.

clear; close all; clc;

N = 10;

[Xwav, Fs] = audioread('speech.wav');
x = Xwav(:, 1);
x = x(:);
L = length(x);
t = (0:L-1).' / Fs;

dft = readmatrix('speech_dft.txt', 'FileType', 'text', 'CommentStyle', '#');
idft = readmatrix('speech_idft.txt', 'FileType', 'text', 'CommentStyle', '#');

seg = dft(:, 1);
Xr = dft(:, 3);
Xi = dft(:, 4);
mag = dft(:, 5);
nSeg = max(seg) + 1;

xr = idft(:, 3);
xi = idft(:, 4);
xhat = xr(1:L);                 % drop the 8 pad samples
pad = numel(xr) - L;

err_c = max(abs(xhat - x));
max_imag = max(abs(xi));

fprintf('Fs = %g Hz, samples = %d, segments = %d, pad = %d\n', ...
    Fs, L, nSeg, pad);
fprintf('Max |C IDFT - speech| = %.3g\n', err_c);
fprintf('Max |imag(C IDFT)| = %.3g\n', max_imag);

% speech and C IDFT
figure(1);
plot(t, x, 'LineWidth', 0.8); hold on;
plot(t, xhat, '--', 'LineWidth', 0.8);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Original speech and C 10-point IDFT');
legend('speech.wav', 'C IDFT (real part)', 'Location', 'best');
saveas(gcf, 'Fig7_speech_vs_c_idft.png');

% IDFT error
figure(2);
plot(t, xhat - x, 'LineWidth', 0.6);
grid on;
xlabel('Time (s)');
ylabel('Error');
title('C IDFT minus original speech');
saveas(gcf, 'Fig8_c_idft_error.png');

% C DFT magnitude
figure(3);
plot(mag, 'LineWidth', 0.5);
grid on;
xlabel('Concatenated bin index (10 bins per segment)');
ylabel('|X[k]| from C DFT');
title('Question 2 C 10-point DFT magnitude');
saveas(gcf, 'Fig9_c_dft_magnitude.png');

% one segment, C vs MATLAB DFT
segShow = 4500;                 % about 1.02 s
rows = dft(seg == segShow, :);
n0 = segShow * N;
blk = x(n0 + (1:N)).';
Xmat = my_dft(blk, N);

figure(4);
stem(0:N-1, rows(:, 5), 'filled'); hold on;
stem(0:N-1, abs(Xmat), 'LineWidth', 1.2);
grid on;
xlabel('Bin k');
ylabel('|X[k]|');
title(sprintf('Segment %d DFT magnitude, C vs MATLAB', segShow));
legend('C DFT (Question 2)', 'MATLAB DFT (Question 1 method)', ...
    'Location', 'best');
saveas(gcf, 'Fig10_segment4500_c_vs_matlab.png');

fprintf('Segment %d max |C - MATLAB| = %.3g\n', ...
    segShow, max(abs(rows(:, 3) + 1j*rows(:, 4) - Xmat.')));
fprintf('Question 1 MATLAB IDFT error was 2.78e-16.\n');
fprintf('Question 3 C IDFT error is %.3g (float32 sample file).\n', err_c);

function X = my_dft(x, N)
    x = x(:).';
    if length(x) < N
        x = [x, zeros(1, N - length(x))];
    else
        x = x(1:N);
    end
    X = zeros(1, N);
    for k = 0:N-1
        acc = 0;
        for n = 0:N-1
            ang = -2*pi*k*n/N;
            acc = acc + x(n+1) * exp(1j*ang);
        end
        X(k+1) = acc;
    end
end