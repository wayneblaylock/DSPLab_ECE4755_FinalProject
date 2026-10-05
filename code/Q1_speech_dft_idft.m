% Question 1
% 10-point DFT and IDFT of speech.wav.
% 10 samples per segment, last segment zero-padded.

clear; close all; clc;

N = 10;
wavFile = 'speech.wav';
if ~isfile(wavFile)
    error(['Cannot find speech.wav. Current folder is %s.\n' ...
           'Copy speech.wav here, then run again.'], pwd);
end

[Xwav, Fs] = audioread(wavFile);
x = Xwav(:, 1);
x = x(:);
L = length(x);

fprintf('speech.wav: %d samples, Fs = %g Hz, duration = %.3f s\n', ...
    L, Fs, L/Fs);

nSeg = ceil(L / N);
pad = nSeg * N - L;
xp = [x; zeros(pad, 1)];

Xseg = zeros(nSeg, N);
xhat = zeros(nSeg * N, 1);

for s = 1:nSeg
    blk = xp((s-1)*N + (1:N)).';
    Xs = my_dft(blk, N);
    Xseg(s, :) = Xs;
    xhat((s-1)*N + (1:N)) = my_idft(Xs, N).';
end

xhat_trim = real(xhat(1:L));
fprintf('Segments = %d, zero-pad on last segment = %d\n', nSeg, pad);
fprintf('Max |real(IDFT) - x| = %g\n', max(abs(xhat_trim - x)));
fprintf('Max |imag(IDFT)| = %g\n', max(abs(imag(xhat))));

t = (0:L-1).' / Fs;

figure(1);
plot(t, x, 'LineWidth', 0.8);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Lab 7 speech.wav — time domain');
saveas(gcf, 'FP1_Fig1_speech_time.png');

segShow = 1;
nIdx = (segShow-1)*N + (0:N-1);
seg = zeros(N, 1);
valid = nIdx < L;
seg(valid) = x(nIdx(valid)+1);
f = (0:N-1) * Fs / N;

figure(2);
subplot(3,1,1);
stem(nIdx, seg, 'filled');
grid on;
xlabel('Sample index n');
ylabel('x[n]');
title(sprintf('Segment %d (10 samples)', segShow));

subplot(3,1,2);
stem(f, abs(Xseg(segShow, :)), 'filled');
grid on;
xlabel('Frequency (Hz)');
ylabel('|X[k]|');
title(sprintf('10-point DFT magnitude, segment %d', segShow));

subplot(3,1,3);
stem(f, angle(Xseg(segShow, :)), 'filled');
grid on;
xlabel('Frequency (Hz)');
ylabel('angle X[k] (rad)');
title(sprintf('10-point DFT phase, segment %d', segShow));
saveas(gcf, 'FP1_Fig2_one_segment_dft.png');

figure(3);
imagesc(f, 1:nSeg, abs(Xseg));
axis xy;
xlabel('Frequency (Hz)');
ylabel('Segment index');
title('|X[k]| of every 10-sample segment');
colorbar;
saveas(gcf, 'FP1_Fig3_all_segment_spectra.png');

magCat = abs(Xseg).';
magCat = magCat(:);
figure(4);
plot(0:length(magCat)-1, magCat, 'LineWidth', 0.6);
grid on;
xlabel('Concatenated bin index (10 bins per segment)');
ylabel('|X[k]|');
title('Combined 10-point DFT spectra of speech.wav');
saveas(gcf, 'FP1_Fig4_combined_dft_mag.png');

figure(5);
plot(t, x, 'LineWidth', 0.8); hold on;
plot(t, xhat_trim, '--', 'LineWidth', 0.8);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('10-point IDFT of each segment, pad samples removed');
legend('Original speech', 'IDFT (real part)', 'Location', 'best');
saveas(gcf, 'FP1_Fig5_idft_vs_original.png');

figure(6);
plot(t, xhat_trim - x, 'LineWidth', 0.6);
grid on;
xlabel('Time (s)');
ylabel('Error');
title('IDFT reconstruction error');
saveas(gcf, 'FP1_Fig6_idft_error.png');

fprintf('Saved FP1_Fig1 through FP1_Fig6 in %s\n', pwd);
fprintf('Bin spacing Fs/10 = %g Hz\n', Fs/N);

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

function x = my_idft(X, N)
    X = X(:).';
    if length(X) < N
        X = [X, zeros(1, N - length(X))];
    else
        X = X(1:N);
    end
    x = zeros(1, N);
    for n = 0:N-1
        acc = 0;
        for k = 0:N-1
            ang = 2*pi*k*n/N;
            acc = acc + X(k+1) * exp(1j*ang);
        end
        x(n+1) = acc / N;
    end
end