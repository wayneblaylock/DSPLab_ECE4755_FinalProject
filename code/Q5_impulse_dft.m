% Question 5
% Same 10-point DFT/IDFT as Question 1, on impulse.wav.

clear; close all; clc;

N = 10;
wavFile = 'impulse.wav';
if ~isfile(wavFile)
    error('Cannot find impulse.wav. Current folder is %s', pwd);
end

[Xwav, Fs] = audioread(wavFile);
x = Xwav(:, 1);
x = x(:);
L = length(x);

fprintf('impulse.wav: %d samples, Fs = %g Hz, duration = %.3f s\n', ...
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
fprintf('Max |real(IDFT) - h| = %g\n', max(abs(xhat_trim - x)));
fprintf('Max |imag(IDFT)| = %g\n', max(abs(imag(xhat))));
fprintf('Bin spacing Fs/10 = %g Hz\n', Fs/N);

t = (0:L-1).' / Fs;

figure(1);
plot(t, x, 'LineWidth', 0.8);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Lab 7 impulse.wav — time domain');
saveas(gcf, 'Fig11_impulse_time.png');

% skip the silent segments at the start
energy = sum(abs(Xseg), 2);
segShow = find(energy > 0, 1, 'first');
if isempty(segShow)
    segShow = 1;
end
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
ylabel('h[n]');
title(sprintf('Impulse segment %d (10 samples)', segShow));
subplot(3,1,2);
stem(f, abs(Xseg(segShow, :)), 'filled');
grid on;
xlabel('Frequency (Hz)');
ylabel('|H[k]|');
title(sprintf('10-point DFT magnitude, segment %d', segShow));
subplot(3,1,3);
stem(f, angle(Xseg(segShow, :)), 'filled');
grid on;
xlabel('Frequency (Hz)');
ylabel('angle H[k] (rad)');
title(sprintf('10-point DFT phase, segment %d', segShow));
saveas(gcf, 'Fig12_impulse_segment_dft.png');

figure(3);
imagesc(f, 1:nSeg, abs(Xseg));
axis xy;
xlabel('Frequency (Hz)');
ylabel('Segment index');
title('|H[k]| of every 10-sample segment');
colorbar;
saveas(gcf, 'Fig13_impulse_all_spectra.png');

magCat = abs(Xseg).';
magCat = magCat(:);
figure(4);
plot(0:length(magCat)-1, magCat, 'LineWidth', 0.6);
grid on;
xlabel('Concatenated bin index (10 bins per segment)');
ylabel('|H[k]|');
title('Combined 10-point DFT spectra of impulse.wav');
saveas(gcf, 'Fig14_impulse_combined_dft.png');

figure(5);
plot(t, x, 'LineWidth', 0.8); hold on;
plot(t, xhat_trim, '--', 'LineWidth', 0.8);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('10-point IDFT of each impulse segment, pad removed');
legend('impulse.wav', 'IDFT (real part)', 'Location', 'best');
saveas(gcf, 'Fig15_impulse_idft_vs_original.png');

figure(6);
plot(t, xhat_trim - x, 'LineWidth', 0.6);
grid on;
xlabel('Time (s)');
ylabel('Error');
title('Impulse IDFT reconstruction error');
saveas(gcf, 'Fig16_impulse_idft_error.png');

fprintf('Saved Fig11 through Fig16. Segment plotted in Fig12 = %d\n', segShow);

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