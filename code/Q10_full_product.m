% Question 10
% Full-file 10-point product, then IDFT.
N = 10;
speechFile = 'speech.wav';
impulseFile = 'impulse.wav';
if ~isfile(speechFile) || ~isfile(impulseFile)
    error('Need speech.wav and impulse.wav in %s', pwd);
end

[xs, Fs] = audioread(speechFile);
[hs, FsH] = audioread(impulseFile);
x = xs(:, 1);
h = hs(:, 1);
fprintf('speech:  %d samples, Fs = %g Hz\n', length(x), Fs);
fprintf('impulse: %d samples, Fs = %g Hz\n', length(h), FsH);

% shorter file is zero-padded so the segments line up
L = max(length(x), length(h));
x = [x; zeros(L - length(x), 1)];
h = [h; zeros(L - length(h), 1)];
nSeg = ceil(L / N);
pad = nSeg * N - L;
x = [x; zeros(pad, 1)];
h = [h; zeros(pad, 1)];

mag_acc = zeros(N, 1);
y_time = zeros(nSeg * N, 1);
for s = 1:nSeg
    xb = x((s-1)*N + (1:N));
    hb = h((s-1)*N + (1:N));
    X = dft10(xb);
    H = dft10(hb);
    Y = X .* H;
    mag_acc = mag_acc + abs(Y);
    y_time((s-1)*N + (1:N)) = real(idft10(Y));
end
mag_avg = mag_acc / nSeg;
y_time = y_time(1:L);

figure;
stem(0:N-1, mag_avg, 'filled');
grid on;
xlabel('DFT bin k');
ylabel('Average |X[k] H[k]|');
title('Question 10: full-file 10-point product magnitude');

figure;
t = (0:length(y_time)-1).' / Fs;
plot(t, y_time, 'LineWidth', 0.6);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Question 10: IDFT of X[k] H[k], full files');

% same excerpts as Q9, for the comparison plot
speechStart = 43968;
impulseStart = 55680;
nEx = 256;
if speechStart + nEx <= length(xs) && impulseStart + nEx <= length(hs)
    xsEx = xs(speechStart + (1:nEx), 1);
    hsEx = hs(impulseStart + (1:nEx), 1);
    mag_ex = zeros(N, 1);
    y_ex = zeros(nEx, 1);
    for s = 1:(nEx / N)
        X = dft10(xsEx((s-1)*N + (1:N)));
        H = dft10(hsEx((s-1)*N + (1:N)));
        Y = X .* H;
        mag_ex = mag_ex + abs(Y);
        y_ex((s-1)*N + (1:N)) = real(idft10(Y));
    end
    mag_ex = mag_ex / (nEx / N);

    figure;
    subplot(2, 1, 1);
    stem(0:N-1, mag_ex, 'filled');
    grid on;
    xlabel('DFT bin k');
    ylabel('|X[k] H[k]|');
    title('Question 10 check: same excerpts as Question 9');
    subplot(2, 1, 2);
    plot(0:nEx-1, y_ex, 'LineWidth', 1);
    grid on;
    xlabel('Sample within excerpt');
    ylabel('Amplitude');
    title('IDFT of the excerpt product (should match fig21)');
end

fprintf('Segments = %d, zero-pad = %d\n', nSeg, pad);
fprintf('Save the figures yourself. This script does not write files.\n');

function X = dft10(x)
    N = 10;
    n = (0:N-1).';
    X = zeros(N, 1);
    for k = 0:N-1
        ang = -2*pi*k*n/N;
        X(k+1) = sum(x(:) .* exp(1j * ang));
    end
end

function x = idft10(X)
    N = 10;
    k = (0:N-1).';
    x = zeros(N, 1);
    for n = 0:N-1
        ang = 2*pi*k*n/N;
        x(n+1) = sum(X(:) .* exp(1j * ang)) / N;
    end
end