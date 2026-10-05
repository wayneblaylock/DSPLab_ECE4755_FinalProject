function X = my_dft(x, N)
% N-point DFT from the sum. Zero-pads if x is short.

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
