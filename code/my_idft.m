function x = my_idft(X, N)
% N-point IDFT from the sum.

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
