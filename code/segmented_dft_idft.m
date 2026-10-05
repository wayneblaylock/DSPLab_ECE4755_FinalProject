function [Xseg, xhat, info] = segmented_dft_idft(x, N)
% Split x into length-N blocks, pad the last one, DFT and IDFT each block.

    x = x(:);
    L = length(x);
    nSeg = ceil(L / N);
    pad = nSeg*N - L;
    xp = [x; zeros(pad, 1)];

    Xseg = zeros(nSeg, N);
    xhat = zeros(nSeg*N, 1);
    for s = 1:nSeg
        blk = xp((s-1)*N + (1:N)).';
        Xs = my_dft(blk, N);
        Xseg(s, :) = Xs;
        xhat((s-1)*N + (1:N)) = my_idft(Xs, N).';
    end

    info.L = L;
    info.N = N;
    info.nSeg = nSeg;
    info.pad = pad;
end
