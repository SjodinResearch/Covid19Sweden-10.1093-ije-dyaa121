function t = radiation_prob_mat(N, s)%radiation_model_reallyvectorized(X,N,s)

    % N is population matrix
    % X is total number of individuals that potentially could travel

    [ni, nj] = ndgrid(N, N);

    t = (ni .* nj) ./ ((ni + s) .* (ni + nj + s));
    t(1:numel(N) + 1:end) = 0;

end
