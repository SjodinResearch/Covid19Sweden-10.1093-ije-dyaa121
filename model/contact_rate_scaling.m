function c = contact_rate_scaling(t, ti, red, c0, pars)

    c = c0 .* ((1 - red) ./ (1 + exp(pars.rate_of_adaptation * (t - ti))) + red);

end
