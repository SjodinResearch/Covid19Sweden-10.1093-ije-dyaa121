function a = alpha_scaling(t, ti, red, pars)

    a = pars.alpha0 .* ((1 - red) ./ (1 + exp(pars.rate_of_adaptation * (t - ti))) + red);

end
