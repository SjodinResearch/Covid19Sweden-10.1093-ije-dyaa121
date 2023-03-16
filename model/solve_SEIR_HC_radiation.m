function [t X var_Z] = solve_SEIR_HC_radiation(pars, str_save)

    p = pars; %parameters_Swe_Corona;

    ic_X0 = reshape(p.X0', [1 numel(p.X0)]); %initial condition

    [t, X] = ode45(@(t, X) SEIR_HC_Radiation(t, X, p), p.tspan, ic_X0, odeset('NonNegative', 1:length(ic_X0)));

    tvec = p.tspan(1):p.tspan(end);
    Sol_mat = interp1(t, X, tvec);

    T = datetime(t, 'ConvertFrom', 'datenum', 'Format', 'dd.MM.yyyy');

    Tz = datetime(tvec, 'ConvertFrom', 'datenum', 'Format', 'dd.MM.yyyy');

    var_H.t = t;
    var_H.T = T;
    var_H.p = p;
    var_Z.t = tvec; var_Z.T = Tz; var_Z.p = p;
    v_name = {'S', 'E', 'I', 'J', 'V', 'IV', 'VR', 'M', 'R', 'RI'};
    model_vars = size(p.X0, 2);
    counter = 0:3:model_vars;

    for kk = 1:length(v_name)

        for jj = 1:3%24
            var_H.(horzcat(v_name{kk}, num2str(jj))) = X(:, counter(kk) + jj:model_vars:end);
            var_Z.(horzcat(v_name{kk}, num2str(jj))) = Sol_mat(:, counter(kk) + jj:model_vars:end);

        end

    end

    %% Saving the results.
    if ~isempty(str_save)
        save(str_save, 'var_H', 'var_Z');
    end

end
