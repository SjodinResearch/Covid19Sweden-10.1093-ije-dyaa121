function solve_Swedish_outbreak(runs, save_str, set_pars)

    addpath('./model')

    %% OMEGA MATRIX age group wise contact rate
    %omgea(:,1)= cr_upto59, omgea(:,2)= cr_upto79; omgea(:,3)= cr_80plus;
    %omgea(:,4)=gamma

    omega = [1 1 1 1;
        0.75 0.5 0.5 1;
        0.75 0.25 0.25 1;
        0.4437 0.0155 0.2821 0.7072;
        0.5 0.1 0.1 0.6];

    for i = runs%1:omega_no
        disp(['Solving ' num2str(i) ' ...'])

        file_name0 = save_str; %'Omega_vec_i3_h3_';

        file_name = horzcat(file_name0, '_', num2str(i));

        p = parameters_Swe_Corona_Radiation(omega(i, :));
        p.omega_v = omega(i, :);

        if ~isempty(set_pars)
            fin = fieldnames(set_pars);

            for j = 1:size(fin, 1)
                p.(fin{j}) = set_pars.(fin{j});
            end

        end

        [t, X] = solve_SEIR_HC_radiation(p, file_name);

    end

end
