function dX = SEIR_HC_Radiation(t, X, pars)

    %%% Parametrization

    %% Variables

    VAR = X2SEIR_Radiation(X);

    S = VAR.S';
    E = VAR.E';
    I = VAR.I';
    J = VAR.J';
    V = VAR.V';
    IV = VAR.IV';
    VR = VAR.VR';
    M = VAR.M';
    R = VAR.R';
    RI = VAR.RI';

    %% epidemiological parameters
    p = pars;
    gamma = p.gamma;

    epsilon = p.proportion_of_I_to_health_care;

    p_J = p.pre_hospitalization_period; %1 ./ p.health_care_rate;;

    p_H = p.health_care_period;

    nu = p.non_hospitalized_mortality_risk;

    chi = p.proportion_of_health_care_to_intensive_health_care;

    mu = p.ihc_proportion_dead;

    p_C = p.ihc_period; % 1 ./ p.intensive_health_care_recovery_rate;

    p_Htilde = p.recovery_period_after_ihc; % 1 ./ p.recovered_ihc_patient_recovery_rate;

    p_E = p.latent_period;

    p_RI = p.post_infectious_period;

    %% contact rate

    p.cr_scaling59 = contact_rate_scaling(t, p.ti, p.cr_age_upto59, p.cr0_age_upto59, p);
    p.cr_scaling79 = contact_rate_scaling(t, p.ti, p.cr_age_upto79, p.cr0_age_upto79, p);
    p.cr_scaling80p = contact_rate_scaling(t, p.ti, p.cr_age_80plus, p.cr0_age_80plus, p);
    %cr_mat = repmat([p.cr_scaling59, p.cr_scaling79, p.cr_scaling80p], length(p.popsize_per_region),1);

    cr_vec = [p.cr_scaling59, p.cr_scaling79, p.cr_scaling80p];
    cr_mat = sqrt(repelem(cr_vec', 1, 3) .* repelem(cr_vec, 3, 1));

    %% N- Matrix

    dmy = reshape(X, [size(p.X0, 2), length(X) / size(p.X0, 2)]); %reshape to find sum
    N_city = sum(dmy, 1);

    N = N_city - sum(V) - sum(IV) - sum(VR) - sum(M);

    % N_city = [sum(dmy(1:3:end,1:end)) ; sum(dmy(2:3:end,1:end)); sum(dmy(3:3:end,1:end))];
    %
    % N_tilde = N_city - V - IV - VR - M;
    %
    % N= sum(cr_vec' .* N_tilde,1);

    %% Travel

    alpha = alpha_scaling(t, p.ti, p.alpha_reduction, p); %contact_rate_scaling(t, p.ti, p.radiation_scaling);

    imports_I = imports_of_X(p.travel_prob, (I + J), alpha); %rad_T_I(1,:);% number of visiting infected individuals per day

    exports_I = exports_of_X(p.travel_prob, (I + J), alpha); %rad_T_I(2,:);% number of exiting infected individuals per day

    exports_S = exports_of_X(p.travel_prob, S, alpha); % rad_T_S(2,:);% number of exiting susceptible individuals per day

    G = cr_mat * (I + J + imports_I - exports_I);

    F = imported_FoI(p.travel_prob, G ./ N); %  FoI = (3 x 290) matrix

    %% SEIR

    % Rates:

    S2E = p.beta .* ((S - exports_S) .* (1 ./ N) .* G + exports_S .* F);

    E2I = (1 - epsilon') .* E / p_E;

    E2J = epsilon' .* E / p_E;

    I2R = I .* gamma .* (1 - nu');

    I2M = I .* gamma .* nu';

    J2V = J .* (1 - chi') ./ p_J;

    J2IV = J .* chi' .* (1 - p.triage') ./ p_J;

    J2M = J .* chi' .* p.triage' ./ p_J;

    V2R = V .* (1 / p_H);

    IV2VR = IV .* (1 - mu') ./ p_C;

    IV2M = IV .* mu' ./ p_C;

    VR2R = VR ./ p_Htilde;

    RI2R = RI ./ p_RI;

    % Equations:

    dX(1:3, :) =- S2E; %% S

    dX(4:6, :) = S2E - E2I - E2J; %% E

    dX(7:9, :) = E2I - I2R - I2M; %- I2V;%% I

    dX(10:12, :) = E2J - J2V - J2IV - J2M; %% J

    dX(13:15, :) = J2V - V2R; % I2V - V2IV - V2M%% V (H)

    dX(16:18, :) = J2IV - IV2VR - IV2M; %V2IV%% IV (C)

    dX(19:21, :) = IV2VR - VR2R; %% VR (\tilde())

    dX(22:24, :) = J2M + IV2M + I2M; %V2M%% M

    dX(25:27, :) = RI2R + V2R + VR2R; %% R (Totally recovered)

    dX(28:30, :) = I2R - RI2R; %% RI ( Recovered but still carrying virus, or in isolation, and not trasmitting virus. )

    dX = reshape(dX, size(X));

end
