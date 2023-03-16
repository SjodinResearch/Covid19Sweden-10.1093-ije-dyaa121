function p = parameters_Swe_Corona_Radiation(omga)

    %% Time

    % p.tspan = [datenum('2020-02-24') datenum('2020-09-01')]; % original t_span
    p.tspan = [datenum('2020-02-24') datenum('2020-09-01')];

    p.ti = datenum('2020-03-29'); % 20e (inflection ~28) Time of intervention

    p.info_ti_date = datetime(p.ti, 'ConvertFrom', 'datenum');

    %% Geography

    p.Geo_Sruct_Name = 'GS'; %% test

    [p.popsize_per_region, p.travel_prob, p.SE_region_names, p.SE_healthcare_region_names, p.SE_age_matrix, p.SE_age_prop_regions] = set_geographical_variables();

    p.regions_names = p.SE_region_names;

    %% Population

    p.age_prop_kommuns = p.SE_age_prop_regions; %[SE_kommun{:,23:25}];

    %% initial conditions

    % corona_cases= 1 * 50;
    %
    % p.info__unreg_cases_in_RSthlm = corona_cases;
    %
    % p.ic_factor= corona_cases/2377000;

    p.ic_factor = ones(290, 1) * (1/100000); %150000
    municipalities = get_municipality_indices_in_region('Region Stockholm', p); % get Region Stocjholm municipality indices
    p.ic_factor(municipalities) = 1/50000; %1 / 30000;
    p.ic = (p.popsize_per_region .* p.ic_factor)'; %(S_pop *ic_factor)' ;

    X_s = (p.SE_age_prop_regions' .* p.popsize_per_region')'; %(p.age_prop_kommuns'.*p.S0)';%susceptibles

    ics = (p.age_prop_kommuns(:, 1)' .* p.ic)'; %seeding for Age group 0-59

    p.X0 = [X_s - ics, zeros(size(X_s, 1), 27)]; %% origianl

    %p.X0(:,7:9)=ics;
    p.X0(:, 7) = ics;

    p.ihcbeds = 100;

    %% parameters from input
    if isempty(omga)
        disp('omga is empty, unsing preset array.')
        p.omega_v = [0.4437 0.0155 0.2821 0.7072]; % from fit ... highbeta_3
    else
        p.omega_v = omga;
    end

    %p.c = omga(1,1);%.25; % contact-rate scaling

    %% Age group wise Contact rates

    p.omega_elements = {'cr_age_upto59', 'cr_age_upto79', 'cr_age_80plus', 'gamma'};

    p.cr_age_upto59 = p.omega_v(1, 1);
    p.cr_age_upto79 = p.omega_v(1, 2);
    p.cr_age_80plus = p.omega_v(1, 3);
    p.cr0_age_upto59 = 1;
    p.cr0_age_upto79 = 1;
    p.cr0_age_80plus = 1;

    %% under reporting and correction factors

    p.underreporting_factor = 4.6286;
    p.info_underreported_infections = [num2str(100 * (1 - (1 / p.underreporting_factor))) ' %'];

    p.icu_proportion_scaling = 1; %7/10;%5 / 4;% to achieve icu:health-care approx. 1:4 instead of 1:5;
    p.heatlhcare_proportion_scaling = 0.5; % to decrease hospitalizations slightly

    %% travel

    p.alpha0 = 0.01;
    p.alpha_reduction = 1; %0.1;

    %% Epidemiology and health care

    % Grisselli tabular structure:
    % [0-20 21-40 41-50 51-60 61-70 71-80 81-90 91-100]
    % [tot; Died in icu; Discherged from icu; Still in icu as of 3/25/2020]
    p.Grisselli_deaths_data = [2 56 142 423 596 340 21 1; ...
                                0 4 16 63 174 136 11 1; ...
                                0 20 35 90 69 40 2 0; ...
                                2 32 91 270 353 164 8 0];

    g = p.Grisselli_deaths_data;

    p.cGdd = [sum(g(2, 1:4), 2) / sum(g(1, 1:4)) sum(g(2, 5:6), 2) / sum(g(1, 5:6)) sum(g(2, 7:8), 2) / sum(g(1, 7:8))];

    p.age_class_dies_if_no_ic = [0 0 1];

    p.age_class_recovers_if_no_ic = [1 1 0];

    p.triage = [0.0 0.1 0.85];

    p.ICU_beds_Swe = 526;

    p.UK_age_class = {'0-9', '10-19', '20-29', '30-39', '40-49', '50-59', '60-69', '70-79', '80+'};

    p.UK_sympt_req_hospital = [0.1 0.3 1.2 3.2 4.9 10.2 16.6 24.3 27.3] / (p.underreporting_factor * 100); % * p.heatlhcare_proportion_scaling / (p.underreporting_factor*100);

    p.UK_hosp_req_crit_care = [5 5 5 5 6.3 12.2 27.4 43.2 70.9] * p.icu_proportion_scaling / 100;

    p.age_class = {'0-59', '60-79', '80+'};

    p.age_matrix = p.SE_age_matrix;

    p.proportion_of_I_to_health_care = average_HC_proportions(p.age_matrix, p.UK_sympt_req_hospital);

    p.proportion_of_health_care_to_intensive_health_care = average_HC_proportions(p.age_matrix, p.UK_hosp_req_crit_care);

    p.pre_hospitalization_period = 3;

    p.health_care_period = 7;

    p.ihc_period = 10; %12

    p.non_hospitalized_mortality_risk = p.cGdd * 0.0069; %0.008 * 1;

    p.ihc_proportion_dead = p.cGdd; %[0.364 0.74 0.857];% <- Based on Grisselli et al.%[0.2 0.49 0.49];

    p.recovery_period_after_ihc = 7;

    % Read Swedish Care unit and ICU data
    p.careunits_data = readtable('SE_ICU_data.csv');

    p.beta = 0.8924; %0.77;%0.7171;%0.8064;%0.91;

    p.latent_period = 4;

    p.infectious_period = p.omega_v(1, 4) * 5; % Assuming that various degrees of isolation will work to reduce this period

    p.gamma = 1 / p.infectious_period;

    % 'post_infectious_period': The number of days, 12 - p.infectious_period,
    % after being infectious, that people still carrying virus and would test
    % positive in PCR.
    p.post_infectious_period = 12 - p.infectious_period;

    p.rate_of_adaptation = 5/18; %5/16

end
