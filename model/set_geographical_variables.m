function [popsize_per_region, travel_prob, SE_region_names, SE_healthcare_region_names, SE_age_matrix, SE_age_prop_regions] = set_geographical_variables()

    %% Read and sort geo data
    SE_regions = readtable('SE_AgeData_Kommuns.csv');
    SE_regions = SE_regions(find(SE_regions.Pop_sum > 2000), :);
    SE_region_names = SE_regions.Name;
    SE_healthcare_region_names = SE_regions.HealthcareRegionName;
    SE_healthcare_region = SE_regions.HealthcareRegionName;

    %% set contact structure
    popsize_per_region = SE_regions.Pop_sum;
    distance_mat = distance_matrix(SE_regions); %distance matrix
    s_matrix = s(distance_mat, popsize_per_region); %radiation_model matrix
    travel_prob = radiation_prob_mat(popsize_per_region', s_matrix);

    %% set and derive age structures per region

    SE_age_matrix = horzcat(...
        SE_regions.Ald0_5, ...
        SE_regions.Ald6_9, ...
        SE_regions.Ald10_15, ...
        SE_regions.Ald16_19, ...
        SE_regions.Ald20_24, ...
        SE_regions.Ald25_29, ...
        SE_regions.Ald30_34, ...
        SE_regions.Ald35_39, ...
        SE_regions.Ald40_44, ...
        SE_regions.Ald45_49, ...
        SE_regions.Ald50_54, ...
        SE_regions.Ald55_59, ...
        SE_regions.Ald60_64, ...
        SE_regions.Ald65_69, ...
        SE_regions.Ald70_74, ...
        SE_regions.Ald75_79, ...
        SE_regions.Ald80_w);

    SE_age_prop_regions = horzcat(SE_regions.Age_gp1, SE_regions.Age_gp2, SE_regions.Age_gp3);

end
