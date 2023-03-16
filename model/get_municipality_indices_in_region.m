function municipalities = get_municipality_indices_in_region(region, pars)
    %GET_MUNICIPALITIES_IN_REGION Array of column numbers for municipalities in
    % a region

    if strcmp(region, 'Sweden')
        municipalities = 1:290;
    else
        municipalities = find(strcmp(region, pars.SE_healthcare_region_names));
    end

end
