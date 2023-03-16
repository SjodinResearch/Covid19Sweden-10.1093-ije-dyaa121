function prod_mat = average_HC_proportions(age_dist_mat_swe, vec_uk)

    indx = {1:2, 3:4, 5:6, 7:8, 9:10, 11:12, 13:14, 15:16, 17};
    id_uk = {1:6, 7:8, 9};
    age_prop_dist = zeros(size(age_dist_mat_swe, 1), length(indx));

    for ii = 1:length(indx)
        age_prop_dist(:, ii) = sum(age_dist_mat_swe(:, indx{ii}), 2);
    end

    age_prop_dist = age_prop_dist ./ age_dist_mat_swe(:, end);

    prod_mat = [vec_prod(age_prop_dist, vec_uk, id_uk{1}), ...
                vec_prod(age_prop_dist, vec_uk, id_uk{2}), ...
                vec_prod(age_prop_dist, vec_uk, id_uk{3})];

end

function prd = vec_prod(va, vb, id)

    prd = sum((va(:, id) ./ sum(va(:, id), 2)) .* vb(:, id), 2);

end
