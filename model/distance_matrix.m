function d = distance_matrix(cities_table)

    % Returns "geodesic" distance in km from lat long ccordinates

    %% Distance using matrices - efficient calculation
    ln_mat = [cities_table.lat, cities_table.long];
    n = length(ln_mat);

    DM = [repelem(ln_mat, n, 1) repmat(ln_mat, n, 1)];

    dm = deg2km(distance(DM(:, 1), DM(:, 2), DM(:, 3), DM(:, 4)));

    d = reshape(dm, [n n]); % reshape

    % %%%%%%%%%%%%%%%%%
    % c = cities_table;
    % d = zeros(size(c,1));
    % for i = 1:size(c,1)
    %     for j = 1:size(c,1)
    %         lati = c(i,:).lat;
    %         loni = c(i,:).long;
    %         latj = c(j,:).lat;
    %         lonj = c(j,:).long;
    %         d(i,j) = distance(lati, loni, latj, lonj);
    %     end
    % end
    % d = deg2km(d); % convert to km
end
