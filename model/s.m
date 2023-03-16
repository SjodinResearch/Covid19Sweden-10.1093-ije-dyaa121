function ss = s(dm, pop_vec)
    %pop_vec=N
    %dm=d
    ss = zeros(size(dm));

    for i = 1:size(dm, 1)

        for j = 1:size(dm, 1)

            if i == j
                ss(i, j) = 0;
            else
                indx = find(dm(i, :) < dm(i, j) & dm(i, :) ~= dm(i, i));

                if (isempty(indx) == 1)
                    ss(i, j) = 0;
                    %T(i,j)= (Ni*Nj)/((Ni+s(i,j))*(Ni+Nj+s(i,j)));
                else
                    ss(i, j) = sum(pop_vec(indx'));
                    %T(i,j)= (Ni*Nj)/((Ni+s(i,j))*(Ni+Nj+s(i,j)));
                end

            end

        end

    end

end
