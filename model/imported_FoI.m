function IFOI = imported_FoI(travel_probability, FOI) %FR = radiation_model_import_FoI(F,N,s)

% Returns import-'force of infection'/beta

% The returned values must accordingly be multiplied by exiting suceptible
% individuals

% X is total number of individuals that potentially could travel

% N is population matrix

% [ni ,nj] = ndgrid(N, N);
% 
% T = (ni.*nj)./((ni+s).*(ni+nj+s));
% T(1:numel(N)+1:end)=0;

IFOI = (travel_probability' * FOI')';


end