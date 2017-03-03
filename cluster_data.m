% -----------------------
% Copyright (c) 2017, Korbinian Breinl
% All rights reserved.
% -----------------------
function ic = data_cluster( data, ic, cluster, period_array, period_array_all, period )
% DATA_CLUSTER
%
% ic = data_cluster(data, ic, cluster, period_array, period_array_all, period)
%
% Cluster the data.
%

% Filter all occurrence vectors of a period
data_period = data(period_array == period, :);

% Cluster amount vectors with normalized amounts via zscores
clustera = cluster;

if cluster - 1 > size(data_period, 1)
    clustera = size(data_period, 1);
end

data_clus = kmeans(data_period>0, clustera - 1, 'distance', 'hamming', ...
    'MaxIter', 100, 'emptyaction', 'singleton', 'replicates', 1);

% Write clusters into observation matrix with condition > 1
% (non-zero amount vectors and month) and adds 1 to cluster
% to keep zero-amount vectors of number 1
ic(ic > 1 & period_array_all == period) = data_clus + 1;

end