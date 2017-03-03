%% Multisite Markov Model (TripleM)

% Copyright (c) 2017, Korbinian Breinl
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE

% The software TripleM uses the subroutine trans_count.m 
% Copyright (c) 2014, Shawn Pethel

%%

clear


%% Read the data
% Import the observed precipitation data and model parameters

disp('Reading the input data...')

% Read in dates (OBS_DT), observed precipitation data (OBS),
% and number of sites (nos)
[obs_dt, obs, NOS] = read_file('rain.csv');

% Read the necessary parameters for simulation
P_THR = read_param('p_thr');
STRT_SIM = read_param('start_sim');
LENGTH = read_param('length');
MAX_DUP = read_param('max_dup');
MC_ORDER = read_param('mc_order');
CL_PERIOD = read_param('cl_period');
PARAM_P = read_param('param_p');
if isequal(PARAM_P, 'on')
    MIN_SAMPLE = read_param('min_sample');
    CORR_RAND = read_param('corr_rand');
    DIST = read_param('dist');
    if isequal(DIST, 'Weibull-GPD')
        PARETO_THR = read_param('pareto_thr');
    end
end

%% Initialize the model
% Perform basic operations to initialize the model

disp('Initializing the model...')

% Transform all values below the predefined
% minimum precipitation threshold (P_THR) to 0
obs(obs <= P_THR) = 0;

% Identify the unique precipitation occurrence vectors
% for the observation period (cluster vector)
[~, ~, cv_obs] = unique(obs > 0, 'rows');

% Get a synthetic date vector for the simulation period
% following the STRT_SIM and LENGTH parameters
sim_dt = (datetime(STRT_SIM, 1, 1) : datetime(STRT_SIM + LENGTH - 1, 12, 31))';

% Derive the corresponding month for the observation and simulation periods
[~, obs_mm, ~] = ymd(obs_dt);
[~, sim_mm, ~] = ymd(sim_dt);

% Establish the corresponding season for each month
% (1: winter, 2: spring, 3: summer, 4: autumn)
month2season = [1; 1; 2; 2; 2; 3; 3; 3; 4; 4; 4; 1];

% Derive the corresponding time (obs_vec, sim_vec) and transition
% (obs_tns, sim_tns) vectors for the observation and simulation
% periods according to the CL_PERIOD parameter. 
if isequal(CL_PERIOD, 'month')
	obs_vec = obs_mm;
	sim_vec = sim_mm;
    obs_tns = obs_mm;
	sim_tns = sim_mm;
elseif isequal(CL_PERIOD, 'season')
	obs_vec = month2season(obs_mm);
	sim_vec = month2season(sim_mm);
    obs_tns = month2season(obs_mm);
	sim_tns = month2season(sim_mm);
else
    error('Error: the value provided for the CL_PERIOD parameter is not valid.')
end


%% Cluster the observed precipitation amount vectors

% Ensure clustering to a minimum number of clusters if user selects 0%
% for clustering (otherwise no clustering is applied)
if MAX_DUP == 0
    MAX_DUP = realmin;
end

% Cluster only if there is more than one site
if NOS > 1
    disp(['Clustering the precipitation occurrence vectors with maximum duplication rate of ' ...
          num2str(round(MAX_DUP, 1)) '%...'])
    
    % Filter all unique amount vectors > 1
    % ("ic = 1" are vectors with zero precipitation at all sites)
    obs_data = obs(cv_obs > 1, :);
    data_dt = obs_dt(cv_obs > 1);
    [~, data_mm, ~] = ymd(data_dt);
    if isequal(CL_PERIOD, 'season') 
        data_sn = month2season(data_mm);
    end

    % Initialize the variables for the clustering process
    uniq = 0;
    cluster = 3;

    while uniq < MAX_DUP / 100
        % Apply clustering
        if isequal(CL_PERIOD, 'month')
            for month = 1:12
                cv_obs = cluster_data(obs_data, cv_obs, cluster, data_mm, obs_vec, month);
            end
        else
            for season = 1:4
                cv_obs = cluster_data(obs_data, cv_obs, cluster, data_sn, obs_vec, season);
            end
        end
     
        % Check the % of duplicated vectors by fitting different
        % transition matrices according to the value of cluster in the loop
        [f, ~, ~, ~] = trans_count(cv_obs, MC_ORDER);
        uniq = sum(nonzeros(f) == 1) / sum(nonzeros(f));
        cluster = cluster + 1;
     
    end
 
elseif NOS == 1
    % Give the user the % of duplicated days when not using clustering
    % of precipitation patterns (FSS version)
    [f, ~, ~, ~] = trans_count(cv_obs, MC_ORDER);
    disp(['Duplication rate is ' ...
          num2str(round(sum(nonzeros(f) == 1) / sum(nonzeros(f)) * 100), 2) ...
          '% on average...'])
   
else
    error('Error: At least one data site needs to be included!')
 
end


%% Generate Markov time series

disp('Simulating Markov process...')

% Simulate an arbitraty January starting point
occv = cv_obs(obs_mm == 1);  % occurrence vector
occv(length(occv)+1 : length(occv)+MC_ORDER) = occv(1:MC_ORDER);
ran = randi([1, length(occv)-MC_ORDER-1]);

% Initialize the cluster vector for the simulation period
cv_sim = zeros(size(sim_dt));
cv_sim(1 : MC_ORDER) = occv(ran : ran+MC_ORDER-1);

% Initialize the Markov chain
[f, w, ~, ~] = trans_count(cv_obs(obs_vec == 1), MC_ORDER);
i = MC_ORDER + 1;

while i < length(cv_sim) + 1
 
    % Derive new transition probabilities either between months or seasons
    if sim_tns(i) ~= sim_tns(i - 1)
        [f, w, ~, ~] = trans_count(cv_obs(obs_tns == sim_tns(i)), MC_ORDER);
    end

    past_state = cv_sim(i-MC_ORDER : i-1);
    k = find(ismember(w, past_state', 'rows'));

    % Draw a new cluster ID if transition has not been observed or if
    % beginning of month or season
    if isempty(k) == 1 || sim_vec(i) ~= sim_vec(i - 1)
        occv = cv_obs(obs_tns == sim_tns(i));
        ran = randi([1, length(occv) - MC_ORDER - 1]);
        cv_sim(i : i+MC_ORDER-1) = occv(ran : ran+MC_ORDER-1);
        i = i + MC_ORDER - 1;

    else
        cs = full(f(k, :));
        
        % Catch rare error of empty transition probabilities when last day
        % of month/season is drawn (so no follower exists)
        if sum(isnan(cs)) > 0
            cs(isnan(cs)) = 0;
            cs(randi(size(cs, 2), 1)) = 1;
        end

        % Random sampling according to frequencies of clusters
        cs = cumsum(cs / sum(cs))';
        ran = rand;
        t = 1;
        while ran > cs(t, 1)
            t = t + 1;
        end
        cv_sim(i) = w(t, MC_ORDER);

    end
    
    i = i + 1;

end


%% Resample the observed precipitation amount vectors

disp('Resampling of observed precipitation values...'),

% Initialize the simulation matrix
sim = zeros(length(sim_dt), NOS);

% Resamples either from the month or season pool depending on model type
for i = 1:length(cv_sim)
    
    if cv_sim(i) > 1
        
        obs_smp = obs(cv_obs == cv_sim(i) & obs_tns == sim_tns(i), :);
        sim(i, :)=obs_smp(randi(size(obs_smp,1)), :);
        cv_sim(i)=0;
        
    elseif cv_sim(i) == 1
        
        sim(i, :) = 0;
        
    end
end


%% Perform a parametric sampling of precipitation amounts

if isequal(PARAM_P, 'on')
 
    disp('Performing parametric sampling of the precipitation values...')
    
    % Save resampled simulations for reshuffling in the bottom
    sim_rs = sim;

    % Generate random numbers for parametric sampling
    ran = rand(size(sim, 1), NOS);
    
    % Pre-allocate arrays for distribution parameters
    if isequal(CL_PERIOD, 'season')
        param_wbl = zeros(NOS, 2, 4) + 9999;
        param_gam = zeros(NOS, 2, 4) + 9999;
        param_gpd = zeros(NOS, 2, 4) + 9999;
        fit_thr = zeros(4, NOS) + 9999;
        periods = 4;
    else
        param_wbl = zeros(NOS, 2, 12) + 9999;
        param_gam = zeros(NOS, 2, 12) + 9999;
        param_gpd = zeros(NOS, 2, 12) + 9999;
        fit_thr = zeros(12, NOS) + 9999;
        periods = 12;
    end

    % Compute a Cholesky factorization if correlated random numbers
    % are specified and if there is more than one site
    if isequal(CORR_RAND, 'on') && NOS > 1

        ran = randn(size(sim, 1), NOS);

        % Write all observations with at least two sites of simultaneous
        % precipitation into a new matrix
        obs_tmp = obs((sum(obs' > 0)) > 1, :);
        obs_tns_tmp = obs_tns((sum(obs' > 0)) > 1, :);

        % Initialize Cholesky matrices for each period (season or month) and
        % preallocate arrays for distribution parameters
        if isequal(CL_PERIOD, 'season')
            ch_mat = zeros(NOS, NOS, 4);
        else
            ch_mat = zeros(NOS, NOS, 12);
        end

        % Generate correlated random numbers per season/month and in case of
        % positive indefinite matrices try to derive from all observations
        for period = 1:periods
            try
                ch_mat(:, :, period) = chol(corr(obs_tmp(obs_tns_tmp == period, :)));
                ran(sim_tns == period, :) = normcdf(ran(sim_tns == period, :) * ...
                                               ch_mat(:, :, period));
            catch
                ch_mat(:, :, period) = chol(corr(obs_tmp));
                ran(sim_tns == period, :) = normcdf(ran(sim_tns == period, :) * ...
                                               ch_mat(:, :, period));
            end
        end
    
    end
    
    % Perform the parametric sampling
    for period = 1:periods
        
        % Initialize temporary arrays for the simulation period
        sim_per = sim(sim_tns == period, :);
        sim_per(:, NOS + 1) = 1:size(sim_per, 1);  % extra col for sorting
        
        sim_rs_per = sim_rs(sim_tns == period, :);
        sim_rs_per(:, NOS + 1) = 1:size(sim_rs_per, 1);
        
        for site=1:NOS
            % Initialize temporary arrays for the observation period and
            % site
            obs_site = obs(obs_tns == period, site);  % observation data (period)
            obs_site = obs_site(obs_site > 0);

            % Derive distribution parameters if the sample size is large enough
            if size(obs_site, 1) >= MIN_SAMPLE
                % Weibull distributions
                if isequal(DIST, 'Weibull') || isequal(DIST, 'Weibull-GPD')
                    a = wblfit(obs_site);
                    param_wbl(site, 1, period) = a(1);
                    param_wbl(site, 2, period) = a(2);
                    if isequal(DIST, 'Weibull-GPD') && ...
                       size(obs_site(obs_site > PARETO_THR, :), 1) >= MIN_SAMPLE
                        obs_site = obs_site(obs_site > PARETO_THR, :);
                        b = gpfit(obs_site(obs_site > 0, :) - PARETO_THR);
                        param_gpd(site, 1, period) = b(1);
                        param_gpd(site, 2, period) = b(2);
                        fit_thr(period, site) = wblcdf(PARETO_THR, a(1), a(2));
                    end
                % Gamma distribution
                elseif isequal(DIST, 'Gamma')
                    a = gamfit(obs_site);
                    param_gam(site, 1, period) = a(1);
                    param_gam(site, 2, period) = a(2);
                end
            else
                if isequal(CL_PERIOD, 'season')
                    disp(['Warning: not enough samples in season ' ...
                          num2str(period) ' at site ' num2str(site) ... 
                          '. Pure resampling applied'])
                else
                    disp(['Warning: not enough samples in month ' ...
                          num2str(period) ' at site ' num2str(site) ...
                          '. Pure resampling applied'])
                end
            end

            % Sample parametric precipitation amounts
            k = sim_tns == period & sim(:, site) > 0;
            % Weibull distributions
            if isequal(DIST, 'Weibull') || isequal(DIST, 'Weibull-GPD')
                if (param_wbl(site, 1, period) ~= 9999)
                    sim(k, site) = wblinv(ran(k, site), ...
                                          param_wbl(site, 1, period), ...
                                          param_wbl(site, 2, period));
                end
                if isequal(DIST, 'Weibull-GPD')
                    if (param_gpd(site, 1, period) ~= 9999)
                        k = sim_tns == period & sim(:, site) > PARETO_THR;
                        sim(k, site) = gpinv((ran(k, site) - ...
                                              fit_thr(period, site)) / ...
                                             (1 - fit_thr(period, site)), ...
                                             param_gpd(site, 1, period), ...
                                             param_gpd(site, 2, period), ...
                                             PARETO_THR);
                    end
                end
            % Gamma distribution
            elseif isequal(DIST, 'Gamma')
                if (param_gam(site, 1, period) ~= 9999)
                    sim(k, site) = gaminv(ran(k, site), ...
                                          param_gam(site, 1, period), ...
                                          param_gam(site, 2, period));
                end
            end

            % Reshuffle precipitation so that inter-site 
            % correlations are maintained
            sim_rs_site = sim_rs_per(sim_rs_per(:, site) > 0, :);
            sim_rs_site = sortrows(sim_rs_site, site);
            uv = unique(sim_rs_site(:, site));  % uv = unique vector
            % Reshuffle within the same amounts to avoid bias (increasing
            % values towards the end of the simulated time series)
            for j = 1:length(uv)
                vec = sim_rs_site(sim_rs_site(:, site) == uv(j), NOS + 1);
                sim_rs_site(sim_rs_site(:, site) == uv(j), NOS + 1) = ...
                vec(randperm(length(vec)), 1);
            end

            sim_site = sim_per(sim_per(:, site) > 0, site);
            sim_site = sort(sim_site);

            sim_rs_site(:, site) = sim_site;
            sim_rs_site = sortrows(sim_rs_site, NOS + 1);

            sim_rs_per(sim_rs_per(:, site) > 0, site) = sim_rs_site(:, site);

        end  % Sites

        sim(sim_tns == period, :) = sim_rs_per(:, 1:NOS);

    end  % Periods
 
end  % Parametric precipitation


%% Write out the generated precipitation time series

disp('Saving the generated precipitation time series...')

sim = round(sim, 1);
writetable(table(sim_dt, sim), 'sim.csv', 'WriteVariableNames', false);

disp('Simulation finished!')

clearvars -except obs sim 





