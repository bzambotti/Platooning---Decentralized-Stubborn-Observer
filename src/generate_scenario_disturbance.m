function dist_hist = generate_scenario_disturbance(params, num_steps, scenario_type)
% GENERATE_SCENARIO_DISTURBANCE
% Generates the outlier sequences used in the simulation scenarios.
%   'spatial_isolation'   : sporadic outliers with random dwell time in [min_dwell, max_dwell]
%   'temporal_relaxation' : outliers with fixed period T_attack
%   'safety_critical'     : periodic positive outliers in an active window (not used in the paper)
% Outlier magnitude: params.magnitude * U[0.5, 1.5], random sign.

n_meas = size(params.C, 1);
dist_hist = zeros(n_meas, num_steps);

if isfield(params, 'target_sensor_idx')
    target_idx = params.target_sensor_idx;
else
    error('params.target_sensor_idx is required');
end

switch scenario_type
    case 'spatial_isolation'
        % Scenario 1: Sporadic high-magnitude spikes.
        
        % Defaults
        if ~isfield(params, 'min_dwell'), min_d = 1; else, min_d = params.min_dwell; end
        if ~isfield(params, 'max_dwell'), max_d = 1; else, max_d = params.max_dwell; end
        
        % --- Safety check ---
        % min_d >= 1
        min_d = max(1, min_d);
        
        % max_d >= min_d (e.g. when n_local = 0)
        if max_d < min_d
            max_d = min_d; 
        end
        % -----------------------------------------------
        
        k = 1;
        while k <= num_steps
            % Magnitude and sign
            mag_factor = 0.5 + 1.0 * rand; 
            magnitude = params.magnitude * mag_factor;
            if rand < 0.5, s = 1; else, s = -1; end
            
            % Apply outlier
            dist_hist(target_idx, k) = magnitude * s;
            
            % Next attack time (dwell time)
            step = randi([min_d, max_d]);
            k = k + step;
        end
        
    case 'temporal_relaxation'
        % Scenario 2: Sporadic spikes with fixed period T_attack.
        if ~isfield(params, 'T_attack'), error('T_attack required'); end
        T_attack = params.T_attack; 
        
        start_k = randi([1 5]); 
        for k = start_k:T_attack:num_steps
             mag_factor = 0.5 + 1.0 * rand;
             magnitude = params.magnitude * mag_factor;
             if rand < 0.5, s = 1; else, s = -1; end
             dist_hist(target_idx, k) = magnitude * s;
        end
        
    case 'safety_critical'
        % Scenario 3: Braking under High-Freq Fire.
        if ~isfield(params, 'T_attack'), error('T_attack required'); end
        T_attack = max(1, round(params.T_attack));
        
        if isfield(params, 'active_window')
            k_start = params.active_window(1);
            k_end = params.active_window(2);
        else
            k_start = 1; k_end = num_steps;
        end
        
        for k = k_start:T_attack:k_end
             % Positive outliers only (spoofed distance)
             mag_factor = 0.8 + 0.4 * rand; 
             magnitude = params.magnitude * mag_factor;
             dist_hist(target_idx, k) = magnitude; 
        end
        
    otherwise
        error('Unknown scenario type: %s', scenario_type);
end
end
