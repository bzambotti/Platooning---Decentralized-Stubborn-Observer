% function dist_hist = generate_scenario_disturbance(params, num_steps, scenario_type)
% % GENERATE_SCENARIO_DISTURBANCE
% % Generates disturbance sequences for specific CDC26 scenarios.
% %
% % params: struct containing:
% %   - C: Output matrix (to determine number of sensors)
% %   - target_sensor_idx: Index of the sensor to attack (e.g., 3 for Veh 2 Dist)
% %   - T_attack: Attack period (for scenarios 2 and 3)
% %   - magnitude: Magnitude of the attack
% %   - min_dwell, max_dwell: For random sporadic attacks (Scenario 1)
% %   - active_window: [start_step, end_step] for Scenario 3
% %
% % num_steps: Total simulation steps
% % scenario_type: 'spatial_isolation', 'temporal_relaxation', 'safety_critical'
% 
% n_meas = size(params.C, 1);
% dist_hist = zeros(n_meas, num_steps);
% 
% if isfield(params, 'target_sensor_idx')
%     target_idx = params.target_sensor_idx;
% else
%     error('params.target_sensor_idx is required');
% end
% 
% 
% switch scenario_type
%     case 'spatial_isolation'
%         % Scenario 1: Sporadic high-magnitude spikes.
%         % Dwell time: Safe for both (e.g., > max(n_local, n_global)).
% 
%         if ~isfield(params, 'min_dwell') || ~isfield(params, 'max_dwell')
%              error('min_dwell and max_dwell are required for spatial_isolation');
%         end
%         min_dwell = params.min_dwell;
%         max_dwell = params.max_dwell;
% 
%         % Start right away
%         k = 1;
%         while k <= num_steps
%             % Inject spike
% 
%             magnitude = params.magnitude * (0.2 + 1.6 * rand);
%             sign_val = sign(randn); if sign_val==0, sign_val=1; end
%             dist_hist(target_idx, k) = magnitude * sign_val;
% 
%             % Wait for dwell time
%             step = randi([min_dwell, max_dwell]);
%             k = k + step;
%         end
% 
%     case 'temporal_relaxation'
%         % Scenario 2: Sporadic spikes with period T_attack.
%         % "Generate sporadic spikes with a period T_attack".
% 
%         if ~isfield(params, 'T_attack')
%              error('T_attack is required for temporal_relaxation');
%         end
%         T_attack = params.T_attack;
% 
%         % Start after some time to let system settle
%         start_k = 1;
% 
%         for k = start_k:T_attack:num_steps
%              % Random sign to avoid constant bias
% 
%              magnitude = params.magnitude * (0.2 + 1.6 * rand);
%              sign_val = sign(randn); if sign_val==0, sign_val=1; end
%              dist_hist(target_idx, k) = magnitude * sign_val;
%         end
% 
%     case 'safety_critical'
%         % Scenario 3: Braking under High-Freq Fire.
%         % "Large positive impulses... contradict the braking".
% 
%         if ~isfield(params, 'T_attack')
%              error('T_attack is required for safety_critical');
%         end
%         T_attack = params.T_attack;
% 
%         % Attack window
%         if isfield(params, 'active_window')
%             k_start = params.active_window(1);
%             k_end = params.active_window(2);
%         else
%             k_start = 1;
%             k_end = num_steps;
%         end
% 
%         for k = k_start:T_attack:k_end
% 
%              magnitude = params.magnitude * (0.2 + 1.6 * rand);
%              dist_hist(target_idx, k) = magnitude; % POSITIVE IMPULSE (Vehicle appears farther)
%         end
% 
%     otherwise
%         error('Unknown scenario type: %s', scenario_type);
% end
% end
function dist_hist = generate_scenario_disturbance(params, num_steps, scenario_type)
% GENERATE_SCENARIO_DISTURBANCE
% Generates disturbance sequences for specific CDC26 scenarios.

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
        
        % Valori default se mancano
        if ~isfield(params, 'min_dwell'), min_d = 1; else, min_d = params.min_dwell; end
        if ~isfield(params, 'max_dwell'), max_d = 1; else, max_d = params.max_dwell; end
        
        % --- SAFETY CHECK (Il fix per il tuo errore) ---
        % Assicuriamo che min_d sia almeno 1
        min_d = max(1, min_d);
        
        % Assicuriamo che max_d sia sempre >= min_d
        % Se params dava max < min (es. calcoli con n_local=0), qui forziamo l'uguaglianza.
        if max_d < min_d
            max_d = min_d; 
        end
        % -----------------------------------------------
        
        k = 1;
        while k <= num_steps
            % Genera Ampiezza e Segno
            mag_factor = 0.5 + 1.0 * rand; 
            magnitude = params.magnitude * mag_factor;
            if rand < 0.5, s = 1; else, s = -1; end
            
            % Applica disturbo
            dist_hist(target_idx, k) = magnitude * s;
            
            % Calcola il salto (Dwell Time)
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
             % Solo positivi o variazione minore per simulare inganno distanza
             mag_factor = 0.8 + 0.4 * rand; 
             magnitude = params.magnitude * mag_factor;
             dist_hist(target_idx, k) = magnitude; 
        end
        
    otherwise
        error('Unknown scenario type: %s', scenario_type);
end
end