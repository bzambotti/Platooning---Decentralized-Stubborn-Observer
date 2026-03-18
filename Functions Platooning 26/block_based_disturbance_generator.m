function dist_hist = block_based_disturbance_generator(Pd, n_part, f_intensities, Tsim, scenario, tau_cent, target_idx)
    % BLOCK_BASED_DISTURBANCE_GENERATOR
    % Genera attacchi sparsi basati su scenari specifici.
    
    if nargin < 5, scenario = 'safe'; end
    if nargin < 6, tau_cent = sum(n_part); end
    if nargin < 7, target_idx = []; end
    
    total_sensors = length(f_intensities);
    dist_hist = zeros(total_sensors, Tsim );
    num_blocks = length(Pd);
    
    fprintf('>> Generazione Disturbi. Scenario: %s\n', upper(scenario));
    
    for j = 1:num_blocks
        sensors_in_block = Pd{j};
        % n_part(j) è la dimensione memoria locale (tau_local)
        tau_local = n_part(j); 
        
        % Se il blocco non ha intensità di attacco, saltalo
        % Nota: controlliamo se TUTTI i sensori del blocco sono a 0
        if all(f_intensities(sensors_in_block) == 0)
             continue; 
        end
        
        % --- DETERMINAZIONE LIMITI DWELL TIME [Min, Max] ---
        switch lower(scenario)
            case 'safe'
                % Deve essere maggiore sia del tempo di recupero locale che globale
                min_dwell = max(tau_local, tau_cent) + 5; 
                max_dwell = min_dwell + 30; 
                
            case 'chaos'
                % Viola il dwell time locale (attacca molto spesso)
                min_dwell = 1;
                % Se tau_local è 1, max deve essere almeno 1
                max_dwell = max(1, tau_local - 1);
                
            case 'central_violation'
                % Viola SOLO il centralizzato (tau_cent), ma rispetta il locale
                % Intervallo: [tau_local + 1, tau_cent - 1]
                min_dwell = tau_local + 2; 
                max_dwell = max(min_dwell, tau_cent - 2);
                
                % Se tau_cent è troppo piccolo rispetto a tau_local, impossibile violare selettivamente
                if min_dwell > max_dwell
                    warning('Scenario Central Violation impossibile per blocco %d (tau_loc=%d, tau_cent=%d). Uso Safe.', j, tau_local, tau_cent);
                    min_dwell = max(tau_local, tau_cent) + 2;
                    max_dwell = min_dwell + 10;
                end
                
            case 'target_violation'
                % Se il blocco contiene il sensore bersaglio, attacca aggressivamente
                if any(ismember(sensors_in_block, target_idx))
                    min_dwell = 1;
                    max_dwell = max(1, tau_local - 1);
                else
                    % Altri blocchi SAFE
                    min_dwell = tau_local + 5;
                    max_dwell = tau_local + 20;
                end
                
            otherwise
                error('Scenario sconosciuto: %s', scenario);
        end
        
        % --- GENERAZIONE SEQUENZA TEMPORALE ---
        % Start casuale per desincronizzare i blocchi
        current_k = randi([1, 10]); 
        
        while current_k <= Tsim
            
            % 1. Inietta disturbo sui sensori del blocco
            % Solo per i sensori che hanno intensità > 0
            for s_id = sensors_in_block
                % Usa s_id come indice lineare se Pd contiene indici assoluti (1...m)
                if s_id > total_sensors
                    continue; % Protezione indici
                end
                
                f_bar = f_intensities(s_id);
                
                if f_bar > 0
                    % Randomizza magnitudine
                    sign_val = sign(randn); if sign_val==0, sign_val=1; end
                    magnitude = f_bar * (0.9 + 0.2 * rand); % +/- 10%
                    
                    dist_hist(s_id, current_k) = sign_val * magnitude;
                end
            end
            
            % 2. Calcola prossimo passo (Dwell Time dinamico)
            % CORREZIONE QUI: Parentesi quadre per il range
            step = randi([min_dwell, max_dwell], 1);
            
            current_k = current_k + step;
       end
    end
end