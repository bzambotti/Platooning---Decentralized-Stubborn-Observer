function [A, B, C] = model_platoon(params)
    % MODEL_PLATOON Computes the discrete-time state-space representation
    % of the platoon system based on Constant Time Headway (CTH) policy.
    %
    % Inputs:
    %   params - Struct containing simulation parameters:
    %            N: Number of followers in the platoon.
    %            dt: Discretization time step.
    %            tau, h, kp, kd, ka: Control and physics parameters.
    %
    % Outputs:
    %   A: Discrete-time state transition matrix.
    %   B: Discrete-time input matrix.
    %   C: Measurement matrix for tracking spacing errors and velocities.

    N = params.N; % Numero di follower
    n_states = 3 * N;
    dt = params.dt;

    % Parametri del veicolo e del controllo
    tau = params.tau;
    kp = params.kp;
    kd = params.kd;
    ka = params.ka; % <-- Aggiunto per il controllo sull'accelerazione
    h = params.h;

    A_cont = zeros(n_states, n_states);

    for i = 1:N
        idx = (i-1)*3 + 1 : i*3;
        d_idx = idx(1); v_idx = idx(2); a_idx = idx(3);

        % 1. Spacing Error dynamics: e_dot = v_prev - v - h*a
        A_cont(d_idx, v_idx) = -1;
        A_cont(d_idx, a_idx) = -h; % CTH term

        % 2. Velocity dynamics: v_dot = a
        A_cont(v_idx, a_idx) = 1;

        % 3. Acceleration dynamics: a_dot = -1/tau * a + 1/tau * u
        A_cont(a_idx, d_idx) = kp/tau;
        A_cont(a_idx, v_idx) = -kd/tau;
        A_cont(a_idx, a_idx) = -(1 + ka)/tau; % <-- Modificato per includere il feedback ka

        % Accoppiamento con il veicolo precedente (i-1)
        if i > 1
            prev_idx = (i-2)*3 + 1 : (i-1)*3;
            A_cont(d_idx, prev_idx(2)) = 1;       % v_prev in e_dot
            A_cont(a_idx, prev_idx(2)) = kd/tau;  % v_prev in u
            A_cont(a_idx, prev_idx(3)) = ka/tau;  % a_prev in u
        end
    end

    % B_cont rappresenta l'ingresso esogeno (il LEADER) sul primo follower
    % U_leader = [v_0; a_0] -> B_cont avra' 2 colonne
    B_cont = zeros(n_states, 2);

    % Colonna 1: ingressi moltiplicati per v_0
    B_cont(1, 1) = 1;       % v_0 entra in e_dot del primo follower
    B_cont(3, 1) = kd/tau;  % v_0 entra nel controllo (a_dot) del primo follower

    % Colonna 2: ingressi moltiplicati per a_0
    B_cont(3, 2) = ka/tau;  % a_0 entra nel controllo (a_dot) del primo follower

    % Matrice di Uscita C (Misura e_i, v_i per ogni follower)
    C = kron(eye(N), [1 0 0; 0 1 0]);

    % Matrice D (Nessun feedthrough diretto)
    D_cont = zeros(size(C, 1), 2);

    % --- DISCRETIZZAZIONE ---
    % Usa direttamente la tua matrice C così il sistema è già formattato bene
    sys_c = ss(A_cont, B_cont, C, D_cont);
    sys_d = c2d(sys_c, dt, 'zoh');
    A = sys_d.A;
    B = sys_d.B;
    % C rimane invariata nella discretizzazione ZOH
end
