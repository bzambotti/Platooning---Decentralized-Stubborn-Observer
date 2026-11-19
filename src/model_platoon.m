% MODEL_PLATOON  Discrete-time (ZOH) model of a homogeneous platoon of N followers
% under a CTH spacing policy and a predecessor-following CACC law.
% State of follower i: [spacing error; velocity; acceleration]. Output: C = I_N kron [1 0 0; 0 1 0].
% Input u = [v_0; a_0] (leader velocity and acceleration).
function [A, B, C] = model_platoon(params)
    N = params.N; % Number of followers
    n_states = 3 * N;
    dt = params.dt;
    
    % Vehicle and controller parameters
    tau = params.tau;
    kp = params.kp; 
    kd = params.kd; 
    ka = params.ka; % acceleration feedback gain
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
        A_cont(a_idx, a_idx) = -(1 + ka)/tau;
        
        % Coupling with the preceding vehicle (i-1)
        if i > 1
            prev_idx = (i-2)*3 + 1 : (i-1)*3;
            A_cont(d_idx, prev_idx(2)) = 1;       % v_prev in e_dot
            A_cont(a_idx, prev_idx(2)) = kd/tau;  % v_prev in u
            A_cont(a_idx, prev_idx(3)) = ka/tau;  % a_prev in u
        end
    end
    
    % B_cont: exogenous input (the leader) acting on the first follower
    % U_leader = [v_0; a_0] -> B_cont has 2 columns
    B_cont = zeros(n_states, 2); 
    
    % Column 1: terms multiplied by v_0
    B_cont(1, 1) = 1;       % v_0 enters e_dot of the first follower
    B_cont(3, 1) = kd/tau;  % v_0 enters the control law (a_dot) of the first follower
    
    % Column 2: terms multiplied by a_0
    B_cont(3, 2) = ka/tau;  % a_0 enters the control law (a_dot) of the first follower
    
    % Output matrix C (measures e_i, v_i of each follower)
    C = kron(eye(N), [1 0 0; 0 1 0]); 
    
    % No direct feedthrough
    D_cont = zeros(size(C, 1), 2);
   
    % --- DISCRETIZZAZIONE ---
    sys_c = ss(A_cont, B_cont, C, D_cont); 
    sys_d = c2d(sys_c, dt, 'zoh'); 
    
    A = sys_d.A; 
    B = sys_d.B; 
    % C is unchanged by the ZOH discretization
end
