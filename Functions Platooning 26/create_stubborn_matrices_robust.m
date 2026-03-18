function [S, G, R] = create_stubborn_matrices_robust(n_vec, m_vec, r_vec)
    % n_vec: vector containing n_j (memory length for each block j)
    % m_vec: vector containing m_j (number of outputs for each block j)
    % r_vec: vector containing r_j (tuning parameter for each block j)
    
    s_blocks = length(m_vec);
    S_list = {}; G_list = {}; R_list = {};
    
    for j = 1:s_blocks
        nj = n_vec(j); 
        mj = m_vec(j); 
        rj = r_vec(j); 
        
        % 1. Construction S_chan (nj x nj)
        % The memory size depends only on the observability index nj,
        % not on the number of outputs mj.
        S_chan = zeros(nj, nj);
        if nj > 1
            idx_sub = sub2ind([nj, nj], 2:nj, 1:nj-1);
            S_chan(idx_sub) = 1;
        end
        
        % 2. Construction G_chan (nj x mj)
        % G must map the mj residuals into the first state of the memory.
        % We sum the errors (or take the norm) into the first shift register element.
        G_chan = zeros(nj, mj);
        if nj > 0
            % The first row collects errors from all mj channels
            G_chan(1, :) = 1; 
        end
        
        % 3. Construction R_chan (mj x nj)
        % R must map the nj memory states to thresholds for mj outputs.
        % Typically, the threshold is broadcast to all outputs in the block.
        % We assume the threshold depends on the sum of the memory buffer.
        R_chan = zeros(mj, nj);
        if nj > 0
             % Each row corresponds to an output channel. 
             % All outputs in this block share the same threshold logic.
             R_chan = rj * ones(mj, nj); 
        end
        
        % Store in lists (No Kronecker expansion for S!)
        S_list{j} = S_chan;
        G_list{j} = G_chan;
        R_list{j} = R_chan;
    end
    
    % Final assembly
    S = blkdiag(S_list{:});
    G = blkdiag(G_list{:});
    R = blkdiag(R_list{:});
end