% Decentralized finite-memory stubborn observer (one threshold per output block).
function [x_hat_next, sigma_next,Delta, y_hat] = stubborn_finitememory(x_hat, sigma, y, u, plant)
    y_hat = plant.C * x_hat + plant.D * u;
    innov = y - y_hat;
    Delta = plant.R * sigma; 
    sat_innov = sign(innov) .* min(abs(innov), Delta); 
    x_hat_next = plant.A * x_hat + plant.B * u + plant.L * sat_innov;
    sigma_next = plant.S * sigma + plant.G * abs(innov); 
end