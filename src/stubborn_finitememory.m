% STUBBORN_FINITEMEMORY Updates the distributed finite memory stubborn observer state.
% It computes the next observer state estimation based on the current measure,
% input, and dynamic saturation threshold computed from its finite memory.
function [x_hat_next, sigma_next, Delta, y_hat] = stubborn_finitememory(x_hat, sigma, y, u, plant)
    y_hat = plant.C * x_hat + plant.D * u;
    innov = y - y_hat;
    Delta = plant.R * sigma;
    sat_innov = sign(innov) .* min(abs(innov), Delta);
    x_hat_next = plant.A * x_hat + plant.B * u + plant.L * sat_innov;
    sigma_next = plant.S * sigma + plant.G * abs(innov);
end