%% Fisher's z-transformation and its inverse
function z = fisher_z(r)
    % Handle NaN values
    z = r;
    valid = ~isnan(r);
    
    % Apply Fisher's z-transformation
    z(valid) = 0.5 * log((1 + r(valid)) ./ (1 - r(valid)));
end

function r = inverse_fisher_z(z)
    % Handle NaN values
    r = z;
    valid = ~isnan(z);
    
    % Apply inverse Fisher's z-transformation
    r(valid) = (exp(2 * z(valid)) - 1) ./ (exp(2 * z(valid)) + 1);
end