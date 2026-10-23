function r = inverse_fisher_z(z)
    r = (exp(2*z) - 1) ./ (exp(2*z) + 1);
end
