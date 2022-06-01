function [a] = calc_semimajor(r, v, mu)
    a = 1 / (2 / norm(r) - norm(v)^2/mu);
end