function [s] = stumpff_s_iter(z, n)
    s = 0;
    for k=0:n
        s = s + (-1)^k * z^k / factorial(2 * k + 3);
    end
end