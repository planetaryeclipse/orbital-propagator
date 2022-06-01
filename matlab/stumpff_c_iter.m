function [c] = stumpff_c_iter(z, n)
    c = 0;
    for k=0:n
        c = c + (-1)^k * z^k / factorial(2 * k + 2);
    end
end