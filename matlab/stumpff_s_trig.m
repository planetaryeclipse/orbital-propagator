function [s] = stumpff_s_trig(z)
    if z < 0
        s = (sinh(sqrt(-z)) - sqrt(-z)) / sqrt(-z)^3;
    elseif z == 0
        s = 1 / 6;
    else % z > 0
        s = (sqrt(z) - sin(sqrt(z))) / sqrt(z)^3;
    end
end