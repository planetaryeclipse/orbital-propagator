function [c] = stumpff_c_trig(z)
    if z < 0
        c = (cosh(sqrt(-z)) - 1) / (-z);
    elseif z == 0
        c = 1 / 2;
    else % z > 0
        c = (1 - cos(sqrt(z))) / z;
    end
end

% fprintf("%.*f\n", 8, stumpff_c_trig(-15))