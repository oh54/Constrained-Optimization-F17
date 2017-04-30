function [ val ] = f( x )
    syms x1 x2 x3 x4 x5
    f = exp(x1*x2*x3*x4*x5) - 1/2 * (x1^3 + x2^3 + 1)^2;
    xCell = num2cell(x);
    [x1, x2, x3, x4, x5] = xCell{:};
    val = eval(subs(f));
end

