function [ val ] = jacobian_c( x )
    syms x1 x2 x3 x4 x5
    c = [ x1^2 + x2^2 + x3^2 + x4^2 + x5^2 - 10 ; x2*x3 - 5*x4*x5 ; x1^3 + x2^3 + 1];
    xCell = num2cell(x);
    [x1, x2, x3, x4, x5] = xCell{:};
    val = eval(subs(jacobian(c)));
end

