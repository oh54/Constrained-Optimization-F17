function [ val ] = gradient_f( x )
    syms x1 x2 x3 x4 x5
    f = exp(x1*x2*x3*x4*x5) - 1/2 * (x1^3 + x2^3 + 1)^2;
    nabla_f = gradient(f, [x1, x2, x3, x4, x5]);
    xCell = num2cell(x);
    [x1, x2, x3, x4, x5] = xCell{:};
    val = eval(subs(nabla_f));
end

