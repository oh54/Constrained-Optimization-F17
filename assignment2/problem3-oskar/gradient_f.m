function [ val ] = gradient_f( x )
    syms x1 x2
    f = (x1^2 + x2 - 11)^2 + (x1 + x2^2 -7)^2;
    nabla_f = gradient(f, [x1, x2]);
    xCell = num2cell(x);
    [x1, x2] = xCell{:};
    val = eval(subs(nabla_f));
end

