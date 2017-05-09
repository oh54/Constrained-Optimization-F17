function [ val ] = f( x )
    syms x1 x2
    f = (x1^2 + x2 - 11)^2 + (x1 + x2^2 -7)^2;
    xCell = num2cell(x);
    [x1, x2] = xCell{:};
    val = eval(subs(f));
end

