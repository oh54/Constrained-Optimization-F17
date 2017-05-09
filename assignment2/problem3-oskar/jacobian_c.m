function [ val ] = jacobian_c( x )
    syms x1 x2
    c = [ (x1 + 2)^2 - x2 ; -4*x1 + 10*x2 ];
    xCell = num2cell(x);
    [x1, x2] = xCell{:};
    val = eval(subs(jacobian(c)));
end

