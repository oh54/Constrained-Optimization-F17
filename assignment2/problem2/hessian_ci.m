function [ val ] = hessian_ci( x, i )
    syms x1 x2 x3 x4 x5
    
    if (i == 1)
        c = x1^2 + x2^2 + x3^2 + x4^2 + x5^2 - 10;
    else if (i == 2)
        c = x2*x3 - 5*x4*x5;
        else if (i == 3)
            c = x1^3 + x2^3 + 1;
            end
        end
    end
    H = hessian(c, [x1, x2, x3, x4, x5]);
    xCell = num2cell(x);
    [x1, x2, x3, x4, x5] = xCell{:};
    val = eval(subs(H));
end

