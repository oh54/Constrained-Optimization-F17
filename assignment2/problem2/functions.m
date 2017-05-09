function [ val ] = f( x )
    syms x1 x2 x3 x4 x5
    f = exp(x1*x2*x3*x4*x5) - 1/2 * (x1^3 + x2^3 + 1)^2;
    xCell = num2cell(x);
    [x1, x2, x3, x4, x5] = xCell{:};
    val = eval(subs(f));
end

function [ val ] = gradient_f( x )
    syms x1 x2 x3 x4 x5
    f = exp(x1*x2*x3*x4*x5) - 1/2 * (x1^3 + x2^3 + 1)^2;
    nabla_f = gradient(f, [x1, x2, x3, x4, x5]);
    xCell = num2cell(x);
    [x1, x2, x3, x4, x5] = xCell{:};
    val = eval(subs(nabla_f));
end

function [ val ] = hessian_f( x )
    syms x1 x2 x3 x4 x5
    f = exp(x1*x2*x3*x4*x5) - 1/2 * (x1^3 + x2^3 + 1)^2;
    nabla2_f = hessian(f, [x1, x2, x3, x4, x5]);
    xCell = num2cell(x);
    [x1, x2, x3, x4, x5] = xCell{:};
    val = eval(subs(nabla2_f));
end

function [ val ] = c( x )
    syms x1 x2 x3 x4 x5
    c1 = x1^2 + x2^2 + x3^2 + x4^2 + x5^2 - 10;
    c2 = x2*x3 - 5*x4*x5;
    c3 = x1^3 + x2^3 + 1;
    xCell = num2cell(x);
    [x1, x2, x3, x4, x5] = xCell{:};
    val = [ eval(subs(c1)) eval(subs(c2)) eval(subs(c3)) ]';
end

function [ val ] = jacobian_c( x )
    syms x1 x2 x3 x4 x5
    c = [ x1^2 + x2^2 + x3^2 + x4^2 + x5^2 - 10 ; x2*x3 - 5*x4*x5 ; x1^3 + x2^3 + 1];
    xCell = num2cell(x);
    [x1, x2, x3, x4, x5] = xCell{:};
    val = eval(subs(jacobian(c)));
end

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


