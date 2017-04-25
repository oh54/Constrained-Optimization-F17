function [ck,Ak] = constraint_fun(x)

c1x1 = 2*x(1)+4;
c1x2 = -1; 
c2x1 = -4; 
c2x2 = 10;

Ak = [c1x1,c2x1;c1x2,c2x2];

c1 = (x(1)+2).^2-x(2);
c2 = -4*x(1)+10*x(2);

ck = [c1;c2];

end