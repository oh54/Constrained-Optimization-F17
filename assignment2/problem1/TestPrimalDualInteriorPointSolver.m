%% Test PrimalDualInteriorPointSolver
%
%  The test script is based on the one provided during classes in week 7
%

% These values are of arbitrary choice and define the size of the problem
n = 200;
m = 50;

% For replicability define the state of random seed
state = 1000;
rand('state',state);

% Create random matrix A
A = randn(m,n);

% Define the solution to the problem 
x = zeros(n,1);

x(1:m,1) = abs(rand(m,1));

% Define multipliers
lambda = zeros(n,1);
lambda(m+1:n,1) = abs(rand(n-m,1));
mu = rand(m,1);

% Generate the objective function
g = A'*mu + lambda;

% And the rhs of the constraints
b = A*x;

% x0 is a starting point
x0 = ones(n,1);

% Call the algorithm
[xlp, mulp, lambdalp, Converged] = PrimalDualInteriorPointSolver(g,A,b,x0);

% If converged calculate the difference between true and calculated
% solutions
if Converged
    X = max(abs(xlp-x))
    M = max(abs(mulp-mu))
    L = max(abs(lambdalp-lambda))
end
