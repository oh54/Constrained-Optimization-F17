clear
clc 
close all

%% Compare performance of LDL, LDL sparse, LU, LU sparse, Null-space, Range-space solvers
u_bar = 0.2;
d_0 = 1;
clc

ns = 10:10:1000;
times = zeros(6, size(ns,2));

name = 'Data/problem_2.mat';
if exist(name, 'file') 
    disp('Loading data..')
    load(name)
else

    disp('Iterating...')
    i = 1;
    for n = ns

        % Print n every 100 iterations
        if mod(n,100) == 0
            str = sprintf('n = %d', n);
            disp(str)
        end

        start_time = cputime;
        [x, lambda] = EqualityQPSolverLDL(n,u_bar,d_0);
        times(1,i) = cputime-start_time;

        start_time = cputime;
        [x, lambda] = EqualityQPSolverLU(n,u_bar,d_0);
        times(2,i) = cputime-start_time;

        start_time = cputime;
        [x, lambda] = NullSpaceQR(n,u_bar,d_0);
        times(3,i) = cputime-start_time;

        start_time = cputime;
        [x, lambda] = RangeSpaceSolver(n,u_bar,d_0);
        times(4,i) = cputime-start_time;    

        start_time = cputime;
        [x, lambda] = EqualityQPSolverLDL(n,u_bar,d_0, true);
        times(5,i) = cputime-start_time;

        start_time = cputime;
        [x, lambda] = EqualityQPSolverLU(n,u_bar,d_0, true);
        times(6,i) = cputime-start_time;
        
        i = i + 1;
    end
    save(name, 'times')
end
disp('Finished')

% Geranate plots
plot(ns,times(1,:), 'ro')
hold on
plot(ns,times(2,:), 'bs')
plot(ns,times(3,:), 'g*')
plot(ns,times(4,:), 'k.')
plot(ns,times(5,:), 'mx')
plot(ns,times(6,:), 'cd')
xlabel('Problem size: n')
ylabel('CPU time')
legend('LDL solver','LU solver','Null Space solver','Range Space solver', 'LDL sparse solver','LU sparse solver')

%% Plot sparsity pattern of KKT matrix
n = 100;
KKT = ConstructKKT(n, u_bar, d_0);
spy(KKT);
