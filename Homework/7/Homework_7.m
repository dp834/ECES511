%% Globals
global a b r T
a = -2;
b = 5;
r = 1;
T = 1/8;

csys = ss(a,b,1,0);
dsys = c2d(csys, T, "zoh")

a_d = exp(-1/4);
b_d = (1-a_d)*5/2;
recursive_solve_to_n(a_d, b_d, 0, 10)


%% Functions
function x = recursive_solve_to_n(a_d, b_d, x0, n)
    x = zeros(n,1);
    x(1) = x0;
    i = 2;
    x = private_recursive_solve_to_n(a_d,b_d, n, x, i);

end

function x = private_recursive_solve_to_n(a_d, b_d, n, x, i)
    if (i > n)
        return
    end
    x(i) = a_d*x(i-1) + b_d;
    x = private_recursive_solve_to_n(a_d,b_d, n, x, i+1);
end
