%% Constants
global J k_t k_b N R L b m g l
J = .02;
k_t = .01;
k_b = k_t;
N = 10;
R = 2;
L = .5;
b = .2;
m = 1;
g = 10;
l = .3;

%% Script
%% Part a)
% Finding equilibrium points
for i = [15, 21.2, 30, 40, 75]
    fprintf("Equilibrium for u(t) = %f\n", i)
    tmp = equilibrium_constant(i);
    fprintf("x_1: %g + %gi\nx_2: %g + %gi\nx_3: %g + %gi\n", real(tmp(1)), imag(tmp(1)),real(tmp(2)), imag(tmp(2)),real(tmp(3)), imag(tmp(3)));
end


%% Part b)
% For accessing the linear system cells

linear_systems = {};
for i = [15, 45]
    fprintf("Linearization for u(t) = %f\n", i)
    [A B C D] = linearize(i);
    linear_systems{end+1} = {ss(A, B, C, D), equilibrium_constant(i), i};
end



%% Part c)
t = 0:0.005:2;
u = zeros(size(t));
system_conditions = {{10, [pi/2 0 0]},{21.2, [0 0 0]},{21.2, [pi/4 0 0]},{21.2, [pi/3 0 0]},{21.2, [pi/2 0 0]},{30, [-pi/2 0 0]},{30, [0 0 0]}};
for cond = system_conditions
    for sys = linear_systems
        % convert coordinates to the local linear approx coordintaes
        u(:,:) = cond{1}{1} - sys{1}{3};
        y = lsim(sys{1}{1}, u, t, cond{1}{2} - sys{1}{2}.');
        plot(t, y(:,1) + sys{1}{2}(1));
        plot_title = sprintf("lsim of system linearized around input = %f \n u(t) = %f  x_0 = %f", sys{1}{3}, cond{1}{1}, cond{1}{2}(1));
        title(plot_title);
        xlabel("Time");
        ylabel("Angle in radians");
        filename = sprintf("images/linearized_%f_input_%f_initial_pos_%f.png", sys{1}{3}, cond{1}{1}, cond{1}{2}(1));
        saveas(gcf, filename);
    end
end

%% Part d)
t = 0:0.005:2;
u = zeros(size(t));
system_conditions = {{10, [pi/2 0 0]},{21.2, [0 0 0]},{21.2, [pi/4 0 0]},{21.2, [pi/3 0 0]},{21.2, [pi/2 0 0]},{30, [-pi/2 0 0]},{30, [0 0 0]}, {5, [pi/2 0 0]}, {50, [0 0 0]}, {80, [0 0 0]}};
for cond = system_conditions
    [ts, y] = ode45(@(t,x) dx(x, @(x) cond{1}{1}, t), [t(1) t(end)], cond{1}{2});
    plot(ts, y(:,1));
    hold on

    legend_strings = ["ode45"];
    for sys = linear_systems
        % convert coordinates to the local linear approx coordintaes
        u(:,:) = cond{1}{1} - sys{1}{3};
        y = lsim(sys{1}{1}, u, t, cond{1}{2} - sys{1}{2}.');
        plot(t, y(:,1) + sys{1}{2}(1));
        legend_strings(end+1) = sprintf("Linearized (u(t)=%f)", sys{1}{3});
    end
    legend("location", "best");
    legend(legend_strings);

    plot_title = sprintf("Comparison of lsim to ode45 \n u(t) = %f  x_0 = %f", cond{1}{1}, cond{1}{2}(1));
    title(plot_title);
    xlabel("Time");
    ylabel("Angle in radians");

    filename = sprintf("images/input_%f_initial_pos_%f_vs_ode45.png", cond{1}{1}, cond{1}{2}(1));
    saveas(gcf, filename);
    hold off
end

%% Part f
% for functions see the bottom of the functions section
system_conditions = {{10, [pi/2 0 0]},{21.2, [0 0 0]},{21.2, [pi/4 0 0]},{21.2, [pi/3 0 0]},{21.2, [pi/2 0 0]},{30, [-pi/2 0 0]},{30, [0 0 0]}};
dts = [.01 .05, .1];
for dt = dts
    t = 0:dt:5;
    for cond = system_conditions
        [ts, y] = ode45(@(t,x) dx(x, @(x) cond{1}{1}, t), [t(1) t(end)], cond{1}{2});
        plot(ts, y(:,1), "r");

        hold on

        y = sim_rect(@(t,x) dx(x, @(x) cond{1}{1}, t), t, cond{1}{2});
        plot(t, y(:,1), "g");

        y = sim_trap(@(t,x) dx(x, @(x) cond{1}{1}, t), t, cond{1}{2});
        plot(t, y(:,1), "b");

        y = sim_runge_kutta(@(t,x) dx(x, @(x) cond{1}{1}, t), t, cond{1}{2});
        plot(t, y(:,1), "m");

        plot_title = sprintf("Comparison of numerical methods to ode45 \n dt = %f u(t) = %f  x_0 = %f", dt, cond{1}{1}, cond{1}{2}(1));
        title(plot_title);

        legend("location", "best");
        legend("ode45", "rectangular", "trapezoidal", "Runge Kutta");
        xlabel("Time");
        ylabel("Angle in radians");

        filename = sprintf("images/numerical_approx_dt_%f_input_%f_initial_pos_%f_vs_ode45.png", t(2)-t(1), cond{1}{1}, cond{1}{2}(1));
        saveas(gcf, filename);

        hold off
    end
end




%% Functions
function x_dot = dx(x,func_impulse, t)
    global J k_t k_b N R L b m g l
    % dx_1 = x_2
    % dx_2 = (1/J) * (N*k_t*x_3-b*x_2-m*g*l*sin(x_1))
    % dx_3 = (1/L) * (u(t) - R*x_3-k_b*N*x_2)
    x_dot = [
        x(2);
        (1/J) * (N*k_t*x(3)-b*x(2)-m*g*l*sin(x(1)))
        (1/L) * (func_impulse(t) - R*x(3)-k_b*N*x(2))
    ];
end

% returns the equilibrium given a constant input function
function eq = equilibrium_constant(u_t)
    global J k_t k_b N R L b m g l
    eq = [
    asin((N*k_t*u_t)/(R*m*g*l));
    0;
    u_t/R;
    ];
end

% Linearize given an constant input function
function [A B C D] = linearize(u_t)
    global J k_t k_b N R L b m g l
    eq = equilibrium_constant(u_t);
    A = [
    0 1 0;
    -m*g*l*cos(eq(1))/J -b/J N*k_t/J;
    0 -k_b*N/L -R/L;
    ];
    B = [
    0;
    0;
    1/L;
    ];
    C = [
    1 0 0;
    0 0 0;
    0 0 0;
    ];
    D = [0];
end

% Approximate the simulations using the rectangular method of integration
function y = sim_rect(func_derivative_t_x, t, initial_condition)
    % Allocate space for all the results
    y = zeros(size(t,2),3);
    y(1, :) = initial_condition;
    for i = 2:size(y,1)
        dt = t(i)-t(i-1);
        y(i,:) = y(i-1,:).' + func_derivative_t_x(t(i), y(i-1,:))*dt;
    end
end

% Approximate the simulations using the rectangular method of integration
function y = sim_trap(func_derivative_t_x, t, initial_condition)
    % Allocate space for all the results
    y = zeros(size(t,2),3);
    y(1, :) = initial_condition;

    delta_x = zeros(2,3);
    for i = 2:size(y,1)
        dt = t(i)-t(i-1);
        delta_x(1,:) = dt*func_derivative_t_x(t(i-1), y(i-1,:));
        delta_x(2,:) = dt*func_derivative_t_x(t(i), y(i-1,:) + delta_x(1,:));
        y(i,:) = y(i-1,:) + (delta_x(1,:) + delta_x(2,:))/2;
    end
end

% Approximate the simulations using the rectangular method of integration
function y = sim_runge_kutta(func_derivative_t_x, t, initial_condition)
    dt = t(2)-t(1);
    y = zeros(size(t,2),3);
    y(1, :) = initial_condition;

    k = zeros(4,3);
    for i = 2:size(y,1)
        k(1,:) = func_derivative_t_x(t(i-1),y(i-1,:));
        k(2,:) = func_derivative_t_x(t(i-1)+dt/2, y(i-1,:) + dt*k(1,:)/2);
        k(3,:) = func_derivative_t_x(t(i-1)+dt/2, y(i-1,:) + dt*k(2,:)/2);
        k(4,:) = func_derivative_t_x(t(i), y(i-1,:) + dt*k(3,:));
        y(i,:) = y(i-1,:) + dt/6 * (k(1,:) + 2*k(2,:) + 2*k(3,:) + k(4,:));
    end
end
