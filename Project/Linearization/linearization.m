%% Constants
global R L m K g
R = 1;
L = 0.01;
m = 0.05;
K = 0.0001;
g = 9.81;

%% Script
%% Part a)
% Finding equilibrium points
for i = [1, 7, 25, 50]
    fprintf("Equilibrium for u(t) = %f\n", i)
    tmp = equilibrium_constant(i);
    fprintf("x_1: %g + %gi\nx_2: %g + %gi\nx_3: %g + %gi\n", real(tmp(1)), imag(tmp(1)),real(tmp(2)), imag(tmp(2)),real(tmp(3)), imag(tmp(3)));
end


%% Part b)
% For accessing the linear system cells

linear_systems = {};
for i = [1, 7, 25, 50]
    fprintf("Linearization for u(t) = %f\n", i)
    [A, B, C, D] = linearize(i);
    linear_systems{end+1} = {ss(A, B, C, D), equilibrium_constant(i), i};
end



%% Part c)
t = 0:0.005:2;
u = zeros(size(t));
system_conditions = {{5, [0 .01 0]}, {7, [0 .005 0]}, {10, [0 .03 0]}, {25, [0 .1 0]}, {linear_systems{2}{3}, linear_systems{2}{2}.'}};
for cond = system_conditions
    for sys = linear_systems
        % convert coordinates to the local linear approx coordintaes
        u(:,:) = cond{1}{1} - sys{1}{3};
        [~, ~, y] = lsim(sys{1}{1}, u, t, cond{1}{2} - sys{1}{2}.');
        y = y + sys{1}{2}.';
        subplot(3,1,1);
        plot(t,y(:,1));
        xlabel("Time");
        ylabel("Current in Amps");
        subplot(3,1,2);
        plot(t,y(:,2));
        xlabel("Time");
        ylabel("Distance in meters");
        subplot(3,1,3);
        plot(t,y(:,3));
        xlabel("Time");
        ylabel("Velocity in m/s");
        plot_title = sprintf("lsim of system linearized around input = %f \n u(t) = %f  x_0 = %f", sys{1}{3}, cond{1}{1}, cond{1}{2}(2));
        sgtitle(plot_title);
        filename = sprintf("images/linearized_%f_input_%f_initial_pos_%f.png", sys{1}{3}, cond{1}{1}, cond{1}{2}(2));
        saveas(gcf, filename);
    end
end
clf();

%% Part d)
t = 0:0.005:2;
u = zeros(size(t));
for cond = system_conditions
    [ts, y] = ode45(@(t,x) dx(x, @(x) cond{1}{1}, t), [t(1) t(end)], cond{1}{2});
    subplot(3,1,1);
    plot(ts, y(:,1));
    hold on
    subplot(3,1,2);
    plot(ts, y(:,2));
    hold on
    subplot(3,1,3);
    plot(ts, y(:,3));
    hold on

    legend_strings = ["ode45"];
    for sys = linear_systems
        legend_strings(end+1) = sprintf("Linearized (u(t)=%f)", sys{1}{3});
        % convert coordinates to the local linear approx coordintaes
        u(:,:) = cond{1}{1} - sys{1}{3};
        [~, ~, y] = lsim(sys{1}{1}, u, t, cond{1}{2} - sys{1}{2}.');
        y = y + sys{1}{2}.';
        subplot(3,1,1);
        plot(t,y(:,1));
        subplot(3,1,2);
        plot(t,y(:,2));
        subplot(3,1,3);
        plot(t,y(:,3));
    end
    subplot(3,1,1);
    legend("location", "best");
    legend(legend_strings);
    xlabel("Time");
    ylabel("Current in Amps");
    subplot(3,1,2);
    legend("location", "best");
    legend(legend_strings);
    xlabel("Time");
    ylabel("Distance in meters");
    subplot(3,1,3);
    legend("location", "best");
    legend(legend_strings);
    xlabel("Time");
    ylabel("Velocity in m/s");

    plot_title = sprintf("Comparison of lsim to ode45 \n u(t) = %f  x_0 = %f", cond{1}{1}, cond{1}{2}(1));
    sgtitle(plot_title);

    filename = sprintf("images/input_%f_initial_pos_%f_vs_ode45.png", cond{1}{1}, cond{1}{2}(1));
    saveas(gcf, filename);
    subplot(3,1,1);
    hold off
    subplot(3,1,2);
    hold off
    subplot(3,1,3);
    hold off
end
clf();

%% Part f
% for functions see the bottom of the functions section
dts = [.001 .01 .05, .1];
for dt = dts
    t = 0:dt:2;
    for cond = system_conditions
        [ts, y] = ode45(@(t,x) dx(x, @(x) cond{1}{1}, t), [t(1) t(end)], cond{1}{2});
        subplot(3,1,1);
        plot(ts, y(:,1), "r");
        hold on
        subplot(3,1,2);
        plot(ts, y(:,2), "r");
        hold on
        subplot(3,1,3);
        plot(ts, y(:,3), "r");
        hold on

        y = sim_rect(@(t,x) dx(x, @(x) cond{1}{1}, t), t, cond{1}{2});
        subplot(3,1,1);
        plot(t,y(:,1), "g");
        subplot(3,1,2);
        plot(t,y(:,2), "g");
        subplot(3,1,3);
        plot(t,y(:,3), "g");

        y = sim_trap(@(t,x) dx(x, @(x) cond{1}{1}, t), t, cond{1}{2});
        subplot(3,1,1);
        plot(t,y(:,1), "b");
        subplot(3,1,2);
        plot(t,y(:,2), "b");
        subplot(3,1,3);
        plot(t,y(:,3), "b");

        y = sim_runge_kutta(@(t,x) dx(x, @(x) cond{1}{1}, t), t, cond{1}{2});
        subplot(3,1,1);
        plot(t,y(:,1), "m");
        subplot(3,1,2);
        plot(t,y(:,2), "m");
        subplot(3,1,3);
        plot(t,y(:,3), "m");

        plot_title = sprintf("Comparison of numerical methods to ode45 \n dt = %f u(t) = %f  x_0 = %f", dt, cond{1}{1}, cond{1}{2}(1));
        sgtitle(plot_title);

        subplot(3,1,1);
        legend("location", "best");
        legend("ode45", "rectangular", "trapezoidal", "Runge Kutta");
        xlabel("Time");
        ylabel("Current in Amps");
        subplot(3,1,2);
        legend("location", "best");
        legend("ode45", "rectangular", "trapezoidal", "Runge Kutta");
        xlabel("Time");
        ylabel("Distance in meters");
        subplot(3,1,3);
        legend("location", "best");
        legend("ode45", "rectangular", "trapezoidal", "Runge Kutta");
        xlabel("Time");
        ylabel("Velocity in m/s");

        filename = sprintf("images/numerical_approx_dt_%f_input_%f_initial_pos_%f_vs_ode45.png", t(2)-t(1), cond{1}{1}, cond{1}{2}(1));
        saveas(gcf, filename);

        subplot(3,1,1);
        hold off
        subplot(3,1,2);
        hold off
        subplot(3,1,3);
        hold off
    end
end
clf();




%% Functions
function x_dot = dx(x,func_impulse, t)
    global R L m K g
    % dx_1 = (u(t) - R*x_1(t))/(L);
    % dx_2 = x_3(t);
    % dx_3 = m*g - (K*x_1^2(t))/(x_2(t));
    x_dot = [
        (func_impulse(t) - R*x(1))/(L);
        x(3);
        m*g - (K*x(1)^2)/(x(2));
    ];
end

% returns the equilibrium given a constant input function
function eq = equilibrium_constant(u_t)
    global R L m K g
    eq = [
    u_t/R;
    K*u_t^2/(m*g*R^2);
    0;
    ];
end

% Linearize given an constant input function
function [A, B, C, D] = linearize(u_t)
    global R L m K g
    eq = equilibrium_constant(u_t);
    A = [
    -R/L 0 0;
    0 0 1;
    -2*K*eq(1)/(m*eq(2)) (K/m)*(eq(1)/(eq(2)))^2 0;
    ];
    B = [
    1/L;
    0;
    0;
    ];
    C = [
    0 0 0;
    0 1 0;
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
