%% Problem 1
A = [-2 0; 1 -1];
b = [1;0];
u = 5;
[ts, x] = ode45(@(t,x) A*x+b*u, [0 10], [0;0]);
plot(ts, x);
hold on;
y = x*[2; 1];
plot(ts, y, "g-");
fplot(@(x) ((5/2)*(3-exp(-2*x)-2*exp(-x))), [ts(1) ts(end)], "b*");
legend("x1", "x2", "Simulation", "Predicted");
saveas(gcf, "images/p1_predicted_vs_simulation.png");

hold off;
cla

%% Problem 2
R = 1.2;
C = 1;
L = .2;
A = [0 1/C; -1/L -R/L];
b = [0;1/L];

% Part 2
u = 0;
[ts, x] = ode45(@(t,x) A*x+b*u, [0 10], [0;5]);
plot(ts, x(:,1), "b");
hold on;
plot(ts, x(:,2), "r");
fplot(@(x)  1.25*exp(-x)-1.25*exp(-5*x), [ts(1) ts(end)], "b*");
fplot(@(x) -1.25*exp(-x)+6.25*exp(-5*x), [ts(1) ts(end)], "r*");
legend("Sim Voltage Capacitor", "Sim Current", "Predicted Voltage Capacitor", "Predicted Current");
saveas(gcf, "images/p2_2_predicted_vs_simulation.png");
hold off;

% Part 3
A = [0 1;-4  -4];
B = [0;1];
[ts, x] = ode45(@(t,x) A*x+B*exp(-2*t)*sin(t), [0 10], [0;0]);
plot(ts, x*[1; 0], "b");
hold on;
fplot(@(x)  exp(-2*x).*(x-sin(x)), [ts(1) ts(end)], "r*");
legend("Sim output", "Predicted output");
saveas(gcf, "images/p3_predicted_vs_simulation.png");
hold off;
