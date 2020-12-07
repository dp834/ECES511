
%% Problem 1
%ts = 0:.01:5;
%gs = arrayfun(@(x) p1_gt(x), ts);
%ps = arrayfun(@(x) p1_pt(x), ts);
%ys = conv(gs,ps);
%size(ys)
%plot(0:.01:10, ys);
%saveas(gcf, "images/p1_convolution.png");


%% Problem 2
A = [1 2 -3;
     3 5 9 ;
     5 0 3];
v1 = A(:,1);
v2 = A(:,2);
v3 = A(:,3);

% 3 since all columns are linearly independent
rank(A)

% Since we are full rank the nullspace is only the zero vector
null(A)

%% Problem 3
A = [2 1; 6 1; 20 1; 30 1; 40 1];
b = [20;18;10;6;2];
x = pinv(A)*b
myerr = 1/5096 * A*[-2448;105056] - b;
myerr = myerr'*myerr
materr = A*x-b;
materr = materr'*materr

%% Problem 4
% Part 1 and 2
% checking that the state space solutions give the same transfer function
% Since state space solutions aren't unique this is the easiest way to verify
y_s = [6];
r_s = [1 5 6];

[A,B,C,D] = tf2ss(y_s, r_s)

A =  [ 0 1 ; -6 -5];
B =  [0; 6];
C = [1 0];
D = [0];
state_space = ss(A,B,C,D);
transfer_func = tf(state_space)

% Symbolically fine the solution instead of comparing graphs
syms s ;
h = 6/(s^2 + 5*s +6);
impulse_response = ilaplace(h)

%% Problem 5
A = [0 1; -2 3];
Ainv = A^-1
A*Ainv

%% Problem 6
A = [-2 0; 0 -4];
B = [4 ; -1];
C = [1 3];
D = 0;
x0 = [4; 5];

syms t;
expm(A*t)

state_space = ss(A,B,C,D);
transfer_func = tf(state_space)
t = 0:0.05:10;
u = zeros(length(t),1);
u(:,1) = 2;
lsim(state_space,u,t,x0);
hold on;
ylim([ -1 20]);
% Zero state
fplot(@(x)  5/2 + 3/2 *exp(-4*x)-4*exp(-2*x), [t(1) t(end)], "g");
% Zero input 
fplot(@(x) 15*exp(-4*x)-4*exp(-2*x), [t(1) t(end)], "r");
% Predicted
fplot(@(x) 5/2 + 33/2 *exp(-4*x), [t(1) t(end)], "black*");
legend("Lsim", "Zero state response", "Zero input response", "Predicted response");
saveas(gcf, "images/p6_lsim.png");
hold off;

%% Problem 7
A = [-2 0; 0 -4];
B = [4 ; -1];
C = [1 3];
D = 0;
x0 = [4; 5];
ts = 1/2;

state_space = ss(A,B,C,D);
t = 0:0.05:10;
u = zeros(length(t),1);
u(:,1) = 2;

[y, ~, x] = lsim(state_space, u, t, x0);
plot(t, y, "b");
hold on;
plot(t, x(:,1), "g");
plot(t, x(:,2), "r");
legend("y", "x1", "x2");
saveas(gcf, "images/p7_lsim.png");
hold off;

discrete_model = c2d(state_space, ts)
[y, ~, x] = lsim(state_space, u, t, x0);
plot(t, y, "b--");
hold on;
plot(t, x(:,1), "g--");
plot(t, x(:,2), "r--");
t = 0:ts:10
u = zeros(length(t),1);
u(:,1) = 2;
[y, ~, x] = lsim(discrete_model, u, t, x0);
plot(t, y, "b");
plot(t, x(:,1), "g");
plot(t, x(:,2), "r");
legend("y Continuous", "x1 Continuous", "x2 Continuous", "y Discrete", "x1 Discrete", "x2 Discrete");
saveas(gcf, "images/p8_lsim.png");
hold off;

%% Problem 8
A = [2 0 0;
     0 2 1;
     0 1 2;
     0 0 0];
Uc = 1/sqrt(2) * [0 sqrt(2) 0 0;
                  -1   0   -1 0;
                  -1   0    1 0;
                  0    0    0 sqrt(2)];
Wc = [3 0 0;
      0 2 0;
      0 0 1;
      0 0 0];
Vtc = 1/sqrt(2) * [0      -1 -1;
                   sqrt(2) 0  0;
                   0      -1  1];
Uc * Wc * Vtc
ata = A'*A
[U, W, V] = svd(A)
eq(Uc,U)
eq(Wc,W)
eq(Vtc,V)



%% Functions

function y = p1_gt(x)
    if( x < 0 | x > 3 )
        y = 0;
    elseif(x <= 1)
        y = x;
    elseif(x <= 2)
        y = 1;
    elseif(x <= 3)
        y = -(x-3);
    end
end

function y = p1_pt(x)
    if(x >=0 && x <=1)
        y = 1;
    else
        y=0;
    end
end
