
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
