a = -2;
b = 5
r = 1;
x0 = 10;

int = Integrate[Exp[a (t-tau)] b r,{tau, 0 , t}];
Print["Convolution: " int];
Print[TeXForm[int]]

x[t_] = Exp[a t] x0 + int //FullSimplify
Print["x(t): ", x[t]]
Print[TeXForm[x[t]]]

ClearAll["Global`*"]


xhat = LaplaceTransform[x[t],t,s]
lt = LaplaceTransform[D[x[t],t] == a x[t] + b, t,s]

xhat = xhat /. Solve[lt, xhat]
xt = FullSimplify[InverseLaplaceTransform[xhat /. {a->-2,b->5,x[0]->10}, s,t]];
Print["x(t): ", xt]
Print[TeXForm[xt]]


ad[T_] = Exp[-2 * T];
Print["A_d(T=1/8): ", ad[1/8]];
Print[TeXForm[ad[1/8]]]
Print[TeXForm[N[ad[1/8]]]]

bd[T_] = Integrate[5Exp[-2 alpha],{alpha, 0, T}];
Print["B_d(T=1/8): ", bd[1/8]];
Print[TeXForm[bd[1/8]]]
Print[TeXForm[N[bd[1/8]]]]

dsys[1] = 0;
dsys[n_] := ad[1/8] dsys[n-1] + bd[1/8];
Print[Table[N[dsys[i]], {i,1,10}]]
