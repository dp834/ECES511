a = {{0,1},{-2,-2}};
b = {{1},{1}};
c = {{2,3}};

eat = MatrixExp[a t] //FullSimplify;
Print["e^(At): ", eat];
Print[TeXForm[eat]];
toIntegrate = MatrixExp[a (t-tau)].b // FullSimplify;
Print["e^(A(t-tau))B: ", toIntegrate];
Print[TeXForm[toIntegrate]]

int = Integrate[toIntegrate, {tau, 0, t}]
Print["Convolution e^(A(t-tau)) B: ", int];
Print[TeXForm[int]]

out = c.int // FullSimplify
Print["y(t): ", out];
Print[TeXForm[out]]


siaInv = Inverse[s IdentityMatrix[2] - a] //FullSimplify
Print["(sI-A)^-1: ", siaInv]
Print[TeXForm[siaInv]]

g = c.siaInv.b //FullSimplify
Print["C(sI-A)^(-1)B: ", g]
Print[TeXForm[g]]
Print[TeXForm[Apart[g]]]

u = {{1}}
uHat = LaplaceTransform[u, t,s]
Print["Laplace(u): ", uHat]
Print[TeXForm[uHat]]

yHat = g.uHat // FullSimplify
Print["yHat: ", yHat]
Print[TeXForm[yHat]]

y = InverseLaplaceTransform[yHat, s, t]//FullSimplify;
Print["y(t): ", y];
Print[TeXForm[y]]
