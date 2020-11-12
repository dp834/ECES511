a = {{0,1},{-2,-2}};
b = {{1},{1}};
c = {{2,3}};

ad = MatrixExp[a T] //FullSimplify
Print["A_d: ", ad]
Print[TeXForm[ad]]


bd = Integrate[MatrixExp[a alpha], {alpha, 0,T}].b //FullSimplify
Print["B_d: ", bd]
Print[TeXForm[bd]]

cd = c;

ad1 = N[ad /. T->1];
Print["A_d(T=1): ", ad1]
Print[TeXForm[ad1]]

bd1 = N[bd /. T->1];
Print["B_d(T=1): ", bd1]
Print[TeXForm[bd1]]

adpi = N[ad /. T->Pi];
Print["A_d(T=pi): ", adpi]
Print[TeXForm[adpi]]

bdpi = N[bd /. T->Pi];
Print["B_d(T=pi): ", bdpi]
Print[TeXForm[bdpi]]
