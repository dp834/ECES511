a = {{2,0,0},{0,2,1},{0,1,2},{0,0,0}};
Print["A: ", a];
Print[TeXForm[a]];

at = Transpose[a];
Print["A^T: ", at];
Print[TeXForm[at]];

ata = at.a;
Print["A^T A: ", ata];
Print[TeXForm[ata]];

cpoly = Factor[CharacteristicPolynomial[ata, lambda]];
Print["Characteristic Polynomial of A^T A: ", cpoly];
Print[TeXForm[cpoly]];

eigenvals = Eigenvalues[ata]
Print["Eigenvaues of A^T A: ", eigenvals];
Print[TeXForm[eigenvals]];

s = SingularValueList[a];
Print["Singular values A: ",s]
Print[TeXForm[s]]


{u,sm,vt} = SingularValueDecomposition[a];
Print["svd A u: ",u]
Print[TeXForm[u]]
Print["svd A s: ",sm]
Print[TeXForm[sm]]
Print["svd A vt: ",Transpose[vt]]
Print[TeXForm[Transpose[vt]]]

Print["Inverse A:", Inverse[a]]
Print[TeXForm[Inverse[a]]]
