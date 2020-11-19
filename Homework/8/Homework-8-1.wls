a = {{1,1,0,1},{0,0,0,1},{1,1,0,0}};
Print["A: ", a];
Print[TeXForm[a]];

at = Transpose[a];
Print["A^T: ", at];
Print[TeXForm[at]];

ata = at.a;
Print["A^T A: ", ata];
Print[TeXForm[ata]];

cpoly = CharacteristicPolynomial[ata, lambda];
Print["Characteristic Polynomial of A^T A: ", cpoly];
Print[TeXForm[cpoly]];

eigenvals = Eigenvalues[ata]
Print["Eigenvaues of A^T A: ", eigenvals];
Print[TeXForm[eigenvals]];

s = SingularValueList[a];
Print["Singular values A: ",s]
Print[TeXForm[s]]
