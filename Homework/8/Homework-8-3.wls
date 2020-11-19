a = {{1,2},{2,1},{3,4},{4,3}};
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

z = a.{{1,-1},{1,1}}
Print["z: ", z];
Print[TeXForm[z]];

{u,sm,vt} = SingularValueDecomposition[a];
Print["svd A u: ",u]
Print[TeXForm[u]]
Print["svd A s: ",sm]
Print[TeXForm[sm]]
Print["svd A vt: ",Transpose[vt]]
Print[TeXForm[Transpose[vt]]]
