// **********************************************************************
// Ultra small LU in C - https://en.wikipedia.org/wiki/LU_decomposition
// Factorize square matrix MAT into LU = MAT, solve x in LU.x = b
// 2006-2026, Franck CHANTELOUP, Institut Laue Langevin, Grenoble, France
// **********************************************************************

typedef double realtype;
//typedef long double realtype;
//typedef float realtype;

// **************************************************************
// Factorize square matrix MAT into LU, store LU in MAT[col][lin]
// **************************************************************
void LU_Doolittle(realtype **MAT, int n)
	for (int i = 0; i < n; i++)
    {   for (int k = i; k < n; k++)
		{   realtype sum = 0.0;
            for (int j = 0; j < i; j++) sum += (MAT[j][i] * MAT[k][j]);
		    MAT[k][i] -= sum; }
        for (int k = i+1; k < n; k++)
		{   realtype sum = 0.0;
            for (int j = 0; j < i; j++) sum += (MAT[j][k] * MAT[i][j]);
	        if (MAT[i][i]!=0.0) MAT[i][k] = (MAT[i][k] - sum) / MAT[i][i]; }
	}
}
// **************************************************************

// **************************************************************
// Solve LU.x = b, store x in b, if unresolve b[i] is set to zero
// **************************************************************
void LU_bSolve(realtype **LU, realtype *b, int n)
{
    for (int i = 0; i < n; i++)
    { realtype sum = 0.0;
      for (int j = 0; j < i; j++) sum += LU[j][i]*b[j] ;
      b[i] -= sum; }
    for (int i = n-1; i>=0; i--)
    { realtype sum = 0.0;
      for (int j = n-1; j > i; j--) sum += LU[j][i]*b[j] ;
      if (LU[i][i]!=0.0) b[i] = (b[i]-sum)/LU[i][i]; else b[i] = 0.0; }
}
// **************************************************************

// **************************************************************
// Solve LU.b = x, store x in b
// **************************************************************
void LU_xSolve(realtype **LU, realtype *b, int n)
{
    for (int i = 0; i < n; i++)
    { realtype sum = 0.0;
      for (int j = i; j < n; j++) sum += LU[j][i]*b[j] ;
      b[i] = sum; }
    for (int i = n-1; i>=0; i--)
    { realtype sum = 0.0;
      for (int j = i-1; j >=0; j--) sum += LU[j][i]*b[j] ;
      b[i] += sum; }
}
// **************************************************************

// *********************************
// Print square matrix MAT, size n*n
// *********************************
void LU_PrintMat(realtype **MAT, int n)
{
    for (int i=0; i < n; i++) {
      for (int j=0; j < n; j++)
        printf("%+8.3e ", MAT[i][j]);
	printf("\n"); }
}
