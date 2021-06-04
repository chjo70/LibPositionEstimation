#include "stdAfx.h"
#include "MatrixConverter.h"

CMatrixConverter::CMatrixConverter(void)
{
}

CMatrixConverter::~CMatrixConverter(void)
{
}

BOOL CMatrixConverter::MatrixDuplicate(double *matrix_in, double *matrix_out, int nRowInMat, int nColInMat)
{
    // assumes matrix is not null.
    /// Matrix ũ�� ����
    /// Matrix ���� ���
    for (int i = 0; i < nRowInMat; ++i) // copy the values
    {
        for (int j = 0; j < nColInMat; ++j)
        {
            *(matrix_out + i * nColInMat + j) = *(matrix_in + i * nColInMat + j);
        }
    }

    return TRUE;
}

BOOL CMatrixConverter::MatrixTranspose(double *matrix_in, double *matrix_out, int nRowInMat, int nColInMat)
{
    /// Matrix ũ�� ����
    /// Matrix Transpose ���
    for (int i = 0; i < nRowInMat; i++)
    {
        for (int j = 0; j < nColInMat; j++) 
        {
            *(matrix_out + j * nRowInMat + i) = *(matrix_in + i * nColInMat + j);
        }
    }

    return TRUE;
}

double CMatrixConverter::MatrixNorm(double *matrix_in, int nColInMat)
{
    double dOutput = 0.0;

    /// Matrix�� Element�� �������� �� 
    for (int i = 0; i < nColInMat; i++)
    {
        dOutput += *(matrix_in + i) * *(matrix_in + i);
    }

    /// ������ ����
    return sqrt(dOutput);
}

// Norm L-oo �޼ҵ� (MatrixExp ���� ���)
double CMatrixConverter::MatrixNormLoo(double *matrix_in, int nRowInMat)
{
    double row_sum;
    double value;

    value = 0.0;

    for (int i = 0; i < nRowInMat; i++ )
    {
        row_sum = 0.0;
        for (int j = 0; j < nRowInMat; j++ )
        {
            row_sum = row_sum + fabs (*(matrix_in + (i + j * nRowInMat)));
        }
        value = max( value, row_sum );
    }

    return value;
}

// Matrix �� ���� ���� ����
BOOL CMatrixConverter::MatrixSquare(double *matrix_in, double *matrix_out, int nRowInMat, int nColInMat)
{
    for (int i = 0; i < nRowInMat; i++)
    {
        for (int j = 0; j < nColInMat; j++)
        {
            *(matrix_out + i * nColInMat + j) = pow(*(matrix_in + i * nColInMat + j), 2.0);
        }
    }

    return TRUE;
}


// Matrix .* ����(Matrix Scale) �޼ҵ� ����
BOOL CMatrixConverter::MatrixScale(double *matrix_in, double *matrix_out, int nRowInMat, int nColInMat, double nScale)
{
    for (int i = 0; i < nRowInMat; i++)
    {
        for (int j = 0; j < nColInMat; j++) 
        {
            *(matrix_out + i * nColInMat + j) = *(matrix_in + i * nColInMat + j) * nScale;
        }
    }
    return TRUE;
}

// Matrix Log2() �޼ҵ� ����
BOOL CMatrixConverter::MatrixLog2(double *matrix_in, double *matrix_out, int nRowInMat, int nColInMat)
{
    for (int i = 0; i < nRowInMat; i++)
    {
        for (int j = 0; j < nColInMat; j++)
        {
            *(matrix_out + i * nColInMat + j) = log( fabs( *(matrix_in + i * nColInMat + j) ) ) / log(2.0);
        }
    }

    return TRUE;
}

BOOL CMatrixConverter::MatrixAdd(double *matrix_fin, double *matrix_lin, double *matrix_out, int nRowInMat, int nColInMat)
{
    /// Matrix ũ�� ����
    /// Matrix ���ϱ� ����
    for (int i = 0; i < nRowInMat; i++)
    {
        for (int j = 0; j < nColInMat; j++) 
        {
            *(matrix_out + i * nColInMat + j) = *(matrix_fin + i * nColInMat + j) + *(matrix_lin + i * nColInMat + j); 
        }
    }

    return TRUE;
}

BOOL CMatrixConverter::MatrixSubstract(double *matrix_fin, double *matrix_lin, double *matrix_out, int nRowInMat, int nColInMat)
{
    /// Matrix ũ�� ����
    /// Matrix ���� ����
    for (int i = 0; i < nRowInMat; i++)
    {
        for (int j = 0; j < nColInMat; j++) 
        {
            *(matrix_out + i * nColInMat + j) = *(matrix_fin + i * nColInMat + j) - *(matrix_lin + i * nColInMat + j); 
        }
    }

    return TRUE;
}


BOOL CMatrixConverter::MatrixMultiply(double *matrix_fin, double *matrix_lin, double *matrix_out, int nRowFInMat, int nColFInMat, int nColLInMat)
{
    /// Matrix ũ�� ����
    /// Matrix ���ϱ� ����
    for (int i = 0; i < nRowFInMat; i++)
    {
        for (int j = 0; j < nColLInMat; j++) 
        {
            *(matrix_out + i * nColLInMat + j) = 0;
            for (int k = 0; k < nColFInMat; k++)
            {
               *(matrix_out + i * nColLInMat + j) += *(matrix_fin + i * nColFInMat + k) * *(matrix_lin + k * nColLInMat + j); 
            }
        }
    }

    return TRUE;
}

BOOL CMatrixConverter::MatrixInverse(double *dpMatrixIn, double *dpMatrixOut, int nRowInMat, int nColInMat)
{
    if (nRowInMat < nColInMat)
    {
        RInverse(dpMatrixIn, dpMatrixOut, nRowInMat, nColInMat);
    }
    else if (nRowInMat > nColInMat)
    {
        LInverse(dpMatrixIn, dpMatrixOut, nRowInMat, nColInMat);
    }
    else    // if (nRowInMat == nColInMat)
    {
        CInverse(dpMatrixIn, dpMatrixOut, nRowInMat);
    }

    return TRUE;
}

/// H' * inv(H * H') * y
BOOL CMatrixConverter::RInverse(double *matrix_in, double *matrix_out, int nRowInMat, int nColInMat)
{
    double *dpLum = NULL;
    double *dpMulCol = NULL;
    double *dpReturn = NULL;

    /// ������ Transposed Matrix ����
    dpLum = new double [nColInMat * nRowInMat];
    dpMulCol = new double [nColInMat * nColInMat];
    dpReturn = new double [nColInMat * nColInMat];

    /// inv(H' * H) * H' * y - except y
    MatrixTranspose(matrix_in, dpLum, nRowInMat, nColInMat);
    MatrixMultiply(matrix_in, dpLum, dpMulCol, nRowInMat, nColInMat, nRowInMat);
    CInverse(dpMulCol, dpReturn, nRowInMat);
    MatrixMultiply(dpLum, dpReturn, matrix_out, nColInMat, nRowInMat, nRowInMat);

    delete [] dpLum;
    delete [] dpMulCol;
    delete [] dpReturn;

    return TRUE;
}

/// inv(H' * H) * H' * y
BOOL CMatrixConverter::LInverse(double *dpMatrixIn, double *dpMatrixOut, int nRowInMat, int nColInMat)
{
    double *dpTrans; // = NULL;
    double *dpMulCol; // = NULL;
    double *dpReturn; //  = NULL;

    /// ������ Transposed Matrix ����
    /// 2���� �迭 �ʱ�ȭ
    dpTrans  = new double [nColInMat * nRowInMat];     // dpTrans[nColInMat][nRowInMat];
    dpMulCol = new double [nColInMat * nColInMat];    // dpMulCol[nColInMat][nColInMat];
    dpReturn = new double [nColInMat * nColInMat];   // dpReturn[nColInMat][nColInMat];

    /// inv(H' * H) * H' * y - except y
    MatrixTranspose(dpMatrixIn, dpTrans, nRowInMat, nColInMat);
    MatrixMultiply(dpTrans, dpMatrixIn, dpMulCol, nColInMat, nRowInMat, nColInMat);
    CInverse(dpMulCol, dpReturn, nColInMat);
    MatrixMultiply(dpReturn, dpTrans, dpMatrixOut, nColInMat, nColInMat, nRowInMat);

    /// Matrix ����
    delete [] dpTrans;
    delete [] dpMulCol;
    delete [] dpReturn;

    return TRUE;
}

/// ������� Inverse Matrix
BOOL CMatrixConverter::CInverse(double *matrix_in, double *matrix_out, int nRowInMat)
{
    int     nRows;      // �Է� ��Ʈ���� Row �� Col ũ��
    int		P[MAX_SVID_NUM];

    double  *pA;        // �Է� ��Ʈ���� ������
    double  *pR;        // ��� ��Ʈ���� ������

    double  MaxPivot;
    double	*pTmpR1;    // For Swap
    double  *pTmpR2;    // For Swap
    double  Temp;       // For Swap
    
    /// Matrix pointer and size copy
    pA = matrix_in;
    pR = matrix_out;
    nRows = nRowInMat;

	// Copy source matrix 'pA' to target matrix 'pR'
    MatrixDuplicate(pA, pR, nRows, nRows);

	/// Gauss reduction - with partial pivoting.
	/// In the permutation vector 'P', vector 'P', element 'i' 
	/// indicates that row 'i' was interchanged with row 'P[i]'.
	/// At the end of the matrix inversion, the rows are unscrambled 
	/// in reverse order.	
	for (int i = 0; i < nRows; i++)
    {
		/// Search for maximum pivot.
		int k = i;
		pTmpR1 = pR + i;

		MaxPivot = (double)fabs(*(pTmpR1 + nRows * i));
		for ( int j = i + 1; j < nRows; j++ )
		{
            double Temp = (double)fabs(*(pTmpR1 + nRows * j));
			if ( Temp >= MaxPivot )
            {
				MaxPivot = Temp;
				k = j;
			}
		}

		/// Maximum pivot is in a different row, so swap rows 'i' & 
		/// 'k' and adjust permutation vector 'P'.
		if ( k != i )
        {
			P[i] = k;
			pTmpR1 = pR + nRows * i;
			pTmpR2 = pR + nRows * k;
			
            for ( int j = 0;j < nRows;j++ )
            {
				Temp = *pTmpR1;
				*pTmpR1++ = *pTmpR2;
				*pTmpR2++ = Temp;
			}
		}
		else
        {
			P[i] = i;
        }

		/// Perform Gaussian Elimination.
		pTmpR1 = pR + nRows * i;
		Temp = *(pTmpR1 + i);
		for ( int j = 0; j < nRows; j++)
        {
			*(pTmpR1 + j) /= Temp;
        }
        *(pTmpR1 + i) = 1.0F / Temp;
		for ( int j = 0; j < nRows; j++ )
		{
			if ( i == j )
            {
				continue;
            }
			pTmpR2 = pR + nRows * j;
			Temp = *(pTmpR2 + i);
			*(pTmpR2 + i) = 0.0;
			for ( int r = 0; r < nRows; r++ )
            {
				*(pTmpR2 + r) -= *(pTmpR1 + r) * Temp;
            }
		}
	}

	///  Unscramble columns in the reverse order of the original permutation.
	for ( int i = nRows - 1; i >= 0; i-- )
	{
		if ( P[i] == i )
        {
			continue;
        }
		pTmpR1 = pR + i;
		pTmpR2 = pR + P[i];
		for ( int j = 0; j < nRows; j++ )
		{
			Temp = *pTmpR1;
			*pTmpR1 = *pTmpR2;
			*pTmpR2 = Temp;
			pTmpR1 += nRows;
			pTmpR2 += nRows;
		}
	}

    return TRUE;
}

/// Least Sqaure �Լ�
BOOL CMatrixConverter::LeastSquare(double *matrix_fin, double *matrix_lin, double *matrix_out, int nRowFInMat, int nColFInMat, int nColLInMat)
{
    double *dpReturn = NULL;

    /// ������ Inversed Matrix ����
    dpReturn = new double [nColFInMat * nRowFInMat];

    /// Least Sqaure ���� �� ��� ���� - inv(X) * Y
    MatrixInverse(matrix_fin, dpReturn, nRowFInMat, nColFInMat);
    MatrixMultiply(dpReturn, matrix_lin, matrix_out, nColFInMat, nRowFInMat, nColLInMat);

    /// Matrix ����
    delete [] dpReturn;

    return TRUE;
}

BOOL CMatrixConverter::HMatrix(double *matrix_fin, double *matrix_lin, double *matrix_out, int nRowFInMat, int nColFInMat)
{
    double *pdVector = NULL;
    double dRetNorm  = 0.0;

    /// ���������� 4�� ���ϸ� FALSE ����
    if (nRowFInMat < 4)
    {
        return FALSE;
    }
    
    /// ������ Temp Vector Matrix ����
    pdVector = new double [nRowFInMat * (nColFInMat + 1)] ;

    /// HMatrix ����
    for (int i = 0; i < nRowFInMat; i++)
    {
        for (int j = 0; j < nColFInMat; j++)
        {
            *(pdVector + i * (nColFInMat + 1) + j) = *(matrix_lin + j) - *(matrix_fin + i * nColFInMat + j);
        }
        
        dRetNorm = MatrixNorm( pdVector + i * (nColFInMat + 1), 3 );

        for (int j = 0; j < nColFInMat + 1; j++)
        {
            *(matrix_out + i * (nColFInMat + 1) + j) = *(pdVector + i * (nColFInMat + 1) + j) / dRetNorm;

            if (j == nColFInMat)
            {
                *(matrix_out + i * (nColFInMat + 1) + j) = 1;
            }
        }
    }

    /// �����迭 ����
    delete [] pdVector;

    return TRUE;
}

// Matrix Exponential �Լ� ���� - Enhanced Kalman Filter���� ���
BOOL CMatrixConverter::MatrixExp(double *matrix_in, double *matrix_out, int nRowInMat, int nRowOutMat) 
{
    int     nMaxResult;
    BOOL    p;
    int     q;
    // double  c = 1.0;
    // double  c = 2.0;
    double  c = 1.0;

    double *pdNormMatrix;
    double *pdScaleMatrix;
    double *pdMultMatrix;
    double *pdIdentMatrixF;
    double *pdIdentMatrixL;

    double dNormResult;
    double dLog2Result;
    double dTimeResult;

    /// Matrix �Է°� ����
    pdNormMatrix  = new double [nRowInMat * nRowInMat];
    pdScaleMatrix = new double [nRowInMat * nRowInMat];
    pdMultMatrix = new double [nRowInMat * nRowInMat];
    pdIdentMatrixF = new double [nRowInMat * nRowInMat];
    pdIdentMatrixL = new double [nRowInMat * nRowInMat];

    /// matrix_in == a, pdNormMatrix == a2
    MatrixDuplicate(matrix_in, pdNormMatrix, nRowInMat, nRowInMat);
    
    /// Matrix Norm ����
    /// pdNormMatrix == a2, dNormResult == a_norm
    dNormResult = MatrixNormLoo(pdNormMatrix, nRowInMat);
    
    /// Norm ����� ���� Matrix Scale ����
    // dNormResult == a_norm, dLog2Result == ee
    // nMaxResult == s, dTimeResult == t
    MatrixLog2(&dNormResult, &dLog2Result, 1, 1);
    nMaxResult = max(0, ((int)dLog2Result + 1) + 1);
    dTimeResult = 1.0 / pow(2.0, nMaxResult);
    // pdScaleMatrix == a2
    MatrixScale(pdNormMatrix, pdScaleMatrix, nRowInMat, nRowInMat, dTimeResult);

    /// ScaleMatrix(a2)�� NormMatrix(x)�� ���� �� ����
    // pdScaleMatrix == a2, pdNormMatrix == x
    MatrixDuplicate(pdScaleMatrix, pdNormMatrix, nRowInMat, nRowInMat);
    
    /// Identity Matrix ���� (e = F, d = L)
    // pdIdentMatrixF == e, pdIdentMatrixL == d
    MatrixIdentity(pdIdentMatrixF, nRowInMat);
    MatrixIdentity(pdIdentMatrixL, nRowInMat);

    /// Matrix Add �ι� ���� (e = F, d = L)
    // pdIdentMatrixF == e, pdScaleMatrix == a2, pdIdentMatrixF == e
    MatrixScale(pdIdentMatrixF, pdIdentMatrixF, nRowInMat, nRowInMat, 1.0);
    // MatrixScale(pdScaleMatrix, pdMultMatrix, nRowInMat, nRowInMat, 0.5);
    MatrixScale(pdScaleMatrix, pdMultMatrix, nRowInMat, nRowInMat, c);
    MatrixAdd(pdIdentMatrixF, pdMultMatrix, pdIdentMatrixF, nRowInMat, nRowInMat);
    // pdIdentMatrixL == d, pdScaleMatrix == a2, pdIdentMatrixL == d
    MatrixScale(pdIdentMatrixL, pdIdentMatrixL, nRowInMat, nRowInMat, 1.0);
    // MatrixScale(pdScaleMatrix, pdMultMatrix, nRowInMat, nRowInMat, -0.5);
    MatrixScale(pdScaleMatrix, pdMultMatrix, nRowInMat, nRowInMat, -c);
    MatrixAdd(pdIdentMatrixL, pdMultMatrix, pdIdentMatrixL, nRowInMat, nRowInMat);

    p = TRUE; q = 6;
    // p = TRUE; q = 12;
    for (int k = 2; k <= q; k++)
    {
        c = c * (double)(q - k + 1) / (double)(k * (2 * q - k + 1));

        // pdMultMatrix == x
        MatrixMultiply(pdScaleMatrix, pdNormMatrix, pdMultMatrix, nRowInMat, nRowInMat, nRowInMat);
        // pdNormMatrix == x
        MatrixDuplicate(pdMultMatrix, pdNormMatrix, nRowInMat, nRowInMat);
        MatrixScale(pdNormMatrix, pdMultMatrix, nRowInMat, nRowInMat, c);
        // pdIdentMatrixF == e, pdMultMatrix == x .* c
        MatrixAdd(pdMultMatrix, pdIdentMatrixF, pdIdentMatrixF, nRowInMat, nRowInMat);

        if (p)
        {
            // pdMultMatrix == x .* c, pdIdentMatrixL == d
            MatrixScale(pdNormMatrix, pdMultMatrix, nRowInMat, nRowInMat, c);
            MatrixAdd(pdMultMatrix, pdIdentMatrixL, pdIdentMatrixL, nRowInMat, nRowInMat);
        }
        else
        {
            // pdMultMatrix == x .* -c, pdIdentMatrixL == d
            MatrixScale(pdNormMatrix, pdMultMatrix, nRowInMat, nRowInMat, -c);
            MatrixAdd(pdMultMatrix, pdIdentMatrixL, pdIdentMatrixL, nRowInMat, nRowInMat);
        }

        if (p == TRUE)
        {
            p = FALSE;
        }
        else
        {
            p = TRUE;
        }
    }

    // (e = F, d = L)
    // r8mat_minvm ( n, n, d, e, e ); // E -> inverse(D) * E
    LeastSquare(pdIdentMatrixL, pdIdentMatrixF, pdScaleMatrix, nRowInMat, nRowInMat, nRowInMat);
    MatrixDuplicate(pdScaleMatrix, pdIdentMatrixF, nRowInMat, nRowInMat);

    // E -> E^(2*S)
    for (int i = 0; i < nMaxResult; i++)
    {
        // pdIdentMatrixF == e
        MatrixMultiply(pdIdentMatrixF, pdIdentMatrixF, pdMultMatrix, nRowInMat, nRowInMat, nRowInMat);
    }
    MatrixDuplicate(pdMultMatrix, matrix_out, nRowInMat, nRowInMat);

    /// �޸� ����
    delete [] pdNormMatrix;
    delete [] pdScaleMatrix;
    delete [] pdMultMatrix;
    delete [] pdIdentMatrixF;
    delete [] pdIdentMatrixL;

    return TRUE;
}

/// (1,1)�� ��û�ϸ� ù��° �� ��ȯ 
double CMatrixConverter::MatrixElement(double *matrix_in, int nRowInMat, int nColInMat , int nReqRowInMat, int nReqColInMat)
{
    int nOffset = ((nRowInMat) * (nReqRowInMat-1) + (nReqColInMat-1));
    return *(double*)(matrix_in + nOffset);    ///double �� �����Ϳ� 1 ���ϸ� 8����Ʈ ������. 
}


///���� ����� �밢 ����� ���� ���Ѵ�. 
double CMatrixConverter::MatrixTrace(double *matrix_in,int nRowInMat)
{
    double dRetVal = 0;
    int uIndex=0;

    ///�밢 ��� ��
    for(uIndex=0; uIndex<nRowInMat; uIndex++ )
    {
        dRetVal += *(matrix_in + uIndex*nRowInMat + uIndex) ;
    }

    return dRetVal;  
}

///���� ����� �밢 ����� ���� ���Ѵ�. ���� ����� �ƴϰų� ��û�� �밢�� ��� ����  �Ѱ踦 �����  0�� ��ȯ �Ѵ�. nDiagonalNum
/// �밢�� ����� start�� end ������ ���� ��ȯ �Ѵ�. 4:4 ����� ��� 1:4, 1:4  �� ��û�ϸ� ��ü �밢 ����� ���� ���Ѵ�. 
double CMatrixConverter::MatrixTrace(double *matrix_in,int nRowInMat, int nColInMat, int nRowFirstNum, int nRowEndNum,int nColFirstNum, int nColEndNum)
{
    double dRetVal = 0;

    if(nRowInMat>= nRowEndNum && nColInMat >=nColEndNum)
    {
        
    
    
    }

    return dRetVal =0;  
}


// 1�� �ʱ�ȭ�� ���
BOOL CMatrixConverter::MatrixOnes(double *dpMatrix, int nRowInMat,int nColInMat)
{
    for (int i = 0; i < nRowInMat; i++)
    {
        for (int j = 0; j < nColInMat; j++)
        {
            *(dpMatrix + i * nColInMat + j) = 1;
        }
    }
    return TRUE;
}

// Identity �޼ҵ� ����
BOOL CMatrixConverter::MatrixIdentity(double *dpMatrix, int nRowInMat)
{
    for (int i = 0; i < nRowInMat; i++)
    {
        for (int j = 0; j < nRowInMat; j++)
        {
            if (i == j)
            {
                *(dpMatrix + i * nRowInMat + j) = 1;
            }
            else
            {
                *(dpMatrix + i * nRowInMat + j) = 0;
            }
        }
    }
    return TRUE;
}

// Identity �޼ҵ� ����(Ư�������� scale)
BOOL CMatrixConverter::MatrixIdentity(double *dpMatrix, int nRowInMat,double dScale)
{
    for (int i = 0; i < nRowInMat; i++)
    {
        for (int j = 0; j < nRowInMat; j++)
        {
            if (i == j)
            {
                *(dpMatrix + i * nRowInMat + j) = dScale;
            }
            else
            {
                *(dpMatrix + i * nRowInMat + j) = 0;
            }
        }
    }
    return TRUE;
}


/// element ���� 
BOOL CMatrixConverter::MatrixReciprocal(double *matrix_in, double *matrix_out, int nRowInMat, int nColInMat)
{
    for (int i = 0; i < nRowInMat; i++)
    {
        for (int j = 0; j < nColInMat; j++) 
        {
            ///0���� ������ �����Ƿ� 
            if(*(matrix_in + i * nColInMat + j) != 0)
            {
                *(matrix_out + i * nColInMat + j) = 1.0 / (*(matrix_in + i * nColInMat + j)) ;
            }
            else
            {
                *(matrix_out + i * nColInMat + j) =0;
            }
        }
    }
    return TRUE;
}


/// (0,0)�� ��û�ϸ� ù��° �� ��ȯ 
double CMatrixConverter::GetMatrixElement(double *matrix_in, int nRowInMat, int nColInMat , int nReqRowInMat, int nReqColInMat)
{
    int nOffset = ((nColInMat) * (nReqRowInMat) + (nReqColInMat));
    return *(double*)(matrix_in + nOffset);    ///double �� �����Ϳ� 1 ���ϸ� 8����Ʈ ������. 
}


double CMatrixConverter::MatrixSum(double *matrix_in,int nRowInMat, int nColInMat)
{
    double dRetVal =0;

    for (int i = 0; i < nRowInMat; i++)
    {
        for (int j = 0; j < nColInMat; j++) 
        {
            dRetVal += *(matrix_in + i * nColInMat + j) ;
        }
    }

    return dRetVal;;
}

///matrix clear
BOOL CMatrixConverter::MatrixClear(double *matrix_in,int nRowInMat, int nColInMat)
{
    memset(matrix_in,0,sizeof(double)*nRowInMat*nColInMat);

    return TRUE;
}

///���� �밢 ��� ���� 
BOOL CMatrixConverter::MakeDiagMatrix(double *matrix_in,double *matrix_out,int nRowInMat)
{
    int uIndex=0;

    MatrixClear(matrix_out,nRowInMat,nRowInMat);

    for(uIndex=0; uIndex<nRowInMat; uIndex++ )
    {
        *(matrix_out + uIndex*nRowInMat + uIndex) = *(matrix_in + uIndex);
    }

    return TRUE;
}

///���� ����� �밢 ��� ����
BOOL CMatrixConverter::GetDiagElementMat(double *matrix_in,double *matrix_out,int nRowInMat)
{
    int uIndex=0;

    for(uIndex=0; uIndex<nRowInMat; uIndex++ )
    {
        *(matrix_out + uIndex) = *(matrix_in + uIndex*nRowInMat + uIndex) ;
    }

    return TRUE;
}

///��Ŀ��� Ư�� ��� ���� ����� ������. 
BOOL CMatrixConverter::GetSubMatrix(double *matrix_in,double *matrix_out, int nRowInMat, int nColInMat, int nRowFirstNum, int nRowEndNum,int nColFirstNum, int nColEndNum)
{

    int nOutputColmSize = 0;

    nOutputColmSize = (nColEndNum - nColFirstNum)+1;

    for (int i = nRowFirstNum; i <= nRowEndNum; i++)
    {
        for (int j = nColFirstNum; j <= nColEndNum; j++)
        {
            *(matrix_out +(i-nRowFirstNum)*nOutputColmSize + (j-nColFirstNum)) = *(matrix_in + i*nColInMat + j);
        }
    }

    return TRUE;
}

///��Ŀ��� Ư�� ��� ����  ���� �Է� subMatrix ������ ������ 
BOOL CMatrixConverter::SetSubMatrix(double *matrix_in,double *matrix_out, int nRowInMat, int nColInMat, int nRowFirstNum, int nRowEndNum,int nColFirstNum, int nColEndNum)
{

    int nOutputColmSize = 0;

    nOutputColmSize = (nColEndNum - nColFirstNum)+1;

    for (int i = nRowFirstNum; i <= nRowEndNum; i++)
    {
        for (int j = nColFirstNum; j <= nColEndNum; j++)
        {
            *(matrix_out + i*nColInMat + j) = *(matrix_in +(i-nRowFirstNum)*nOutputColmSize + (j-nColFirstNum));
        }
    }

    return TRUE;
}