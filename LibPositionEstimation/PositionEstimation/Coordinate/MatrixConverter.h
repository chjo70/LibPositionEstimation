/********************************************************
 @file      MatrixConverter.h
 @author    ���ؿ� (junwoo.jung@lignex1.com)
 @brief     SBAS ��Ʈ���� ó�� Ŭ����
 @detail    SBAS ��Ʈ���� ó�� Ŭ����
            ��Ʈ���� ����, ����, ����, ����� �� ����
 @version   0.1.0a
 @date      2014/07/23
 @history
  - 2014/07/04 ���ؿ� ������Ʈ ����
  - 2014/07/15 ���ؿ� �⺻���� ��� ����, ����, ����, ����� �ϼ�
  - 2014/07/22 ���ؿ� Least Sqaure ���� �޼ҵ� �߰�
  - 2014/07/23 ���ؿ� Least Sqaure ���꿡�� �޸� ���� ���� ����
  - 2014/07/28 ���ؿ� Norm ���� �޼ҵ� �߰�
  - 2014/07/29 ���ؿ� Enhanced Kalman Filter�� Matrix ���� �߰� ���� (MatrixExp ��)
 @Copyright copyright by LIG Nex1 Co. Ltd
********************************************************/
#ifndef _MATRIX_CONVERTER_H_
#define _MATRIX_CONVERTER_H_

#include <math.h>
#include <malloc.h>

#define MAX_SVID_NUM 3000

class CMatrixConverter
{
///�������
public:

protected:
    // int m_cols;
    // int m_rows;
    // double **m_MatrixOutput;

private:
///����Լ�


public:
    CMatrixConverter(void);
    virtual ~CMatrixConverter(void);

public:

    static BOOL MatrixDuplicate(double *matrix_in, double *matrix_out, int nRowInMat, int nColInMat);
    static BOOL MatrixTranspose(double *matrix_in, double *matrix_out, int nRowInMat, int nColInMat);
    
    static double MatrixNorm(double *matrix_in, int nRowInMat);
    static double MatrixNormLoo(double *matrix_in, int nRowInMat);
    static BOOL   MatrixIdentity(double *dpMatrix, int nRowInMat);
    static BOOL   MatrixIdentity(double *dpMatrix, int nRowInMat,double dScale);/// Identity �޼ҵ� ����(Ư�������� scale)
    
    static BOOL MatrixScale(double *matrix_in, double *matrix_out, int nRowInMat, int nColInMat, double nScale);
    static BOOL MatrixLog2(double *matrix_in, double *matrix_out, int nRowInMat, int nColInMat);
    static BOOL MatrixSquare(double *matrix_in, double *matrix_out, int nRowInMat, int nColInMat);
    
    static BOOL MatrixAdd(double *matrix_fin, double *matrix_lin, double *matrix_out, int nRowInMat, int nColInMat);
    static BOOL MatrixSubstract(double *matrix_fin, double *matrix_lin, double *matrix_out, int nRowInMat, int nColInMat);
    static BOOL MatrixMultiply(double *matrix_fin, double *matrix_lin, double *matrix_out, int nRowFInMAt, int nColFInMat, int nColLInMat);
    
    static BOOL MatrixInverse(double *matrix_in, double *matrix_out, int nRowInMat, int nColInMat);
    static BOOL RInverse(double *matrix_in, double *matrix_out, int nRowInMat, int nColInMat);
    static BOOL LInverse(double *matrix_in, double *matrix_out, int nRowInMat, int nColInMat);
    static BOOL CInverse(double *matrix_in, double *matrix_out, int nRowInMat);

    static BOOL LeastSquare(double *matrix_fin, double *matrix_lin, double *matrix_out, int nRowFInMat, int nColFInMat, int nColLInMat);
    static BOOL HMatrix(double *matrix_fin, double *matrix_lin, double *matrix_out, int nRowFInMat, int nColFInMat);

    static BOOL MatrixExp(double *matrix_in, double *matrix_out, int nRowInMat, int nRowOutMat);

    static double MatrixElement(double *matrix_in, int nRowInMat, int nColInMat , int nReqRowInMat, int nReqColInMat);
    
    static double MatrixTrace(double *matrix_in,int nRowInMat);
    static double MatrixTrace(double *matrix_in,int nRowInMat, int nColInMat, int nRowFirstNum, int nRowEndNum,int nColFirstNum, int nColEndNum);

    static double GetMatrixElement(double *matrix_in, int nRowInMat, int nColInMat , int nReqRowInMat, int nReqColInMat);

    static BOOL MatrixReciprocal(double *matrix_in, double *matrix_out, int nRowInMat, int nColInMat);

    static BOOL MatrixOnes(double *dpMatrix, int nRowInMat,int nColInMat);
    static BOOL MatrixClear(double *matrix_in,int nRowInMat, int nColInMat);
    static double MatrixSum(double *matrix_in,int nRowInMat, int nColInMat);

    static BOOL MakeDiagMatrix(double *matrix_in,double *matrix_out,int nRowInMat);   ///���� �밢 ��� ���� 
    static BOOL GetDiagElementMat(double *matrix_in,double *matrix_out,int nRowInMat);   ///�밢 ��Ҹ� ���� 
    static BOOL GetSubMatrix(double *matrix_in,double *matrix_out, int nRowInMat, int nColInMat, int nRowFirstNum, int nRowEndNum,int nColFirstNum, int nColEndNum);  ///Get sub matrix
    static BOOL SetSubMatrix(double *matrix_in,double *matrix_out, int nRowInMat, int nColInMat, int nRowFirstNum, int nRowEndNum,int nColFirstNum, int nColEndNum);  ///Set sub matrix   

};

#endif // _MATRIX_CONVERTER_H_