//////////////////////////////////////////////////////////////////////////
/*!
 * @file      Quadratic.cpp
 * @brief     Quadratic ��ġ ���� �˰���
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @date      2013-09-09 ���� 4:07 
 * @warning   
 */

#include "pch.h"

#define _USE_MATH_DEFINES

#include <math.h>

#include "DistanceLeastSquare.h"

#include "../UTM/UTM.h"

#include "../Matrix/Matrix.h"



/**
 * @brief		CDistanceLeastSquare
 * @param		void
 * @return		
 * @author		��ö�� (churlhee.jo@lignex1.com)
 * @version		0.0.1
 * @date		2021/03/16 14:37:18
 * @warning		
 */
CDistanceLeastSquare::CDistanceLeastSquare(void)
{
}


/**
 * @brief		~CDistanceLeastSquare
 * @param		void
 * @return		
 * @author		��ö�� (churlhee.jo@lignex1.com)
 * @version		0.0.1
 * @date		2021/03/16 14:37:23
 * @warning		
 */
CDistanceLeastSquare::~CDistanceLeastSquare(void)
{
}

//////////////////////////////////////////////////////////////////////////
/*!
 * @brief     
 * @param     SELPositionEstimationResult * pResult
 * @param     SELUTMTIME * pSensor
 * @param     int nEle
 * @return    void
 * @version   0.0.1
 * @exception ���ܻ����� �Է����ְų� '�ش���� ����' ���� ���ּ���.
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @date      2013-09-09 ���� 7:04 
 * @warning   
 */
bool CDistanceLeastSquare::Run( SELPE_RESULT *pResult, double *pUTMX, double *pUTMY, double *pLob, unsigned int nLob )
{
	int i;
	double a, b, c;
	double A, D, E;
	double B1, B2;
	double dDiv;
	double sumA, sumB, sumB1, sumB2, sumC, sumD, sumE, sumG;

	double dTheta;

	double *ppUTMX, *ppUTMY;
	//double dDistX, dDistY;

	bool bRet;

    m_nLob = nLob;

	if (pResult != NULL) {
		memset(pResult, 0, sizeof(SELPE_RESULT));
		pResult->dEEP_theta = -1;

		// 	for( i=0 ; i < nLob ; ++i ) {
		// 		printf( "\n [%3d] ����[%f], �浵[%f], ����[%f]" , i, pLatitude[i], pLongitude[i], 90.-pLob[i] );
		// 	}

		m_pUTMX = pUTMX;
		ppUTMX = pUTMX;
		ppUTMY = pUTMY;
		m_pUTMY = pUTMY;

		sumA = 0.0;
		sumB = 0.0;
		sumB1 = 0.0;
		sumB2 = 0.0;
		sumC = 0.0;
		sumD = 0.0;
		sumE = 0.0;
		sumG = 0.0;

		// 
		for (i = 0; i < nLob; ++i) {
			dTheta = 90.0 - *pLob;
			dTheta = (dTheta * M_PI) / 180.;
			a = sin(dTheta);
			b = -cos(dTheta);
			c = *pUTMX * a + *pUTMY * b;

			A = a * b;
			B1 = a * a;
			B2 = b * b;
			D = b * c;
			E = a * c;

			sumA = A + sumA;
			sumB = b + sumB;
			sumB1 = B1 + sumB1;
			sumB2 = B2 + sumB2;
			sumC = c + sumC;
			sumD = D + sumD;
			sumE = E + sumE;

			++pUTMX;
			++pUTMY;

			++pLob;
		}

		//sumA2 = sumA * sumA;
		dDiv = (sumB1 * sumB2) - (sumA * sumA);

		pResult->dEEP_major_axis = -1.0;
		pResult->dEEP_minor_axis = -1.0;
		pResult->dEEP_theta = 0.0;
		pResult->dCEP_error = -1.0;

		if (dDiv > 0. || dDiv < 0.) {
			double dDistX, dDistY;

			pResult->bResult = true;

			pResult->dEasting = ((sumB2 * sumE) - (sumA * sumD)) / dDiv;
			pResult->dNorthing = ((sumB1 * sumD) - (sumA * sumE)) / dDiv;

			dDistX = fabs(pResult->dEasting - ppUTMY[0]);
			dDistY = fabs(pResult->dNorthing - ppUTMX[0]);

			if ( /* ( dDistX < 1.00 && dDistY < 1.00 ) || dDistX > 1000000.0 || dDistY > 1000000.0 ) ||  */
				pResult->dEasting == 0 || pResult->dNorthing == 0 /* || pResult->dEasting < 0 */ ) {
				pResult->dEasting = -1;
				pResult->dNorthing = -1;
				pResult->dBLongitude = -1;
				pResult->dBLatitude = -1;
				pResult->dLongitude = -1;
				pResult->dLatitude = -1;
				pResult->bResult = false;
			}
			else {
				CalAnalyticNonlinear(pResult);
			}

		}
		else {
			pResult->bResult = false;

		}

		bRet = pResult->bResult;
	}
	else {
		bRet = false;

	}

	return bRet;


}

#define RADIUS_ERATH			(6378137)		// [m]
/**
 * @brief		CalCEP
 * @param		SELPE_RESULT * pResult
 * @param		SELABTDATA_EXT * pABTExtData
 * @return		bool
 * @author		��ö�� (churlhee.jo@lignex1.com)
 * @version		0.0.1
 * @date		2021/02/18 16:54:40
 * @warning		
 */
bool CDistanceLeastSquare::CalCEP( SELPE_RESULT *pResult, SELABTDATA_EXT *pABTExtData )
{
	bool bRet=true;

	int i, iNumPE;

	double dValue;
	double dEasting, dNorthing;

	double dENorthing, dEEasting;
	double dSSNorthing, dSSEasting, dSLL, dSNorthing, dSEasting;

	double dR1, dR2;
	double dSumSS, dSumA, dC;

	double dDiffNE2, dLamda1, dLamda2, dSLL2;

	if( pABTExtData == NULL || pResult->bResult == false || ( pABTExtData->bFullOfPE == false && pABTExtData->uiPE <= 1 ) ) {
		bRet = false;
	}
	// ��ġ ���� ������ 2�� �̻��� �� CEP�� ����Ѵ�.
	else {
		LatLonToUTMXY( pResult->dBLatitude, pResult->dBLongitude, UTM_ZONE, dEasting, dNorthing );

		if( pABTExtData->bFullOfPE == true ) {
			iNumPE = MAX_OF_LOBS_PE;
		}
		else {
			iNumPE = pABTExtData->uiPE;
		}


		// Northing ��� ���ϱ�
		dENorthing = pResult->dNorthing;
		for( i=0 ; i < iNumPE ; ++i ) {
			dENorthing += pABTExtData->dNorthing[i];
		}
		dENorthing /= (iNumPE+1);
			
		// Easting ��� ���ϱ�
		dEEasting = pResult->dEasting;
		for( i=0 ; i < iNumPE ; ++i ) {
			dEEasting += pABTExtData->dEasting[i];
		}
		dEEasting /= (iNumPE+1);

		// Northing �л� ���ϱ�
		dSSNorthing = dENorthing - pResult->dNorthing;
		dSSNorthing *= dSSNorthing;
		for( i=0 ; i < iNumPE ; ++i ) {
			dValue = dENorthing - pABTExtData->dNorthing[i];
			dSSNorthing += ( dValue * dValue );
		}
		dSSNorthing /= (iNumPE);
		dSNorthing = sqrt( dSSNorthing );			// -> devY2 (dSSNorthing)

		// Easting �л� ���ϱ�
		dSSEasting = dEEasting - pResult->dEasting;
		dSSEasting *= dSSEasting;
		for( i=0 ; i < iNumPE ; ++i ) {
			dValue = dEEasting - pABTExtData->dEasting[i];
			dSSEasting += ( dValue * dValue );
		}
		dSSEasting /= (iNumPE);
		dSEasting = sqrt( dSSEasting );				// -> devX2 ( dSSEasting )

		// Easting/Northing �л� ���ϱ�
		dSLL = ( dENorthing - pResult->dNorthing ) * ( dEEasting - pResult->dEasting );
		for( i=0 ; i < iNumPE ; ++i ) {
			dSLL += ( ( dENorthing - pABTExtData->dNorthing[i] ) * ( dEEasting - pABTExtData->dEasting[i] ) );
		}
		dSLL /= (iNumPE);
		dSLL2 = dSLL * dSLL;

		// CEP ���ϱ�
		dC = -2. * log( 1 - 0.9 );
		dDiffNE2 = ( dSSNorthing - dSSEasting ) * ( dSSNorthing - dSSEasting );
		dLamda1 = ( ( dSSNorthing + dSSEasting ) + sqrt( dDiffNE2 + 4 * dSLL2 ) ) / 2.0;
		dLamda2 = ( ( dSSNorthing + dSSEasting ) - sqrt( dDiffNE2 + 4 * dSLL2 ) ) / 2.0;
		pResult->dCEP_error = RADIUS_ERATH * 0.75 * sqrt( dC * dLamda1 + dC * dLamda2 );
		pResult->dCEP_error = 0.75 * sqrt( dSSNorthing + dSSEasting );

		//Log( enNormal, "CEP=%.2f[m]" , pResult->cep_error );

		// EEP ���� ���ϱ�
		pResult->dEEP_theta = 0.5 * atan( ( 2. * dSNorthing * dSEasting ) / ( dSNorthing - dSEasting ) );
		pResult->dEEP_theta = ( RADIAN2DEGREE( pResult->dEEP_theta ) - 0. ) + 180.;
		pResult->dEEP_theta = fmod( pResult->dEEP_theta + 360.0, 360.0 );

		// EEP ����/���� ���ϱ�
		dSumSS = dSSNorthing + dSSEasting;
		dSumA = ( dSumSS * dSumSS ) - 4.0 * ( ( dSSNorthing * dSSEasting ) - ( dSLL * dSLL ) );
		dR1 = 0.5 * ( dSumSS + sqrt( dSumA ) );
		dR2 = 0.5 * ( dSumSS - sqrt( dSumA ) );

		pResult->dEEP_major_axis = sqrt( dR1 * dC ) * 1.;		
		pResult->dEEP_minor_axis = sqrt( dR2 * dC ) * 1.;		
	}

	return bRet;
}

//////////////////////////////////////////////////////////////////////////
/*!
 * @brief     ���� ����(��ġ ����) ���� ���� 24�� ����
 * @param     SELPositionEstimationResult * pResult
 * @param     SELUTMTIME * pEmitterXY
 * @param     SELUTMTIME * pSensorXY
 * @param     double * pTrueLob
 * @param     int nEle
 * @return    void
 * @version   0.0.1
 * @exception ���ܻ����� �Է����ְų� '�ش���� ����' ���� ���ּ���.
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @date      2014-02-03 ���� 4:42 
 * @warning   
 */
void CDistanceLeastSquare::CalCEP( SELPositionEstimationResult *pResult, SELPE_RESULT *pEmitterXY, SELPE_RESULT *pSensorXY, double *pTrueLob, int nEle )
{
	int i;

	CMatrix H, Ht;
	CMatrix W;
	CMatrix Q;

	double *pLob2;

	bool bret=true;

	// H ���� ����
	pLob2 = pTrueLob;
	H = CMatrix( nEle, 2 );
 	for( i=1 ; i <= nEle ; ++i ) {
 		H( i, 1 ) = sin( *pLob2 );
 		H( i, 2 ) = -cos( *pLob2 );
 
 		++ pLob2;
 	}

	// Ht ���� ����
	pLob2 = pTrueLob;
	Ht = CMatrix( 2, nEle );
	for( i=1 ; i <= nEle ; ++i ) {
		Ht( 1, i ) = sin( *pLob2 );
		Ht( 2, i ) = -cos( *pLob2 );

		++ pLob2;
	}

	// W ���� ����
	pLob2 = pTrueLob;
	W = Diag( nEle );

	double aoaVariance = AOA_VARIANCE * M_PI * M_PI / 180. / 180.;

	//TRACE2( "\nXt:%f, Xt:%f" , pEmitterXY->x, pEmitterXY->y );
		
	for( i=1 ; i <= nEle ; ++i ) {
		double val;
		double diff_x, diff_y;

		diff_x = ( pEmitterXY->dBLatitude - pSensorXY->dBLongitude ) / 1.;
		diff_y = ( pEmitterXY->dBLatitude - pSensorXY->dBLongitude ) / 1.;

		val = 1. / ( ( ( diff_x * diff_x ) + ( diff_y * diff_y ) ) * aoaVariance );
		W( i, i ) = val;

		//TRACE2( "\nXi:%f, Yi:%f" , pSensorXY->x, pSensorXY->y );

		++ pSensorXY;
	}

	// Q���� ����
	/*! \debug  �ŷڼ�: 3���� ���ϱ� ������ 2�ٷ� ������ ����. 
			\author ��ö�� (churlhee.jo@lignex1.com)
			\date 	2015-10-5 16:32:45
	*/
	// Q = Ht * W * H;
	Q = Ht * W;
	Q = Q * H;
	
	Q = Inv( Q, & bret );
	// Q.Print();

	// ��� �ʱ�ȭ
	pResult->cep_error = -1.0;
	pResult->eep_major_axis = -1.0;
	pResult->eep_minor_axis = -1.0;
	pResult->eep_theta = -1.0;

	if( bret == true ) {
		double c;
		double randa_1, randa_2;
		double sigma_x_square2;

		double sigma_x_square = Q( 1, 1 );
		double sigma_y_square = Q( 2, 2 );
		double rho_xy = Q( 1, 2 );
		double rho_xy2 = rho_xy * rho_xy;

		c = -2.0 * log( 1.0 - ( PROBABILITY_OF_BEING_INSIDE / 100. ) );

		sigma_x_square2 = (sigma_x_square - sigma_y_square);
		sigma_x_square2 *= sigma_x_square2;
		randa_1 = 0.5 * ( ( sigma_x_square + sigma_y_square ) + sqrt( sigma_x_square2 + ( 4 * rho_xy2 ) ) );
		randa_2 = 0.5 * ( ( sigma_x_square + sigma_y_square ) - sqrt( sigma_x_square2 + ( 4 * rho_xy2 ) ) );

		/*! \bug  	CEP/EEP ���� -infinity ������ ���� ����
		    \author ��ö�� (churlhee.jo@lignex1.com)
		    \date 	2014-03-22 14:57:22
		*/
		if( randa_1 < 0 || randa_2 < 0 ) {
		}
		else {
			// EEP / CEP ��� ����
			pResult->eep_major_axis = ( sqrt( randa_1 ) * c ) / 1000.;
			pResult->eep_minor_axis = ( sqrt( randa_2 ) * c ) / 1000.;

			pResult->eep_theta = ( 0.5 * atan( ( 2 * rho_xy ) / ( sigma_x_square - sigma_y_square ) ) ) * ( 180. / M_PI );

			/*! \debug  EEP ��ȯ
				\author ��ö�� (churlhee.jo@lignex1.com)
				\date 	2015-07-20 17:01:47
			*/
			pResult->eep_theta = pResult->eep_theta - 90.;

			if( pResult->eep_theta < 0 )
				pResult->eep_theta += 360.;
			else if( pResult->eep_theta >= 360.0 )
				pResult->eep_theta -= 360.;
			else {
			}

			pResult->cep_error = 0.75 * sqrt( randa_1 + randa_2 ) / 1000.;
		}
	}

}

