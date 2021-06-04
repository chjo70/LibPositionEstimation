//////////////////////////////////////////////////////////////////////////
/*!
 * @file      AnalyticCEP.cpp
 * @brief     
 * @author    Á¶Ã¶Èñ (churlhee.jo@lignex1.com)
 * @date      2013-09-09 ¿ÀÈÄ 4:07 
 * @warning   
 */

#include "pch.h"

#define _USE_MATH_DEFINES

#include <math.h>

#include "../Matrix/Matrix.h"
#include "../UTM/UTM.h"
//#include "../Coordinate/Coordinate.h"

#include "AnalyticCEP.h"


/**
 * @brief		CAnalyticCEP
 * @param		void
 * @return		
 * @author		Á¶Ã¶Èñ (churlhee.jo@lignex1.com)
 * @version		0.0.1
 * @date		2021/03/16 14:36:56
 * @warning		
 */
CAnalyticCEP::CAnalyticCEP(void)
{
}

/**
 * @brief		~CAnalyticCEP
 * @param		void
 * @return		
 * @author		Á¶Ã¶Èñ (churlhee.jo@lignex1.com)
 * @version		0.0.1
 * @date		2021/03/16 14:36:53
 * @warning		
 */
CAnalyticCEP::~CAnalyticCEP(void)
{
}


/**
 * @brief		CalAnalyticNonlinear
 * @param		SELPE_RESULT * pResult
 * @return		bool
 * @author		Á¶Ã¶Èñ (churlhee.jo@lignex1.com)
 * @version		0.0.1
 * @date		2021/03/16 14:36:47
 * @warning		
 */
bool CAnalyticCEP::CalAnalyticNonlinear( SELPE_RESULT *pResult )
{
    int i;
    bool bRet=true;

    double dSigma;

    double dM, dN, dDiv;

    double *pUTMX, *pUTMY;

    pUTMX = m_pUTMX;
    pUTMY = m_pUTMY;

    pResult->dEEP_major_axis = -1;
    pResult->dEEP_minor_axis = -1;
    pResult->dEEP_theta = -1;
    pResult->dCEP_error = -1;

	if (m_nLob > 1) {
		try
		{
			CMatrix theX(2, m_nLob), theInvW(m_nLob, m_nLob);

			for (i = 1; i <= m_nLob; ++i) {
				dM = pResult->dEasting - *pUTMX;
				dN = pResult->dNorthing - *pUTMY;

				dDiv = (dM * dM) + (dN * dN);
				theX(1, i) = dN / dDiv;
				theX(2, i) = dM / dDiv;

				++pUTMX;
				++pUTMY;
			}

			theX.Print();

			//theInvW = ( 1. / sigma_doa * M_PI / 180. ) * theInvW.Ident(m_nLob, m_nLob );
			//theInvW.Print();
			dSigma = DEGREE2RADIAN(_SIGMA_OF_DOA_ * _SIGMA_OF_DOA_);

			CMatrix theQ, theXt;

			theXt = theX.Transpose();
			theQ = theX * theXt;
			theQ.Print();

			theQ = Inv(theQ, &bRet);
			theQ.Print();

			theQ = dSigma * theQ;
			theQ.Print();

			double dPe = -2.0 * log(1 - 0.5);
			double dSigma_x_square, dSigma_y_square, dRho_xy, dA_square, dB_square, dSqaure;

			dSigma_x_square = theQ.get(1, 1);
			dSigma_y_square = theQ.get(2, 2);
			dRho_xy = theQ.get(2, 1);

			dSqaure = (dSigma_x_square - dSigma_y_square);
			dSqaure *= dSqaure;
			dSqaure += (4.0 * dRho_xy * dRho_xy);
			dA_square = 0.5 * (dSigma_x_square + dSigma_y_square + sqrt(dSqaure));
			dB_square = 0.5 * (dSigma_x_square + dSigma_y_square - sqrt(dSqaure));

			pResult->dEEP_theta = RADIAN2DEGREE(0.5 * atan((2. * dRho_xy) / (dSigma_x_square - dSigma_y_square)));
			pResult->dEEP_theta = 90.0 - pResult->dEEP_theta;

			if (dA_square > dB_square) {
				pResult->dEEP_major_axis = sqrt(dA_square * dPe) / 2. / 1000.;
				pResult->dEEP_minor_axis = sqrt(dB_square * dPe) / 2. / 1000.;
			}
			else {
				pResult->dEEP_major_axis = sqrt(dB_square * dPe) / 2. / 1000.;
				pResult->dEEP_minor_axis = sqrt(dA_square * dPe) / 2. / 1000.;
			}

			pResult->dCEP_error = 0.75 * sqrt(dA_square + dB_square) / 2. / 1000.;
		}
		catch (Exception err) {
			bRet = false;
			printf("Error: %s\n", err.msg);
		}
		catch (...) {
			bRet = false;
			printf("An error occured...\n");
		}
	}

    return bRet;

}