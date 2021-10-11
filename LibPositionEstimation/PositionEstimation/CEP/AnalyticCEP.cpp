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

    double *pUTMX, *pUTMY;

    double dAoaSquare;

    pUTMX = m_pUTMX;
    pUTMY = m_pUTMY;

    /* Å×½ºÆ® */
//     pUTMX[0] = -100000;
//     pUTMY[0] = 0;
//     pUTMX[1] = 0;
//     pUTMY[1] = 0;
//     pUTMX[2] = -100000;
//     pUTMY[2] = 100000;
//     pResult->dEasting = 33019;
//     pResult->dNorthing = 103360;
// 
//     pResult->dEEP_major_axis = -1;
//     pResult->dEEP_minor_axis = -1;
//     pResult->dEEP_theta = -1;
//     pResult->dCEP_error = -1;

	if (m_nLob > 1 ) {
		try
		{
            double divM, dMHat, dNHat, dU;
            CMatrix theH( m_nLob, 2 );
            CMatrix theQ, theW, theInvW, theTransH;

            dAoaSquare = DEGREE2RADIAN( AOA_VARIANCE );
            dAoaSquare = dAoaSquare * dAoaSquare;

            for( i=1 ; i <= m_nLob ; ++i ) {
                dMHat = pResult->dEasting - *pUTMX;
                dNHat = pResult->dNorthing - *pUTMY;
                dU = dNHat / dMHat;

                divM = ( 1 + dU * dU ) * dMHat;
                theH( i, 1 ) = - dU / divM;
                theH( i, 2 ) = 1. / divM;

                ++pUTMX;
                ++pUTMY;
            }

            // theH.Print();
            
            theW = theW.Ident( m_nLob, m_nLob ) * dAoaSquare;
            theInvW = Inv( theW, & bRet );
            theTransH = theH.Transpose();
            theQ = theTransH * theInvW * theH;
            theQ = Inv( theQ, & bRet );

            //theQ.Print();

			double dC = -2.0 * log(1 - 0.99);
			double dSigma_x_square, dSigma_y_square, dRho_xy, dRanda1, dRanda2, dSqaure, dSigma_xy_square;

 			dSigma_x_square = theQ.get(1, 1);
 			dSigma_y_square = theQ.get(2, 2);
 			dRho_xy = theQ.get(2, 1);

			dSigma_xy_square = (dSigma_x_square + dSigma_y_square);
            dSqaure = (dSigma_x_square - dSigma_y_square);
			dSqaure = ( dSqaure * dSqaure ) + ( 4 * dRho_xy * dRho_xy );
 			dRanda1 = 0.5 * (dSigma_xy_square + sqrt(dSqaure));
			dRanda2 = 0.5 * (dSigma_xy_square - sqrt(dSqaure));

			//char szBuffer[100];
			pResult->dEEP_theta = RADIAN2DEGREE( 0.5 * atan2( (2. * dRho_xy), (dSigma_x_square - dSigma_y_square) ) );
			//sprintf_s( szBuffer, "\n ±â¿ï±â : %.2f", pResult->dEEP_theta );
			//::OutputDebugString(szBuffer);

// 			if (dSigma_x_square - dSigma_y_square < 0) {
// 			    pResult->dEEP_theta = 90.0 - pResult->dEEP_theta;
// 			}
// 			else {
// 				pResult->dEEP_theta = 90.0 - pResult->dEEP_theta + 90.0;
// 			}

			pResult->dEEP_major_axis = sqrt( dRanda1 * dC ) / 1000.;
			pResult->dEEP_minor_axis = sqrt( dRanda2 * dC ) / 1000.;

			pResult->dCEP_error = 0.75 * sqrt( dRanda1 + dRanda2 ) / 2. / 1000.;
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