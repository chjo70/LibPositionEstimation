//////////////////////////////////////////////////////////////////////////
/*!
 * @file      PositionEstimationAlg.cpp
 * @brief     ��ġ ������ ����ϱ� ���� ���� Ŭ����.
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @date      2013-09-09 ���� 3:54
 * @warning
 */

#include "pch.h"

#include <math.h>
#include <float.h>

#include "IsNumber.h"

#include "PositionEstimationAlg.h"


#include "./UTM/UTM.h"

#include "ELUtil.h"

#include "./InverseMethod/VincentyParam.h"


#define UDIV( A, B )            (UINT) ( (float) (A) / (float) (B) + 0.5 )

#define	IS_VALID_LL( A, B )			        ( ( ( IsNumber(A) == true ) && ( IsNumber(B) == true ) && ( !( ( A > 360. ) || ( A < -360. ) || ( B > 360. ) || ( B < -360. ) ) ) ) == true )
#define	IS_NOT_ZERO_LL( A, B )	( ( ( A > 0 || A < 0 ) == true ) && ( ( B > 0 || B < 0 ) == true ) )

#define MAX_LOBS							(100000)

#define THRESHOLD_OF_LOB_FOR_PE				(20)

#define ITEMS_FOR_TWO_ECLPISE					(5)


#if defined(_ENU_POSITION_)

//#define ORG_ENU_LATITUDE                    DEGREE2RADIAN(37.0)
//#define ORG_ENU_LONGITUDE                   DEGREE2RADIAN(127.0)
#define ORG_ENU_HEIGHT                      (100)


#endif


/**
 * @brief     ��ġ ���� ������ Ŭ����: �ʱ�ȭ �� �޸� �Ҵ��Ѵ�.
 * @param     void
 * @return    void
 * @exception
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2017-03-07, ���� 2:38
 * @warning
 */
CPositionEstimationAlg::CPositionEstimationAlg(void)
{

// 	m_pdCoVar = (double **) malloc( 2 * sizeof(double *) );
// 
// 	if( m_pdCoVar != NULL ) {
// 		*(m_pdCoVar+0) = (double * ) malloc( 2 * sizeof(double) );
// 		*(m_pdCoVar+1) = (double * ) malloc( 2 * sizeof(double) );
// 	}

//#if defined(_ENU_POSITION_)
    //m_stOrgLlh.lat = ORG_ENU_LATITUDE;
    //m_stOrgLlh.lon = ORG_ENU_LONGITUDE;
    //m_stOrgLlh.hgt = ORG_ENU_HEIGHT;

//#endif
//     int i;
//     SELPE_RESULT stResult;
// 
//     char szBuffer[1000];
// 
//     stResult.dBLatitude = 38.0;
//     stResult.dBLongitude = 120.0;
// 
//     TRACE( "S52" );
//     for( i=110 ; i < 140 ; ++i ) {
//         stResult.dBLongitude = (double) i;
//         LatLonToUTMXY( stResult.dBLatitude, stResult.dBLongitude, (int) (52), stResult.dEasting, stResult.dNorthing );
// 
//         sprintf( szBuffer, "\n ����:%.2f, �浵:%.2f -> X : %.2f, Y : %.2f" , stResult.dBLatitude, stResult.dBLongitude, stResult.dEasting, stResult.dNorthing );
//         TRACE( szBuffer );
//     }

}


/**
 * @brief     ~CPositionEstimationAlg
 * @param     void
 * @return    void
 * @exception
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2017-03-07, ���� 2:38
 * @warning
 */
CPositionEstimationAlg::~CPositionEstimationAlg(void)
{
	// �Ҹ��ڿ��� ��ü�� �޸𸮿��� ������.

}

/**
 * @brief     Ŭ������ �ν��Ͻ��� �����Ѵ�.
 * @param     void
 * @return    void
 * @exception
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2017-03-07, ���� 2:47
 * @warning
 */
void CPositionEstimationAlg::ReleaseInstance()
{

}

/**
 * @brief     ������ �� ��ġ ���� ����� LOB ���� �̿��Ͽ� ��ġ ���� ����� ����Ѵ�.
 * @param     SELPositionEstimationResult * pResult
 * @param     SELABTDATA_EXT * pABTExtData
 * @param     std::vector<STR_LOBS>* pVecLOB
 * @return    void
 * @exception
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2016-03-17, ���� 3:47
 * @warning
 */
void CPositionEstimationAlg::RunPositionEstimation(SELPE_RESULT *pSELPE_RESULT, int nLob, STR_LOBS *pstrLOB)
{
	double *pLatitude, *pLongitude, *pLob;
	double dMinLatitude = 360., dMaxLatitude = -360.;
	double dMinLongitude = 360.0, dMaxLongitude = -360.;

	SELPositionEstimationResult result;

	STR_LOBS *ppVecLOB;

	// 1. ���� ��ǥ�� ���� �޸� �Ҵ�
	AllocSensors(nLob);

	// 2. �װ��� ��ġ�� LOB ����
	pLatitude = m_Sensor.pLatitude;
	pLongitude = m_Sensor.pLongitude;
	pLob = m_Sensor.pDOA;
	if (nLob - 1 > 0) {
		ppVecLOB = &pstrLOB[nLob - 1];

		m_pR1 = ppVecLOB;
		m_pR2 = ppVecLOB;
		m_pR3 = ppVecLOB;
		m_pR4 = ppVecLOB;
		for (UINT i = 0; i < m_Sensor.n; ++i) {
			*pLatitude = (double)ppVecLOB->fLatitude;
			*pLongitude = (double)ppVecLOB->fLongitude;
			*pLob = (double)ppVecLOB->fDoa;

#if defined(_ENU_POSITION_)
			if (i >= 1) {
				if (m_pR1->fLatitude > *pLatitude) {
					m_pR1 = ppVecLOB;
				}
				if (m_pR2->fLatitude < *pLatitude) {
					m_pR2 = ppVecLOB;
				}
				if (m_pR3->fLongitude > *pLongitude) {
					m_pR3 = ppVecLOB;
				}
				if (m_pR4->fLongitude < *pLongitude) {
					m_pR4 = ppVecLOB;
				}

			}

			dMinLatitude = min(*pLatitude, dMinLatitude);
			dMaxLatitude = max(*pLatitude, dMinLatitude);

			dMinLongitude = min(*pLongitude, dMinLongitude);
			dMaxLongitude = max(*pLongitude, dMaxLongitude);
#endif

			++pLatitude;
			++pLongitude;
			++pLob;

			--ppVecLOB;

		}

	}
	else {
		;// ppVecLOB = &pstrLOB[0];
	}


#if defined(_ENU_POSITION_)
	if( m_pR1 == m_pR2 ) {
		m_pR2 = m_pR3;
		m_pR3 = m_pR4;
	}
	else if( m_pR1 == m_pR3 ) {
		m_pR3 = m_pR4;
	}
	else if( m_pR2 == m_pR3 ) {
		m_pR3 = m_pR4;
	}
	else {

	}	

	m_stOrgLlh.lat = ( dMinLatitude + dMaxLatitude ) / 2.0;
	m_stOrgLlh.lon = ( dMinLongitude + dMaxLongitude ) / 2.0;
	m_stOrgLlh.hgt = ORG_ENU_HEIGHT;
#endif

	// 3. ��ġ ����
	CommonRunPositionEstimation( pSELPE_RESULT );

	// 4. �޸� ����
	FreeSensors();

}

/**
 * @brief     VerifyOfPositionEstimation
 * @param     SELPE_RESULT * pResult
 * @param     SELSensorPosition * pSensor
 * @return    bool
 * @exception
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2021-04-30, 15:31
 * @warning
 */
bool CPositionEstimationAlg::VerifyOfPositionEstimation( SELPE_RESULT *pResult, SELSensorPosition *pSensor )
{
	UINT uI;
	bool bRet=false;
	double dDiffLatitude, dDiffLongitude;

	double *pSensorLatitude, *pSensorLongitude, *pSensorLob;

	if (pResult != NULL) {
		pSensorLob = pSensor->pDOA;
		pSensorLatitude = pSensor->pLatitude;
		pSensorLongitude = pSensor->pLongitude;
		for (uI = 0; uI < pSensor->n; ++uI) {
			dDiffLatitude = *pSensorLatitude - pResult->dBLatitude;
			dDiffLongitude = *pSensorLongitude - pResult->dBLongitude;

			// �浵�� 0 �ʰ��϶� ������ 3/4 ��и鿡 �־�� �Ѵ�.
			if (dDiffLongitude > 0) {
				if (!(*pSensorLob >= 180 && *pSensorLob <= 360)) {
					pResult->bResult = false;
					pResult->dBLongitude = -1;
					pResult->dBLatitude = -1;
					break;
				}

			}
			// �浵�� 0 �̸��϶� ������ 1/2 ��и鿡 �־�� �Ѵ�.
			else if (dDiffLongitude < 0) {
				if (!(*pSensorLob >= 0 && *pSensorLob <= 180)) {
					pResult->bResult = false;
					pResult->dBLongitude = -1;
					pResult->dBLatitude = -1;
					break;
				}
			}
			else {
				;
			}

			++pSensorLatitude;
			++pSensorLongitude;

			++pSensorLob;

		}

		bRet = pResult->bResult;
	}

	return bRet;

}

/**
 * @brief     ��ġ ���� ����� �װ��� ���� ��ġ ����� ���Ͽ� ��ġ ���� ��� ���θ� ���� �����Ѵ�.
 * @param     SELPositionEstimationResult * pResult
 * @param     SELSensorPosition * pSensor
 * @return    void
 * @exception
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2015-06-23, ���� 4:27
 * @warning
 */
void CPositionEstimationAlg::VerifyOfPositionEstimation( SELPE_RESULT *pResult )
{
	SELDISTLOB distlob;
    double dDist;

	if( ! ( pResult->dCEP_error < 0 ) && m_Sensor.n-1 >= 0 ) {
        dDist = m_theInverseMethod.GCDistance( m_Sensor.pLatitude[m_Sensor.n-1], m_Sensor.pLongitude[m_Sensor.n-1], pResult->dBLatitude, pResult->dBLongitude );

		//ST_IMA->VincentyInverse( & distlob, pResult->latitude, pResult->longitude, m_Sensor.pLatitude[m_Sensor.n-1], m_Sensor.pLongitude[m_Sensor.n-1] );

		//
		if( dDist <= ERROR_CHECK_OF_MIN_DISTANCE || dDist >= ERROR_CHECK_OF_MAX_DISTANCE ) {
            pResult->bResult = false;
			pResult->dCEP_error = -1.0;

		}
	}
}

#define STEP_LONGITUDE			0.1 // (0.0020)
#define STEP_LATITUDE			0.1 // (0.0020)
#define	TRY_OF_COMPENSATION	    (400)
#define THRESHOLD_STEP_LL       (0.000001)
#define THRESHOLD_DELTA_MEASURE (0.00001)

int compDM( const void *pA, const void *pB )
{
    int iRet=0;

    //_COMPENSATION_LL_ *ppA = ( _COMPENSATION_LL_ * ) pA;
    //_COMPENSATION_LL_ *ppB = ( _COMPENSATION_LL_ * ) pB;

	if (((_COMPENSATION_LL_ *)pA)->dM < ((_COMPENSATION_LL_ *)pB)->dM) {
		iRet = -1;
	}
	else if (((_COMPENSATION_LL_ *)pA)->dM > ((_COMPENSATION_LL_ *)pB)->dM) {
		iRet = 1;
	}
    else {

    }

    return iRet;
}

/**
 * @brief     CompensationOfPositionEstimation
 * @param     pSELPE_RESULT
 * @return    void
 * @exception
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2021-04-23, 11:12
 * @warning
 */
void CPositionEstimationAlg::CompensationOfPositionEstimation( SELPE_RESULT *pSELPE_RESULT )
{
    unsigned int i;

    double dM;
    _COMPENSATION_LL_ ll[TRY_OF_COMPENSATION+1];

    if( pSELPE_RESULT != NULL && pSELPE_RESULT->bResult == true ) {
        m_dStepLongitude = STEP_LONGITUDE;
        m_dStepLatitude = STEP_LATITUDE;

        //�ʱ� �� ����
        ll[0].dLatitude = pSELPE_RESULT->dBLatitude;
        ll[0].dLongitude = pSELPE_RESULT->dBLongitude;

        // ��ġ ������ �������� M �� ���Ѵ�.
        m_dRefMeasure = CalcAllDeltaTheta( pSELPE_RESULT->dBLongitude, pSELPE_RESULT->dBLatitude, true );

        // �ּ����� ã�´�.
        m_TryCompensation = 0;
        for( i=0 ; i < (unsigned int) TRY_OF_COMPENSATION ; ++i ) {
            ++ m_TryCompensation;

            dM = SelfCall_2Compensate( & ll[i] );

            if( dM <= THRESHOLD_DELTA_MEASURE || m_dStepLongitude <= THRESHOLD_STEP_LL ) {
                break;
            }

            ll[i+1].dLatitude = ll[i].dLatitude;
            ll[i+1].dLongitude = ll[i].dLongitude;

		    if( IsCheckMoreDetail( ll, i ) ) {
			    m_dStepLongitude *= 0.5;
			    m_dStepLatitude *= 0.5;
		    }
        }

        // ã�� �� �߿��� ���� ���� ���� ���Ѵ�.
        qsort( ll, i+1, sizeof(_COMPENSATION_LL_), compDM  );

        pSELPE_RESULT->dLongitude = ll[0].dLongitude;
        pSELPE_RESULT->dLatitude = ll[0].dLatitude;

    }

}

/**
 * @brief     IsCheckMoreDetail
 * @param     _COMPENSATION_LL_ * pLL
 * @param     int iIndex
 * @return    bool
 * @exception
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2021-04-29, 10:53
 * @warning
 */
bool CPositionEstimationAlg::IsCheckMoreDetail(_COMPENSATION_LL_ *pLL, int iIndex )
{
	bool bRet = false;

	if ( iIndex >= 3 ) {
		if( pLL[iIndex - 2].iQuardant == pLL[iIndex].iQuardant && pLL[iIndex - 3].iQuardant == pLL[iIndex - 1].iQuardant && 
			IsEqual( pLL[iIndex].dM, pLL[iIndex - 2].dM ) && IsEqual( pLL[iIndex-1].dM, pLL[iIndex - 3].dM ) ) {
			bRet = true;
		}
	}

	return bRet;
}

/**
 * @brief     IsEqual
 * @param     double a
 * @param     double b
 * @return    bool
 * @exception
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2021-04-29, 10:53
 * @warning
 */
bool CPositionEstimationAlg::IsEqual(double a, double b)
{
	bool bRet;

	if (a > b || a < b ) {
		bRet = false;
	}
	else {
		bRet = true;
	}

	return bRet;
}

//static char g_szDirection[9][10]={ "�»�", "����", "����", "��  ", "����", "�Ʒ�", "���", "����", "����" } ;

/**
 * @brief     SelfCall_2Compensate
 * @param     double dLongitude
 * @param     double dLatitude
 * @return    void
 * @exception
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2021-04-23, 11:52
 * @warning
 */
double CPositionEstimationAlg::SelfCall_2Compensate( _COMPENSATION_LL_ *pll )
{
    int i, j;
    //char szConsole[300];
    
    double dMinDm, dM;

    double dLongitude=pll->dLongitude, dLatitude=pll->dLatitude;

    // ��� 8 �и鿡 ���ؼ� ����Ѵ�.
    dMinDm = DBL_MAX;
    for( i=-1 ; i <= 1 ; ++i ) {
        for( j=-1 ; j <= 1 ; ++j ) {
            dM = CalcAllDeltaTheta( dLongitude + ( i * m_dStepLongitude ), dLatitude + ( j * m_dStepLatitude ) );
            if( dMinDm > dM ) {
                pll->dLongitude = dLongitude + ( i * m_dStepLongitude );
                pll->dLatitude = dLatitude + ( j * m_dStepLatitude );

                pll->iQuardant = ( 3 * i ) + j + 4;
                pll->dM = dM;
                dMinDm = dM;
            }
        }
    }

    // ���� ���
	//i = 0;
	//i += sprintf_s(szConsole, sizeof(szConsole), "\n#%03d [%.5f/%.5f][%s %.6f,%.6f] %8.4f ", m_TryCompensation, m_dStepLongitude, m_dStepLatitude, g_szDirection[pll->iQuardant], pll->dLatitude, pll->dLongitude, pll->dM);

	//szConsole[i] = 0;
	//::OutputDebugString(szConsole);

    return pll->dM;

}

/**
 * @brief     [deg]/[km] �� ����Ѵ�.
 * @param     double dLongitude
 * @param     double dLatitude
 * @return    double
 * @exception
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2021-04-23, 11:30
 * @warning
 */
double CPositionEstimationAlg::CalcAllDeltaTheta( double dLongitude, double dLatitude, bool bInit )
{
    UINT i;
    double dM=0.0;

    double *pLatitude, *pLongitude, *pDOA;

    pDOA = m_Sensor.pDOA;
    pLongitude = m_Sensor.pLongitude;
    pLatitude = m_Sensor.pLatitude;
    for( i=0 ; i < m_Sensor.n ; ++i ) {
        double dDTheta, dDiv;

        dDTheta = m_theInverseMethod.GCAzimuth( *pLatitude, *pLongitude, dLatitude, dLongitude, true );
        dDiv = ( *pDOA - dDTheta );
        if( dDiv < -180.0 ) {
            dDiv = 360.0 + dDiv;
        }
		//�ŷڼ�. �Ľ� ������ ���� �ּ�ó����. ���� �ּ� Ǯ��. ����ö.
		dDiv = sin(FuncD2R(dDiv)) * m_theInverseMethod.GCDistance( *pLatitude, *pLongitude, dLatitude, dLongitude );
        dM += ( dDiv * dDiv );
        ++ pLongitude;
        ++ pLatitude;
        ++ pDOA;
    }

    if( bInit == false ) 
	{
        dM /= m_dRefMeasure;
    }
    return dM;
}
double CPositionEstimationAlg::FuncD2R(double i_dbDegree)
{
	return i_dbDegree * M_PI / 180.0;
}
/**
 * @brief     �� Ŭ�������� ������ ��ġ ���� ���� �޸� �Ҵ�
 * @param     SELSensorPosition * pSensor
 * @param     int nLob
 * @return    bool
 * @exception
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2015-06-25, ���� 11:20
 * @warning
 */
void CPositionEstimationAlg::AllocSensors( int nLob )
{
	/*! \debug  �ŷڼ�: LOB ��/�Ѱ� ����
			\author ��ö�� (churlhee.jo@lignex1.com)
			\date 	2015-10-6 9:54:05
	*/
	m_Sensor.n = (UINT) min( nLob, MAX_OF_LOBS_PE );

    if( nLob != 0 ) {
	    m_Sensor.pLatitude = new double[m_Sensor.n];
	    m_Sensor.pLongitude = new double[m_Sensor.n];
	    m_Sensor.pAltitude = new double[m_Sensor.n];

	    m_Sensor.pX = new double[m_Sensor.n];
	    m_Sensor.pY = new double[m_Sensor.n];
        m_Sensor.pH = new double[m_Sensor.n];

	    m_Sensor.pDOA = new double[m_Sensor.n];
	    m_Sensor.pTime = new time_t[m_Sensor.n];
	    m_Sensor.pValid = new bool[m_Sensor.n];
    }

}

/**
 * @brief     �� Ŭ�������� ������ ��ġ ���� ���� �޸� ����
 * @param     void
 * @return    void
 * @exception
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2015-06-25, ���� 11:21
 * @warning
 */
void CPositionEstimationAlg::FreeSensors()
{
    if( m_Sensor.n != 0 ) {
	    delete[] m_Sensor.pLatitude;
	    delete[] m_Sensor.pLongitude;
	    delete[] m_Sensor.pAltitude;

	    delete[] m_Sensor.pX;
	    delete[] m_Sensor.pY;
        delete[] m_Sensor.pH;

	    delete[] m_Sensor.pDOA;
	    delete[] m_Sensor.pTime;
	    delete[] m_Sensor.pValid;
    }

}

/**
 * @brief     C/E/F ���� ��ġ ���� �˰����� �����Ѵ�. ��ġ ���� ����� ���� ó���Ͽ� �װ��� �� LOB�� �ɷ�����.
 * @param     SELPositionEstimationResult * pResult
 * @param     SELABTDATA_EXT * pABTExtData
 * @return    void
 * @exception
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2015-06-25, ���� 11:26
 * @warning
 */
void CPositionEstimationAlg::CommonRunPositionEstimation( SELPE_RESULT *pSELPE_RESULT, STR_POSITION_ESTIMATION *pPEInfo  )
{

	// ��ġ ����
	// 1. ��ġ ���� �˰����� ����Ѵ�.
	RunPositionEstimation( pSELPE_RESULT, DISTANCE_LEAST_SQUARE, STOP, pPEInfo );

	// 2. ����ó���� üũ�Ͽ� ��ġ ���� ��� ���θ� �����Ѵ�.
	// 2.1 ���� ó�� #1: �װ��� ���� ��ġ �ٹ濡 ��ġ ���� ����� ���������� �����Ѵ�.
	VerifyOfPositionEstimation( pSELPE_RESULT );

	// 2.2 ���� ó�� #2: LOB ����� ��ġ ���� ��� ������ ���ؼ� ��ġ ���� ����� �����Ѵ�.
	//VerifyOfLOB( pABTData );

    // 2.2 ��ġ�� �̺й��� �̿��Ͽ� ���� �Ѵ�.
    CompensationOfPositionEstimation( pSELPE_RESULT );

}

//////////////////////////////////////////////////////////////////////////
/*!
 * @brief     ��ġ ���� �˰����� �����Ͽ� ������ ��ġ ���� ���� ������.
 * @param     SELPositionEstimationResult * pResult
 * @param     EL_POSIOTN_ESTIMATION_ALGORITHM_OPTION eOption
 * @param     EL_TARGET_STATE_OPTION eTargetState
 * @param     SELABTDATA_EXT * pABTExtData
 * @return    void *
 * @version   0.0.1
 * @exception ���ܻ����� �Է����ְų� '�ش���� ����' ���� ���ּ���.
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @date      2013-09-09 ���� 5:08
 * @warning
 */
void CPositionEstimationAlg::RunPositionEstimation( SELPE_RESULT *pSELPE_RESULT, EL_POSIOTN_ESTIMATION_ALGORITHM_OPTION eOption, EL_TARGET_STATE_OPTION eTargetState, STR_POSITION_ESTIMATION *pPEInfo )
{
	time_t firstToa=0;

	bool *pValid;

	double dTargetLocDeg[2];
	double dInitPosRad[2];
	double dEEPData[4];

	m_nLob = m_Sensor.n;

	// ��ġ ���� �ʱ�ȭ
	//pResult->cep_error = -2.0;

	// LOB ������ ������ �� ������ �����Ѵ�.
	if( m_nLob >= 2 && pSELPE_RESULT != NULL && IsVerifyLOB() ) { //DTEC_Else
		// ������ LOB�� �����Ͽ� ��ġ ������ �Է� �����͸� �ɷ�����.
		// FilteredByCensorPosition();

		// ��ġ ���� �����ϱ� ���� �޸� �Ҵ�
		bool bResult;
		UINT uiFlagInitPos=0;
		double *dTemp=NULL;

		pValid = m_Sensor.pValid;

		dTargetLocDeg[0] = 0.;
		dEEPData[0] = 0.;

		// ��ġ ������ �ʱⰪ �����Ѵ�.
		if( pPEInfo != nullptr && pPEInfo->enValid == E_EL_PESTAT_SUCCESS ) {
			uiFlagInitPos = 1;

			dInitPosRad[0] = (double) pPEInfo->fLatitude;
			dInitPosRad[1] = (double) pPEInfo->fLongitude;
		}
		else {
			uiFlagInitPos = 0;

			dInitPosRad[0] = 0.;
			dInitPosRad[1] = 0.;
		}


		ConvertLatLong2( m_nLob, & m_Sensor );

		bResult = m_theDistanceLeastSquare.Run( pSELPE_RESULT, m_Sensor.pX, m_Sensor.pY, m_Sensor.pDOA, m_nLob );        

		if( bResult == true ) {
#if defined(_UTM_POSITION_)
            UTMXYToLatLon( pSELPE_RESULT->dEasting, pSELPE_RESULT->dNorthing, (int) UTM_ZONE, false, pSELPE_RESULT->dBLatitude, pSELPE_RESULT->dBLongitude );
            pSELPE_RESULT->dBLatitude = RadToDeg( pSELPE_RESULT->dBLatitude );
            pSELPE_RESULT->dBLongitude = RadToDeg( pSELPE_RESULT->dBLongitude );
			pSELPE_RESULT->dAltitude = 0.0;

#elif defined(_ENU_POSITION_)
            SLlhPos stLlhPos;
            SEnuPos stEnuPos;

            stEnuPos.east = pSELPE_RESULT->dEasting;
            stEnuPos.north = pSELPE_RESULT->dNorthing;
            stEnuPos.up = EstimatedAltitude( & stEnuPos );

            CCoordinate::ConvertENU2LLH( stEnuPos, m_stOrgLlh, & stLlhPos );
            pSELPE_RESULT->dBLatitude = stLlhPos.lat*RAD2DEG;
            pSELPE_RESULT->dBLongitude = stLlhPos.lon*RAD2DEG;	
			pSELPE_RESULT->dAltitude = stLlhPos.hgt;


#elif defined(_TM_POSITION_)
            ////m_theGeoCoordConv.SetSrcType( kWgs84, kTmWest );
            ////m_theGeoCoordConv.SetDstType( kWgs84, kGeographic );
            ////m_theGeoCoordConv.Conv( pSELPE_RESULT->dEasting, pSELPE_RESULT->dNorthing, pSELPE_RESULT->dBLongitude, pSELPE_RESULT->dBLatitude );
			pSELPE_RESULT->dAltitude = 0.0;
            //VerifyOfPositionEstimation( pResult, & m_Sensor );
#else

#endif

		}

		// ��ġ ���� ������ �޸� ����
		//ReleaseBuffer();
	}
	else {
		if (pSELPE_RESULT != NULL) {
			pSELPE_RESULT->bResult = false;
		}
		//pABTData->iPEValid = _spZero;
		//pABTData->fCEP = -1.0;
	}

	return;
}

/**
 * @brief		IsVerifyLOB
 * @return		bool
 * @author		��ö�� (churlhee.jo@lignex1.com)
 * @version		0.0.1
 * @date		2021/02/18 16:20:12
 * @warning		
 */
bool CPositionEstimationAlg::IsVerifyLOB()
{
	int i;
	bool bRet=false;

	double dLatitude, dLongitude;
	double *pLatitude, *pLongitude;

	pLatitude = & m_Sensor.pLatitude[0];
	pLongitude = & m_Sensor.pLongitude[0];
	dLatitude = *pLatitude;
	dLongitude = *pLongitude;
	++ pLatitude;
	++ pLongitude;
	for( i=1 ; i < (int) m_Sensor.n ; ++i ) {
		//if( *pLatitude != dLatitude || *pLongitude != dLongitude ) {
		if ( IsEqual( *pLatitude, dLatitude ) == false || IsEqual( *pLongitude, dLongitude) == false ) {
			bRet = true;
			break;
		}
	
		++ pLatitude;
		++ pLongitude;

	}

	return bRet;
}

void CPositionEstimationAlg::ConvertLatLong2( unsigned int nLob, SELSensorPosition *pSensor )
{
	unsigned int i;
	double *pX, *pY, *pH;
	double *pLat, *pLong;


	pX = pSensor->pX;
	pY = pSensor->pY;
    pH = pSensor->pH;
	pLat = pSensor->pLatitude;
	pLong = pSensor->pLongitude;
	for( i=0 ; i < nLob ; ++i ) {
#if defined(_UTM_POSITION_)
		LatLonToUTMXY( *pLat, *pLong, (int) UTM_ZONE, *pX, *pY );
#elif defined(_ENU_POSITION_)
        SLlhPos stLlhPos;
        SEnuPos stEnuPos;

        stLlhPos.lat = DEGREE2RADIAN( *pLat );
        stLlhPos.lon = DEGREE2RADIAN( *pLong );
        stLlhPos.hgt = 100;
        CCoordinate::ConvertLLH2ENU( stLlhPos, m_stOrgLlh, & stEnuPos );

        *pX = stEnuPos.east;
        *pY = stEnuPos.north;
        *pH = stEnuPos.up;

#elif defined(_TM_POSITION_)
		//m_theGeoCoordConv.SetSrcType(kWgs84, kGeographic );
		//m_theGeoCoordConv.SetDstType(kWgs84, kTmWest );
		//m_theGeoCoordConv.Conv( *pLong, *pLat, *pX, *pY );
#else

#endif

		++ pX;
		++ pY;
        ++ pH;

		++ pLat;
		++ pLong;
		
	}
}

/**
 * @brief     �װ��� ��ġ�� ����Ͽ� LOB ��ȿ���� �����Ͽ� �ݿ����ִ� �Լ�
 * @param     void
 * @return    void
 * @exception
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2015-06-23, ���� 7:10
 * @warning
 */
#define MAX_OF_CHECK_LOB								(10)
#define THRESHOLD_OF_DISTANCE_OFSENSORS	(50)			// ���� �� m
void CPositionEstimationAlg::FilteredByCensorPosition()
{
	int i;
	int startLOBIndex;

	bool flag;

	bool *pValid;					// LOB ��ȿ �÷���
	double *pLatitude, *pLongitude;

	SELDISTLOB distlob;

	// ��ȿ �÷��� Ŭ����
	memset( m_Sensor.pValid, 0, m_nLob * sizeof(bool) );

	if( m_Sensor.n <= MAX_OF_LOB ) {
		startLOBIndex = 0;
	}
	else {
		startLOBIndex = 0;
		startLOBIndex = (int) m_Sensor.n - MAX_OF_LOB;
	}

	pValid = & m_Sensor.pValid[startLOBIndex];
	pLatitude = & m_Sensor.pLatitude[startLOBIndex];
	pLongitude = & m_Sensor.pLongitude[startLOBIndex];
	for( i=startLOBIndex ; i < (int) m_Sensor.n-1 ; ++i ) {
		// ���浵 ��ǥ ����
		flag = true;
		if( IS_VALID_LL( *pLatitude, *pLongitude ) && IS_NOT_ZERO_LL( *pLatitude, *pLongitude ) ) {
// 			double *pLatitude2, *pLongitude2;
//
// 			pLatitude2 = pLatitude + 1;
// 			pLongitude2 = pLongitude + 1;
//
// 			/*! \todo   ������ �����Ϳ� ���ؼ� ���ؾ� �ϴµ� ��ŷ�� ��� �ϱⰡ ��ƴ�.
// 									����� �⺻������ �ѱ������� ��ŷ�ϴ� ������ ��.
// 									// for( ; j <= MAX_OF_CHECK_LOB && (i+j) < m_Sensor.n ; ++j ) {
// 			    \author ��ö�� (churlhee.jo@lignex1.com)
// 			    \date 	2015-07-20 17:39:25
// 			*/
// 			if( IS_VALID_LL( *pLatitude2, *pLongitude2 ) ) {
// 				ST_IMA->VincentyInverse( & distlob, *pLatitude, *pLongitude, *pLatitude2, *pLongitude2 );
// 				if( distlob.distance < THRESHOLD_OF_DISTANCE_OFSENSORS ) { //DTEC_Else
// 					flag = false;
// 				}
// 			}

			*pValid = flag;
		}
		else { //DTEC_Else
			*pValid = false;
			TRACE( "\n��/�浵 ���� : %.3f/%.3f", *pLatitude, *pLongitude );
		}

		// ��� ������ �̵�
		++ pValid;
		++ pLatitude;
		++ pLongitude;

	}

	if( IS_VALID_LL( m_Sensor.pLatitude[m_Sensor.n-1], m_Sensor.pLongitude[m_Sensor.n-1] ) && IS_NOT_ZERO_LL( *pLatitude, *pLongitude ) ) {
		m_Sensor.pValid[ m_Sensor.n-1 ] = true;
	}

}

float CPositionEstimationAlg::M2Map( int iEEPTiltAngle )
{
	int iTheta=0;
	float fTheta=0.0;

	// 1�� ��ȯ
	iTheta = - iEEPTiltAngle + 900;

	// 2�� ��ȯ (������ �Ѹ��� ���ؼ��� ���� ������ ��ȯ)
	if( iTheta < 0 ) { //EX_DTEC_Else
		iTheta += 7200;
	}

	iTheta %= 1800;

	fTheta = (float) ( iTheta / 10. );

	return fTheta;

}