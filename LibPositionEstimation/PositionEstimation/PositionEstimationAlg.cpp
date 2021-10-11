//////////////////////////////////////////////////////////////////////////
/*!
 * @file      PositionEstimationAlg.cpp
 * @brief     위치 산출을 계산하기 위한 메인 클래스.
 * @author    조철희 (churlhee.jo@lignex1.com)
 * @date      2013-09-09 오후 3:54
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
 * @brief     위치 산출 생성자 클래스: 초기화 및 메모리 할당한다.
 * @param     void
 * @return    void
 * @exception
 * @author    조철희 (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2017-03-07, 오후 2:38
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
//         sprintf( szBuffer, "\n 위도:%.2f, 경도:%.2f -> X : %.2f, Y : %.2f" , stResult.dBLatitude, stResult.dBLongitude, stResult.dEasting, stResult.dNorthing );
//         TRACE( szBuffer );
//     }

}


/**
 * @brief     ~CPositionEstimationAlg
 * @param     void
 * @return    void
 * @exception
 * @author    조철희 (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2017-03-07, 오후 2:38
 * @warning
 */
CPositionEstimationAlg::~CPositionEstimationAlg(void)
{
	// 소멸자에세 객체를 메모리에서 해지함.

}

/**
 * @brief     클래스의 인스턴스를 해지한다.
 * @param     void
 * @return    void
 * @exception
 * @author    조철희 (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2017-03-07, 오후 2:47
 * @warning
 */
void CPositionEstimationAlg::ReleaseInstance()
{

}

/**
 * @brief     이전의 빔 위치 산출 결과와 LOB 들을 이용하여 위치 산출 결과를 계산한다.
 * @param     SELPositionEstimationResult * pResult
 * @param     SELABTDATA_EXT * pABTExtData
 * @param     std::vector<STR_LOBS>* pVecLOB
 * @return    void
 * @exception
 * @author    조철희 (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2016-03-17, 오후 3:47
 * @warning
 */
void CPositionEstimationAlg::RunPositionEstimation(SELPE_RESULT *pSELPE_RESULT, int nLob, STR_LOBS *pstrLOB)
{
	double *pLatitude, *pLongitude, *pLob;
	double dMinLatitude = 360., dMaxLatitude = -360.;
	double dMinLongitude = 360.0, dMaxLongitude = -360.;

	SELPositionEstimationResult result;

	STR_LOBS *ppVecLOB;

	// 1. 센서 좌표에 대한 메모리 할당
	AllocSensors(nLob);

	// 2. 항공기 위치와 LOB 저장
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

	// 3. 위치 산출
	CommonRunPositionEstimation( pSELPE_RESULT );

	// 4. 메모리 해지
	FreeSensors();

}

/**
 * @brief     VerifyOfPositionEstimation
 * @param     SELPE_RESULT * pResult
 * @param     SELSensorPosition * pSensor
 * @return    bool
 * @exception
 * @author    조철희 (churlhee.jo@lignex1.com)
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

			// 경도가 0 초과일때 방위는 3/4 사분면에 있어야 한다.
			if (dDiffLongitude > 0) {
				if (!(*pSensorLob >= 180 && *pSensorLob <= 360)) {
					pResult->bResult = false;
					pResult->dBLongitude = -1;
					pResult->dBLatitude = -1;
					break;
				}

			}
			// 경도가 0 미만일때 방위는 1/2 사분면에 있어야 한다.
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
 * @brief     위치 산출 결과와 항공가 최종 위치 결과를 비교하여 위치 산출 결과 여부를 최종 결정한다.
 * @param     SELPositionEstimationResult * pResult
 * @param     SELSensorPosition * pSensor
 * @return    void
 * @exception
 * @author    조철희 (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2015-06-23, 오후 4:27
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
 * @author    조철희 (churlhee.jo@lignex1.com)
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

        //초기 값 설정
        ll[0].dLatitude = pSELPE_RESULT->dBLatitude;
        ll[0].dLongitude = pSELPE_RESULT->dBLongitude;

        // 위치 산출한 지점에서 M 을 구한다.
        m_dRefMeasure = CalcAllDeltaTheta( pSELPE_RESULT->dBLongitude, pSELPE_RESULT->dBLatitude, true );

        // 최소점을 찾는다.
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

        // 찾은 것 중에서 제일 작은 값을 취한다.
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
 * @author    조철희 (churlhee.jo@lignex1.com)
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
 * @author    조철희 (churlhee.jo@lignex1.com)
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

//static char g_szDirection[9][10]={ "좌상", "좌측", "좌하", "위  ", "현재", "아래", "우상", "우측", "우하" } ;

/**
 * @brief     SelfCall_2Compensate
 * @param     double dLongitude
 * @param     double dLatitude
 * @return    void
 * @exception
 * @author    조철희 (churlhee.jo@lignex1.com)
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

    // 모든 8 분면에 대해서 계산한다.
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

    // 값을 출력
	//i = 0;
	//i += sprintf_s(szConsole, sizeof(szConsole), "\n#%03d [%.5f/%.5f][%s %.6f,%.6f] %8.4f ", m_TryCompensation, m_dStepLongitude, m_dStepLatitude, g_szDirection[pll->iQuardant], pll->dLatitude, pll->dLongitude, pll->dM);

	//szConsole[i] = 0;
	//::OutputDebugString(szConsole);

    return pll->dM;

}

/**
 * @brief     [deg]/[km] 로 계산한다.
 * @param     double dLongitude
 * @param     double dLatitude
 * @return    double
 * @exception
 * @author    조철희 (churlhee.jo@lignex1.com)
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
		//신뢰성. 파싱 에러가 떠서 주석처리함. 추후 주석 풀것. 윤현철.
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
 * @brief     주 클래스에서 생성한 위치 산출 관련 메모리 할당
 * @param     SELSensorPosition * pSensor
 * @param     int nLob
 * @return    bool
 * @exception
 * @author    조철희 (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2015-06-25, 오후 11:20
 * @warning
 */
void CPositionEstimationAlg::AllocSensors( int nLob )
{
	/*! \debug  신뢰성: LOB 상/한값 설정
			\author 조철희 (churlhee.jo@lignex1.com)
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
 * @brief     주 클래스에서 생성한 위치 산출 관련 메모리 해지
 * @param     void
 * @return    void
 * @exception
 * @author    조철희 (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2015-06-25, 오후 11:21
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
 * @brief     C/E/F 공용 위치 산출 알고리즘을 수행한다. 위치 산출 결과를 예외 처리하여 항공기 및 LOB를 걸러낸다.
 * @param     SELPositionEstimationResult * pResult
 * @param     SELABTDATA_EXT * pABTExtData
 * @return    void
 * @exception
 * @author    조철희 (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2015-06-25, 오후 11:26
 * @warning
 */
void CPositionEstimationAlg::CommonRunPositionEstimation( SELPE_RESULT *pSELPE_RESULT, STR_POSITION_ESTIMATION *pPEInfo  )
{

	// 위치 산출
	// 1. 위치 산출 알고리즘을 계산한다.
	RunPositionEstimation( pSELPE_RESULT, DISTANCE_LEAST_SQUARE, STOP, pPEInfo );

	// 2. 예외처리를 체크하여 위치 산출 결과 여부를 결정한다.
	// 2.1 예외 처리 #1: 항공기 최종 위치 근방에 위치 산출 결과가 나오는지를 검증한다.
	VerifyOfPositionEstimation( pSELPE_RESULT );

	// 2.2 예외 처리 #2: LOB 방향과 위치 산출 결과 방향을 비교해서 위치 산출 결과를 검증한다.
	//VerifyOfLOB( pABTData );

    // 2.2 위치를 이분법을 이용하여 보정 한다.
    CompensationOfPositionEstimation( pSELPE_RESULT );

}

//////////////////////////////////////////////////////////////////////////
/*!
 * @brief     위치 산출 알고리즘을 수행하여 최적이 위치 추정 값을 리턴함.
 * @param     SELPositionEstimationResult * pResult
 * @param     EL_POSIOTN_ESTIMATION_ALGORITHM_OPTION eOption
 * @param     EL_TARGET_STATE_OPTION eTargetState
 * @param     SELABTDATA_EXT * pABTExtData
 * @return    void *
 * @version   0.0.1
 * @exception 예외사항을 입력해주거나 '해당사항 없음' 으로 해주세요.
 * @author    조철희 (churlhee.jo@lignex1.com)
 * @date      2013-09-09 오후 5:08
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

	// 위치 정보 초기화
	//pResult->cep_error = -2.0;

	// LOB 개수가 부족할 때 에러로 리턴한다.
	if( m_nLob >= 2 && pSELPE_RESULT != NULL && IsVerifyLOB() ) { //DTEC_Else
		// 센서와 LOB를 검증하여 위치 산출의 입력 데이터를 걸러낸다.
		// FilteredByCensorPosition();

		// 위치 산출 실행하기 전에 메모리 할당
		bool bResult;
		UINT uiFlagInitPos=0;
		double *dTemp=NULL;

		pValid = m_Sensor.pValid;

		dTargetLocDeg[0] = 0.;
		dEEPData[0] = 0.;

		// 위치 산출의 초기값 지정한다.
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

		// 위치 산출 실행후 메모리 해지
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
 * @author		조철희 (churlhee.jo@lignex1.com)
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
 * @brief     항공기 위치를 고려하여 LOB 유효성을 검증하여 반영해주는 함수
 * @param     void
 * @return    void
 * @exception
 * @author    조철희 (churlhee.jo@lignex1.com)
 * @version   0.0.1
 * @date      2015-06-23, 오후 7:10
 * @warning
 */
#define MAX_OF_CHECK_LOB								(10)
#define THRESHOLD_OF_DISTANCE_OFSENSORS	(50)			// 단위 는 m
void CPositionEstimationAlg::FilteredByCensorPosition()
{
	int i;
	int startLOBIndex;

	bool flag;

	bool *pValid;					// LOB 유효 플레그
	double *pLatitude, *pLongitude;

	SELDISTLOB distlob;

	// 유효 플레그 클리어
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
		// 위경도 좌표 검증
		flag = true;
		if( IS_VALID_LL( *pLatitude, *pLongitude ) && IS_NOT_ZERO_LL( *pLatitude, *pLongitude ) ) {
// 			double *pLatitude2, *pLongitude2;
//
// 			pLatitude2 = pLatitude + 1;
// 			pLongitude2 = pLongitude + 1;
//
// 			/*! \todo   나머지 데이터에 대해서 비교해야 하는데 마킹을 어디에 하기가 어렵다.
// 									현재는 기본적으로 한군데에서 마킹하는 것으로 함.
// 									// for( ; j <= MAX_OF_CHECK_LOB && (i+j) < m_Sensor.n ; ++j ) {
// 			    \author 조철희 (churlhee.jo@lignex1.com)
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
			TRACE( "\n위/경도 에러 : %.3f/%.3f", *pLatitude, *pLongitude );
		}

		// 대상 포인터 이동
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

	// 1차 변환
	iTheta = - iEEPTiltAngle + 900;

	// 2차 변환 (지도에 뿌리기 위해서는 양의 값으로 변환)
	if( iTheta < 0 ) { //EX_DTEC_Else
		iTheta += 7200;
	}

	iTheta %= 1800;

	fTheta = (float) ( iTheta / 10. );

	return fTheta;

}