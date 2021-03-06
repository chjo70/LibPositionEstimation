//////////////////////////////////////////////////////////////////////////
/*!
 * @file      PositionEstimationAlg.h
 * @brief     위치 산출 알고리즘
 * @author    조철희 (churlhee.jo@lignex1.com)
 * @date      2013-09-09 오후 3:53 
 * @warning   
 */

#pragma once

using namespace std;
#include <vector>

using namespace std;

#include "./PositionEstimationDefine.h"

//#include "./Quadratic/Quadratic.h"
#include "./DistanceLeastSquare/DistanceLeastSquare.h"
#include "./LinearLS/LinearLS.h"
#include "./NonlinearLS/NonLinearLS.h"


//#include "./Coordinate/Coordinate.h"

#include "./InverseMethod/CInverseMethod.h"


//#include "../Library/EGeoLocCommonDll.h"

// 프로젝트에 따라서 이전 위치 산출 추정치 이력 관련 구조체를 선언해야 한다.
// ELINT용 위치 산출 정보 구초체를 추가함.
//#include "./ElintRsltDataTypedb.h"
#include "./ELEmitterDataType.h"


#define _POSITION_ESTIMATION_OPTION					(1)


//#define D2R_DEGREE(degree) (degree * M_PI / 180.0)
//#define R2D(radian) (radian * 180.0 / M_PI)



// UTM 계로 위치 산출
//#define _UTM_POSITION_
//#define _TM_POSITION_

/*!
 * @def				LOB2AOA
 * @brief			LOB 값을 AOA 값으로 변환해주는 것
 * @author    조철희 (churlhee.jo@lignex1.com)
 * @date      2013-09-14 오전 10:28 
 */
#define LOB2AOA(A)		                    (double) ( ( 360.0 - (double) A ) + 90. )

/**
 * @def       UPPER_LIMIT_OF_LOB
 * @brief     EEP를 계산하기 위한 최소 위치 타겟 개수
 * @author    조철희 (churlhee.jo@lignex1.com)
 */
#define UPPER_LIMIT_OF_LOB					(30)						// 단위 : 개수

/**
 * @def       MAX_OF_LOB
 * @brief     위치 산출하기 위한 최대 LOB 개수
 * @author    조철희 (churlhee.jo@lignex1.com)
 */
#define MAX_OF_LOB							(1000)

/**
 * @def       ERROR_CHECK_OF_MAX_DISTANCE
 * @brief     최대 유효 거리
 * @author    조철희 (churlhee.jo@lignex1.com)
 */
#define ERROR_CHECK_OF_MAX_DISTANCE			(3000000)          // 단위 : [m]

/**
 * @def       ERROR_CHECK_OF_MIN_DISTANCE
 * @brief     최소 유효 거리
 * @author    조철희 (churlhee.jo@lignex1.com)
 */
#define ERROR_CHECK_OF_MIN_DISTANCE			(100)           // 단위 : [m]

/**
* [식별자 : CLS-GMU-EL-L-PEA]
*
* [추적성 : SRS-G-SAFR-012]
*
* @class	CPositionEstimationAlg
* @brief	위치 산출해주는 클래스
*
* (1) 클래스 설명
* - 본 클래스는 LOB 들에 대해서 위치 산출하여 CED/EEP 결과를 리턴해주는 클래스이다.
*
* (2) 설계결정사항
* - 위치 산출 라이브러리는 Geolocation_CLobsDll.dll 이다.
*
* (3) 제한 및 예외처리
* - 해당사항 없음
*/
class CPositionEstimationAlg
{
private:
    int m_TryCompensation;
	//static CPositionEstimationAlg *m_pInstance;					///< 인스턴스 샛체

    // 위치 산출 알고리즘 
	CDistanceLeastSquare m_theDistanceLeastSquare;
	//CQuadratic m_theQuadratric;

    CInverseMethod m_theInverseMethod;

	UINT m_nLob;																				///< 위치 산출할 LOB 개수
	SELPE_RESULT m_estEmitterXY;													///< 위치 산출 결과

	double *m_pLob;																			///< 위치 산출 라이브러리에 사용할 LOB 데이터
	SELSensorPosition m_Sensor;													///< 위치 산출 라이브러리에 사용할 항공기 센서 좌표

	//CGeoCoordConv m_theGeoCoordConv;

	STR_LOBS *m_pR1;
	STR_LOBS *m_pR2;
	STR_LOBS *m_pR3;
	STR_LOBS *m_pR4;

    double m_dRefMeasure;
    double m_dStepLongitude;
    double m_dStepLatitude;

    // ENUN 일때 기본 좌표계
#if defined(_ENU_POSITION_)
    SLlhPos m_stOrgLlh;
#endif

public:
	CPositionEstimationAlg(void);
	~CPositionEstimationAlg(void);
	//static CPositionEstimationAlg* GetInstance();

	virtual BOOL Finalize() { 
		ReleaseInstance(); 
		return TRUE; 
	};

	static void ReleaseInstance();

	void ConvertLatLong2( unsigned int nLob, SELSensorPosition *pSensor );

    void RunPositionEstimation( SELPE_RESULT *pSELPE_RESULT, int nLob, STR_LOBS *pstrLOB );
	//void RunPositionEstimation( STR_POSITION_ESTIMATION *pPEInfo, SELABTDATA_EXT *pABTExtData, std::vector<STR_LOBS> *pVecLOB );

	void RunPositionEstimation( SELPE_RESULT *pSELPE_RESULT, std::vector<STR_LOBS> *pVecLOB );

	//double EstimatedAltitude( SEnuPos *pstEnuPos );
	//////////////////////////////////////////////////////////////////////////

private:
	void RunPositionEstimation( SELPE_RESULT *pSELPE_RESULT, EL_POSIOTN_ESTIMATION_ALGORITHM_OPTION eOption, EL_TARGET_STATE_OPTION eTargetState=STOP, STR_POSITION_ESTIMATION *pPEInfo=NULL );

	bool VerifyOfPositionEstimation( SELPE_RESULT *pResult, SELSensorPosition *pSensor );
	void VerifyOfPositionEstimation( SELPE_RESULT *pResult );

    void CompensationOfPositionEstimation( SELPE_RESULT *pSELPE_RESULT );
	bool IsCheckMoreDetail(_COMPENSATION_LL_* pLL, int iIndex );
    double SelfCall_2Compensate( _COMPENSATION_LL_ *pll );
    double CalcAllDeltaTheta( double dLongitude, double dLatitude ,bool bInit=false );

	void VerifyOfLOB( SRxABTData *pABTData );

	void FilteredByCensorPosition();

	//bool AllocateBuffer( int isize );
	//void ReleaseBuffer();
	double FuncD2R(double i_dbDegree);
	void AllocSensors( int nLob );
	void CommonRunPositionEstimation( SELPE_RESULT *pSELPE_RESULT, STR_POSITION_ESTIMATION *pPEInfo=NULL );
	void FreeSensors();
	float M2Map( int iEEPTiltAngle );

	bool IsVerifyLOB();

	bool IsEqual(double a, double b);

};

/*!
 * @def				ST_PEA
 * @brief			인스턴스 객체를 얻어온다.
 */
#define ST_PEA	CPositionEstimationAlg::GetInstance()