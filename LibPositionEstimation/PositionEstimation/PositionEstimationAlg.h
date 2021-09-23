//////////////////////////////////////////////////////////////////////////
/*!
 * @file      PositionEstimationAlg.h
 * @brief     ��ġ ���� �˰���
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @date      2013-09-09 ���� 3:53 
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

// ������Ʈ�� ���� ���� ��ġ ���� ����ġ �̷� ���� ����ü�� �����ؾ� �Ѵ�.
// ELINT�� ��ġ ���� ���� ����ü�� �߰���.
//#include "./ElintRsltDataTypedb.h"
#include "./ELEmitterDataType.h"


#define _POSITION_ESTIMATION_OPTION					(1)


//#define D2R_DEGREE(degree) (degree * M_PI / 180.0)
//#define R2D(radian) (radian * 180.0 / M_PI)



// UTM ��� ��ġ ����
//#define _UTM_POSITION_
//#define _TM_POSITION_

/*!
 * @def				LOB2AOA
 * @brief			LOB ���� AOA ������ ��ȯ���ִ� ��
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @date      2013-09-14 ���� 10:28 
 */
#define LOB2AOA(A)		                    (double) ( ( 360.0 - (double) A ) + 90. )

/**
 * @def       UPPER_LIMIT_OF_LOB
 * @brief     EEP�� ����ϱ� ���� �ּ� ��ġ Ÿ�� ����
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 */
#define UPPER_LIMIT_OF_LOB					(30)						// ���� : ����

/**
 * @def       MAX_OF_LOB
 * @brief     ��ġ �����ϱ� ���� �ִ� LOB ����
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 */
#define MAX_OF_LOB							(1000)

/**
 * @def       ERROR_CHECK_OF_MAX_DISTANCE
 * @brief     �ִ� ��ȿ �Ÿ�
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 */
#define ERROR_CHECK_OF_MAX_DISTANCE			(3000000)          // ���� : [m]

/**
 * @def       ERROR_CHECK_OF_MIN_DISTANCE
 * @brief     �ּ� ��ȿ �Ÿ�
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 */
#define ERROR_CHECK_OF_MIN_DISTANCE			(100)           // ���� : [m]

/**
* [�ĺ��� : CLS-GMU-EL-L-PEA]
*
* [������ : SRS-G-SAFR-012]
*
* @class	CPositionEstimationAlg
* @brief	��ġ �������ִ� Ŭ����
*
* (1) Ŭ���� ����
* - �� Ŭ������ LOB �鿡 ���ؼ� ��ġ �����Ͽ� CED/EEP ����� �������ִ� Ŭ�����̴�.
*
* (2) �����������
* - ��ġ ���� ���̺귯���� Geolocation_CLobsDll.dll �̴�.
*
* (3) ���� �� ����ó��
* - �ش���� ����
*/
class CPositionEstimationAlg
{
private:
    int m_TryCompensation;
	//static CPositionEstimationAlg *m_pInstance;					///< �ν��Ͻ� ��ü

    // ��ġ ���� �˰��� 
	CDistanceLeastSquare m_theDistanceLeastSquare;
	//CQuadratic m_theQuadratric;

    CInverseMethod m_theInverseMethod;

	UINT m_nLob;																				///< ��ġ ������ LOB ����
	SELPE_RESULT m_estEmitterXY;													///< ��ġ ���� ���

	double *m_pLob;																			///< ��ġ ���� ���̺귯���� ����� LOB ������
	SELSensorPosition m_Sensor;													///< ��ġ ���� ���̺귯���� ����� �װ��� ���� ��ǥ

	//CGeoCoordConv m_theGeoCoordConv;

	STR_LOBS *m_pR1;
	STR_LOBS *m_pR2;
	STR_LOBS *m_pR3;
	STR_LOBS *m_pR4;

    double m_dRefMeasure;
    double m_dStepLongitude;
    double m_dStepLatitude;

    // ENUN �϶� �⺻ ��ǥ��
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
 * @brief			�ν��Ͻ� ��ü�� ���´�.
 */
#define ST_PEA	CPositionEstimationAlg::GetInstance()