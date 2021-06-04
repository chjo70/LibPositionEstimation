//////////////////////////////////////////////////////////////////////////
/*!
 * @file      PositionEstimation.h
 * @brief     ��ġ ���� �˰��� ����ü ����
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @date      2013-09-09 ���� 3:53 
 * @warning   
 */

#pragma once

// ��Ž �л갪
// Analytic CEP�� ����ϱ� ���� ������
#define AOA_VARIANCE									(1)						// ���� : degree

// CEP/EEP Ȯ�� ����
#define PROBABILITY_OF_BEING_INSIDE		(50)					// ���� : %

/*!
 * @typedef   SELPositionEstimation
 * @brief			��ġ ���� ����� ������ ����ü ����
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @date      2013-09-09 ���� 4:54
 */
struct SELPositionEstimationResult {
	// ��ġ ���� ���� ���
	double latitude;			// ����
	double longitude;			// �浵

	// CEP
	double cep_error;						// ���� : m

	// EEP ���� ����
	// ����
	double eep_major_axis;			// ���� : m
	// ����
	double eep_minor_axis;			// ���� : m
	// theta
	double eep_theta;						// ���� : degree

	// EOB �Ǵ� �ٸ� ��ġ ���� �˰��򿡼� ���� ��ġ ���� �ʱⰪ
	double reallatitude;					
	double reallongitude;
	
	SELPositionEstimationResult()
	{
		latitude=0.0;			// ����
		longitude=0.0;			// �浵
		//latitude_x=0.0;
		//longitude_y=0.0;

		 cep_error=0.0;						// ���� : m		
		 eep_major_axis=0.0;			// ���� : m		
		 eep_minor_axis=0.0;			// ���� : m		
		 eep_theta=0.0;						// ���� : degree
		 reallatitude=0.0;					
		 reallongitude=0.0;
	}

}  ;

/*!
 * @typedef   SELSensorPosition
 * @brief			�װ��� ��ġ, ���浵 ��ǥ�� �װ��⿡�� ������ Ÿ���� LOB �� �׸��� ��Ÿ ������ ������ ����ü ����
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @date      2013-09-09 ���� 5:07
 */
struct SELSensorPosition {
	double *pLatitude;		// ���� : degree
	double *pLongitude;		// ���� : degree
	double *pAltitude;		// ���� : meter

	int iCollectorID;

	double *pX;
	double *pY;
    double *pH;

	double *pDOA;				// ���� : degree
	time_t *pTime;				// ���� : ��

	bool *pValid;					// LOB ��ȿ �÷���

	UINT n;

}  ;

struct _COMPENSATION_LL_ {
    double dLongitude;			// �浵 [��]
    double dLatitude;			// ���� [��]

    double dM;

    int iQuardant;
} ;


#ifndef _SELPE_RESULT
#define _SELPE_RESULT
/*!
 * @typedef   SELUTM
 * @brief			UTM�� ����ü ����
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @date      2013-09-09 ���� 5:44
 */
typedef struct {
	double dLongitude;			// �浵 [��], ��ġ ������ �浵 ��.
	double dLatitude;			// ���� [��], ��ġ ������ ���� ��.
	double dAltitude;			// �� [m]

	double dBLongitude;			// �浵 [��], �˰��򿡼� �м��� 1�� �浵 ��
	double dBLatitude;			// ���� [��], �˰��򿡼� �м��� 1�� ���� ��

	double dEasting;
	double dNorthing;

	//unsigned int time_sec;

	double dEEP_major_axis;     // [m]
	double dEEP_minor_axis;     // [m]
	double dEEP_theta;          // [��]
	double dCEP_error;          // [m]

	bool bResult;				// true ����, false ����

} SELPE_RESULT ;
#endif

/*!
 * @typedef   VECTOR
 * @brief			2�����迡���� ���� ���� ǥ��
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @date      2013-10-15 ���� 5:03
 */
struct VECTOR {
	double x;
	double y;

} ;

/*!
 * @typedef   SELINIT
 * @brief			Nonlinear ��ġ ������ ����ϱ� ���� �ʱ� �� ����ü ����
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @date      2013-10-15 ���� 5:03
 */
struct SELINIT {
	VECTOR pe;
	VECTOR speed;
	VECTOR ramp;

}  ;

//! ��ġ ���� �˰��� ����
enum EL_POSIOTN_ESTIMATION_ALGORITHM_OPTION { QUADRATIC=0, DISTANCE_LEAST_SQUARE, LINEAR_LS, NONLINEAR_LS, ADD_PE, AUTO } ;

//! ��ȣ�� ����/���/��� �̵��� ����
enum EL_TARGET_STATE_OPTION { STOP=0, CONSTANT, RAMP } ;

