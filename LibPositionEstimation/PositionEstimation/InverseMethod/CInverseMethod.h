/**
 * @brief 두 지점 -> 거리, 방위각 산출
 */

#pragma once

#define _USE_MATH_DEFINES

#include <math.h>

#include "VincentyParam.h" 	 


using namespace std;


#define DE2RA	((double) 0.01745329252 )
#define RA2DE	((double) 57.2957795129 )
#define ERAD	((double) 6378.137 )
#define ERADM	((double) 6378137.0 )
#define AVG_ERAD	((double) 6371.0 )
#define FLATTENING	((double) (1.000000/298.257223563))			// Earth flattening (WGS84)
#define	EPS			((double) 0.000000000005 )
#define	KM2MI		((double) 0.621371 )
#define	GEOSTATIONARY_ALT	((double) 35786.0 )
	
/**
* [식별자 : CLS-GMU-EL-L-PEA]
*
* [추적성 : SRS-G-SAFR-012]
*
* @class	CInverseMethod
* @brief	지표 관련 클래스
*
* (1) 클래스 설명
* - 본 클래스는 두 좌표간의 거리, 방위각 등을 계산해주는 클래스이다.
*
* (2) 설계결정사항
* - 위치 산출 라이브러리는 Geolocation_CLobsDll.dll 이다.
*
* (3) 제한 및 예외처리
* - 해당사항 없음
*/
class CInverseMethod
{
	
public:
	CInverseMethod( );
	~CInverseMethod( );

private:
	//static CInverseMethod *m_pInstance;				///< 객체 포인터
	double m_dDistance;												///< 두 지점간의 거리 [km]
	double m_dFwdAz;													///< 방위각
	double m_dRevAz;													///< 방위각
	
public:
	bool VincentyInverse( sEllipsoid *e, double lat1, double lon1, double lat2, double lon2 );
	int VincentyInverse( SELDISTLOB *pResult, double lat1, double lon1, double lat2, double lon2 );

	//두 좌표간의 거리
	double GetDistance( ){ return m_dDistance; }

	//두 좌표에서 서로 바라보는 방위각
	//신뢰성. 파싱 에러가 떠서 주석처리함. 추후 주석 풀것. 윤현철.
	double GetFwdAz() { return m_dFwdAz * (double)(180.0 / M_PI); }
	double GetRevAz() { return m_dRevAz * (double)(180.0 / M_PI); }

	double GCAzimuth(double lat1, double lon1, double lat2, double lon2, bool bInitialBearing=false );
	double GCDistance(double lat1, double lon1, double lat2, double lon2);
	double EllipsoidDistance(double lat1, double lon1, double lat2, double lon2);
};

/*!
 * @def				ST_IMA
 * @brief			인스턴스 객체를 얻어온다.
 */
//#define ST_IMA	CInverseMethod::GetInstance()