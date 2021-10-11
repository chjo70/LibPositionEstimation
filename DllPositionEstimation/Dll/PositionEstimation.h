#pragma once

#ifdef MATHFUNCSDLL_EXPORTS
#define MATHFUNCSDLL_API __declspec(dllexport)

#else
#define MATHFUNCSDLL_API __declspec(dllimport)
#endif


#ifdef __cplusplus
extern "C" {
#endif

#ifndef _STR_LOBS
#define _STR_LOBS
// 위치 산출 하기 위한 구조체 정의
typedef struct {
	float fDoa;
	float fLatitude;
	float fLongitude;
	float fAltitude;

	int iCollectorID;

	unsigned int uiLOBID;				// 클러스터링 하면서 구분을 하기 위한 변수 추가

} STR_LOBS;
#endif

#ifndef _SELPE_RESULT
#define _SELPE_RESULT
/*!
	* @typedef	SELUTM
	* @brief	UTM계 구조체 정의
	* @author  조철희 (churlhee.jo@lignex1.com)
	* @date    2013-09-09 오후 5:44
	*/
typedef struct {
	double dLongitude;			// 경도 [도], 위치 보정한 경도 값.
	double dLatitude;			// 위도 [도], 위치 보정한 위도 값.
	double dAltitude;			// 고도 [m]

	double dBLongitude;			// 경도 [도], 알고리즘에서 분석한 1차 경도 값
	double dBLatitude;			// 위도 [도], 알고리즘에서 분석한 1차 위도 값

	double dEasting;
	double dNorthing;

	//unsigned int time_sec;

	double dEEP_major_axis;     // [m] 반장축
	double dEEP_minor_axis;     // [m] 반단축
	double dEEP_theta;          // [도] 기울기, 반장축을 수평각으로 해서 기울기
	double dCEP_error;          // [m], 운형 반지름

	bool bResult;

} SELPE_RESULT ;
#endif

namespace PE
{

	// This class is exported from the MathFuncsDll.dll
	class CPositionEstimation
	{
	private:
		

	public: 
		static MATHFUNCSDLL_API void RunPositionEstimation( SELPE_RESULT *pstSELPE_RESULT, int nLob, STR_LOBS *pstrLOB );
		
	};
}

#ifdef __cplusplus
}
#endif 
