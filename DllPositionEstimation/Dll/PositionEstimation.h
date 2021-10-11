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
// ��ġ ���� �ϱ� ���� ����ü ����
typedef struct {
	float fDoa;
	float fLatitude;
	float fLongitude;
	float fAltitude;

	int iCollectorID;

	unsigned int uiLOBID;				// Ŭ�����͸� �ϸ鼭 ������ �ϱ� ���� ���� �߰�

} STR_LOBS;
#endif

#ifndef _SELPE_RESULT
#define _SELPE_RESULT
/*!
	* @typedef	SELUTM
	* @brief	UTM�� ����ü ����
	* @author  ��ö�� (churlhee.jo@lignex1.com)
	* @date    2013-09-09 ���� 5:44
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

	double dEEP_major_axis;     // [m] ������
	double dEEP_minor_axis;     // [m] �ݴ���
	double dEEP_theta;          // [��] ����, �������� �������� �ؼ� ����
	double dCEP_error;          // [m], ���� ������

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
