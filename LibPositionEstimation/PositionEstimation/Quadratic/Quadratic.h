//////////////////////////////////////////////////////////////////////////
/*!
 * @file      Quadratic.h
 * @brief     Qudratic �˰����� �̿��� ��ġ ����
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @date      2013-09-09 ���� 4:06 
 * @warning   
 */

#pragma once

#include "../PositionEstimationDefine.h"
#include "../ELEmitterDataType.h"

#include "../CEP/AnalyticCEP.h"


class CQuadratic : public CAnalyticCEP
{
private:


public:
	CQuadratic(void);
	~CQuadratic(void);

	bool Run( SELPE_RESULT *pResult, double *pUTMX, double *pUTMY, double *pLob, int nLob );
	void CalCEP( SELPE_RESULT *pResult, SELABTDATA_EXT *pABTExtData );

};

