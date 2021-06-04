//////////////////////////////////////////////////////////////////////////
/*!
 * @file      DistanceLeastSquare.h
 * @brief     DistanceLeastSquare �˰����� �̿��� ��ġ ����
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @date      2013-09-09 ���� 4:06 
 * @warning   
 */

#pragma once

#include "../PositionEstimationDefine.h"
#include "../ELEmitterDataType.h"

#include "../Matrix/Matrix.h"

#include "../CEP/AnalyticCEP.h"

class CDistanceLeastSquare : public CAnalyticCEP
{
private:

public:
	CDistanceLeastSquare(void);
	~CDistanceLeastSquare(void);

	bool Run( SELPE_RESULT *pResult, double *pLatitude, double *pLongitude, double *pLob, unsigned int nLob );
	bool CalCEP( SELPE_RESULT *pResult, SELABTDATA_EXT *pABTExtData );

	void CalCEP( SELPositionEstimationResult *pResult, SELPE_RESULT *pEmitterXY, SELPE_RESULT *pSensorXY, double *pTrueLob, int nEle );

};

