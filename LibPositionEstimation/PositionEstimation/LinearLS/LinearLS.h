//////////////////////////////////////////////////////////////////////////
/*!
 * @file      CLinearLS.h
 * @brief     CLinearLS �˰����� �̿��� ��ġ ����
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @date      2013-09-09 ���� 4:06 
 * @warning   
 */

#pragma once

#include "../PositionEstimationDefine.h"

class CLinearLS
{
private:


public:
	CLinearLS(void);
	~CLinearLS(void);
	void RunLinearLS( SELSensorPosition *pResult, SELPE_RESULT *pSensor, double *pLob, int nEle, EL_TARGET_STATE_OPTION eTargetState=STOP );

};

