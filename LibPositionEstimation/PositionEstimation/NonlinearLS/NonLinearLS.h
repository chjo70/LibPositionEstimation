//////////////////////////////////////////////////////////////////////////
/*!
 * @file      CNonLinearLS.h
 * @brief     CNonLinearLS �˰����� �̿��� ��ġ ����
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @date      2013-09-09 ���� 4:06 
 * @warning   
 */

#pragma once

#include "../PositionEstimationDefine.h"

#include "../Matrix/Matrix.h"


class CNonLinearLS
{
private:
	CMatrix Dk;
	CMatrix Dv;

public:
	CNonLinearLS(void);
	~CNonLinearLS(void);
	void RunNonLinearLS( SELSensorPosition *pResult, SELINIT *pInit, SELPE_RESULT *pSensor, double *pLob, int nEle, EL_TARGET_STATE_OPTION eTargetState );
	int GetTargetState( EL_TARGET_STATE_OPTION eTargetState );
	void CalCEP( SELPositionEstimationResult *pResult, SELPE_RESULT *pEmitterXY, SELPE_RESULT *pSensorXY, double *pTrueLob, int nEle );

};

