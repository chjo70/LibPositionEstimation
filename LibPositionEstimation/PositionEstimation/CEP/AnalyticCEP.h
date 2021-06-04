#pragma once

#include "../PositionEstimationDefine.h"

#ifdef _GD_PROJECT_
#define	_SIGMA_OF_DOA_			(1.0)
#else
#define	_SIGMA_OF_DOA_			(1.0)

#endif

class CAnalyticCEP
{
protected:
    int m_nLob;

    double *m_pUTMX;
    double *m_pUTMY;

public:
    CAnalyticCEP(void);
    ~CAnalyticCEP(void);

    bool CalAnalyticNonlinear( SELPE_RESULT *pResult );
};

