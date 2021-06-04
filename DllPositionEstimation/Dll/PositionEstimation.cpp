#include "stdafx.h"

#define MATHFUNCSDLL_EXPORTS

#include "../../LibPositionEstimation/PositionEstimation/PositionEstimationDefine.h"
#include "PositionEstimation.h"

#include "../../LibPositionEstimation/PositionEstimation/PositionEstimationAlg.h"


CPositionEstimationAlg *gpPositionEstimationAlg;


namespace PE
{

	void CPositionEstimation::RunPositionEstimation( SELPE_RESULT *pstSELPE_RESULT, int nLob, STR_LOBS *pstrLOB )
	{
		if( gpPositionEstimationAlg == NULL ) {
			gpPositionEstimationAlg = new CPositionEstimationAlg;
		}

		if( gpPositionEstimationAlg != NULL ) {
			gpPositionEstimationAlg->RunPositionEstimation( pstSELPE_RESULT, nLob, pstrLOB );

			delete gpPositionEstimationAlg;
			gpPositionEstimationAlg = NULL;

		}

	}

}

