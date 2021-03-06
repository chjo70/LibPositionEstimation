// Util.h: interface for the CUtil class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_UTIL_H__43C28AC5_E619_4079_AEF1_4C680C7CA6B5__INCLUDED_)
#define AFX_UTIL_H__43C28AC5_E619_4079_AEF1_4C680C7CA6B5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#define	KM_PER_DEGREE_FOR_LATITUDE						(111.111)				// [km/deg]
#define	KM_PER_DEGREE_FOR_LONGITUDE						(88.884)				// [km/deg]

BOOL ELCompSwitchLevel( int *pSeries1, int *pSeries2, int coSeries, int margin );
BOOL CompTOAMeanDiff( unsigned long long int x, unsigned long long int y, int thresh );
BOOL CompMarginDiffFFF( float x, float fy1, float fy2, float thresh );
BOOL CompAoaDiff( float x, float y, float fthresh );
float AoaDiff( float x, float y);
void LogPrint( char *format, ... );

bool Is_FZero( const float &value );
bool Is_DZero( const double &dvalue );
bool Is_FNotZero( const float &value );
bool Is_DNotZero( const double &dvalue );

int CalOverlapSpace( int hgh1, int low1, int hgh2, int low2 );
float CalOverlapSpaceF( float hgh1, float low1, float hgh2, float low2 );

template <typename T>
T _diffabs( T x, T y)
{

    if (x > y) {
        return x - y;
    }
    else {
        return y - x;
    }

}

template <typename T>
BOOL CompMeanDiff( T x, T y, T thresh )
{
	T diff;
	BOOL bRet;

	diff = _diffabs<T>( x, y );

	if( diff <= thresh ) {
		bRet = TRUE;
	}
	else {
		bRet = FALSE;
	}

	return bRet;
}

template <typename T>
BOOL ELCompSwitchLevel( T *pSeries1, T *pSeries2, int coSeries, T margin )
{
	int j, k; 
	T iDiff;
	int index1;

	T *pSeries;

	BOOL bRet=FALSE;

	if( coSeries == 0 ) {
		// bRet = FALSE;
	}
	else {
		for( j=0 ; j < coSeries ; ++j ) {
			pSeries = pSeries2;
			bRet = TRUE;
			iDiff = 0;
			for( k=j ; k < coSeries+j ; ++k ) {
				index1 = k % coSeries;
				bRet = CompMeanDiff<T>( pSeries1[index1], *pSeries, margin );
				if( FALSE == bRet ) {
					break;
				}
				iDiff += _diffabs<T>( pSeries1[index1], *pSeries );
				++ pSeries;
			}

			if( TRUE == bRet ) {
				// bRet = TRUE;
				break;
			}

		}

	}

	return bRet;

}

template <typename T>
BOOL CompSwitch2Level( T *pSeries1, T *pSeries2, int coSeries1, int coSeries2, T margin )
{
	int i, j;
	int index1;

	T *pSeries;

	BOOL bRet=FALSE;

	if( coSeries1 == 0 || coSeries2 == 0 ) {
		// bRet = FALSE;
	}

	else {
		if( coSeries1 == coSeries2 ) {
			bRet = ELCompSwitchLevel<T>( pSeries1, pSeries2, coSeries1, margin );
		}
		else if( coSeries1 < coSeries2 ) {
			for( i=0 ; i < coSeries2-coSeries1 ; ++i ) {
				pSeries = pSeries2;
				index1 = i;
				for( j=0 ; j < coSeries1 ; ++j ) {
					bRet = CompMeanDiff( pSeries1[index1], *pSeries, margin );
					if( bRet == FALSE ) {
						break;
					}
					++ pSeries;
					++ index1;
				}

				if( bRet == TRUE ) {
					break;
				}
			}
		}
		else {
			for( i=0 ; i < coSeries1-coSeries2 ; ++i ) {
				pSeries = pSeries2;
				index1 = i;
				for( j=0 ; j < coSeries2 ; ++j ) {
					bRet = CompMeanDiff( pSeries1[index1], *pSeries, margin );
					if( bRet == FALSE ) {
						break;
					}
					++ pSeries;
					++ index1;
				}

				if( bRet == TRUE ) {
					break;
				}
			}
		}

	}

	return bRet;

}

//////////////////////////////////////////////////////////////////////////
/*! \brief    ???? ?? ?????? ???????? ?????? ???? ???? ???????? TRUE, ?????? ?????? FALSE?? ????????.
		\author   ??????
		\param    x ???????? int
		\param    y1 ???????? int
		\param    y2 ???????? int
		\param    thresh ???????? int
		\return   BOOL
		\version  0.0.1
		\date     2008-02-19 17:44:30
		\warning
*/
template <typename T>
BOOL CompMarginDiff( T x, T iy1, T iy2, T fthresh )
{
	BOOL bRet;

	if( ( x >= iy1-fthresh ) && ( x <= iy2+fthresh ) ) {
		bRet = TRUE;
	}
	else {
		bRet = FALSE;
	}
	return bRet;
}

template <typename T>
T CalOverlapSpace(T hgh1, T low1, T hgh2, T low2)
{
	if( hgh1 < low2 || hgh2 < low1 )				// |-------|		|--------|
		return 0;															//			 |---|

	if( hgh1 == low2 || hgh2 == low1 )			// debug, 99-12-22 09:36:19
		return 1;

	if( low1 >= low2 &&	hgh1 <= hgh2 ) 			//          |--------|
		return hgh1 - low1 + 1;								//    |------------------|

	if( low1 <= low2 && hgh1 >= hgh2 )			//    |------------------|
		return hgh2 - low2 + 1;			 					//          |--------|

	if( low1 <= hgh2 && low2 <= low1 )			//          |------------|
		return ( hgh2 - low1 + 1);     			//    |-----------|

	if( hgh1 >= low2 && hgh1 <= hgh2 )   		//    |-----------|
		return ( hgh1 - low2 + 1 );						//          |------------|

	return 0;
}


// #ifdef __cplusplus
// }
// #endif

#endif // !defined(AFX_UTIL_H__43C28AC5_E619_4079_AEF1_4C680C7CA6B5__INCLUDED_)
