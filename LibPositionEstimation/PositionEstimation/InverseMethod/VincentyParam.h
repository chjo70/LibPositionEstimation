/**
 * @brief Vincenty's Formulae를 위한 파라미터들
 */
#pragma once

#include <iostream>

struct sEllipsoid
{
	double dMajorAxis; 	//장축(m)
	double dMinorAxis; 	//단축(m)
	double dFlatness; 	//편평도

	sEllipsoid() : dMajorAxis(0.), dMinorAxis(0.), dFlatness(0.)
	{

	}
};

#ifndef M_PI
#define M_PI ( 3.14159265358979323846 )
#endif

//#define D2R		(double) ( M_PI / 180.0 ) 
//#define R2D		(double) ( ( 180.0 / M_PI ) )

//Ellipsoid 종류에 따른 변환 파라미터
#define	WGS84_MAJOR (double) (6378137.0 )
#define WGS84_FLATNESS	(double) ( (double)1. / (double)298.257223563 )
#define WGS84_MINOR		(double) ( WGS84_MAJOR * ( (double) 1.0 - WGS84_FLATNESS ) )

/*!
 * @typedef   SELDISTLOB
 * @brief			Inverse method 를 이용한 거리와 방위각
 * @author    조철희 (churlhee.jo@lignex1.com)
 * @date      2013-09-09 오후 9:12
 */
struct SELDISTLOB {
	double dDistance;
	double fwdlob;
	double revlob;

	SELDISTLOB() {
		dDistance = 0.;
		fwdlob = 0.;
		revlob = 0.;
	}

} ;
