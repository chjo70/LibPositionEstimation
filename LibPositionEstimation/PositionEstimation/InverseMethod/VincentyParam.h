/**
 * @brief Vincenty's Formulae�� ���� �Ķ���͵�
 */
#pragma once

#include <iostream>

struct sEllipsoid
{
	double dMajorAxis; 	//����(m)
	double dMinorAxis; 	//����(m)
	double dFlatness; 	//����

	sEllipsoid() : dMajorAxis(0.), dMinorAxis(0.), dFlatness(0.)
	{

	}
};

#ifndef M_PI
#define M_PI ( 3.14159265358979323846 )
#endif

//#define D2R		(double) ( M_PI / 180.0 ) 
//#define R2D		(double) ( ( 180.0 / M_PI ) )

//Ellipsoid ������ ���� ��ȯ �Ķ����
#define	WGS84_MAJOR (double) (6378137.0 )
#define WGS84_FLATNESS	(double) ( (double)1. / (double)298.257223563 )
#define WGS84_MINOR		(double) ( WGS84_MAJOR * ( (double) 1.0 - WGS84_FLATNESS ) )

/*!
 * @typedef   SELDISTLOB
 * @brief			Inverse method �� �̿��� �Ÿ��� ������
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 * @date      2013-09-09 ���� 9:12
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
