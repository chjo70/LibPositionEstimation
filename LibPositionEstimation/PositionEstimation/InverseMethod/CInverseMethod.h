/**
 * @brief �� ���� -> �Ÿ�, ������ ����
 */

#pragma once

#define _USE_MATH_DEFINES

#include <math.h>

#include "VincentyParam.h" 	 


using namespace std;

namespace GEO {
	const double PIOVER2 = M_PI/2.0;
	const double TWOPI = 6.28318530718;
	const double DE2RA = 0.01745329252;
	const double RA2DE = 57.2957795129;
	const double ERAD = 6378.137;
	const double ERADM = 6378137.0;
	const double AVG_ERAD = 6371.0;
	const double FLATTENING = 1.000000/298.257223563;// Earth flattening (WGS84)
	const double EPS = 0.000000000005;
	const double KM2MI = 0.621371;
	const double GEOSTATIONARY_ALT = 35786.0; // km - approximate value
}

/**
* [�ĺ��� : CLS-GMU-EL-L-PEA]
*
* [������ : SRS-G-SAFR-012]
*
* @class	CInverseMethod
* @brief	��ǥ ���� Ŭ����
*
* (1) Ŭ���� ����
* - �� Ŭ������ �� ��ǥ���� �Ÿ�, ������ ���� ������ִ� Ŭ�����̴�.
*
* (2) �����������
* - ��ġ ���� ���̺귯���� Geolocation_CLobsDll.dll �̴�.
*
* (3) ���� �� ����ó��
* - �ش���� ����
*/
class CInverseMethod
{
	
public:
	CInverseMethod( );
	~CInverseMethod( );

private:
	static CInverseMethod *m_pInstance;				///< ��ü ������
	double m_dDistance;												///< �� �������� �Ÿ� [km]
	double m_dFwdAz;													///< ������
	double m_dRevAz;													///< ������
	
public:
	static CInverseMethod* GetInstance();
	void Finalize();

	bool VincentyInverse( sEllipsoid *e, double lat1, double lon1, double lat2, double lon2 );
	int VincentyInverse( SELDISTLOB *pResult, double lat1, double lon1, double lat2, double lon2 );

	//�� ��ǥ���� �Ÿ�
	double GetDistance( ){ return m_dDistance; }

	//�� ��ǥ���� ���� �ٶ󺸴� ������
	double GetFwdAz( ){ return m_dFwdAz * R2D; }
	double GetRevAz( ){ return m_dRevAz * R2D; }

	double GCAzimuth(double lat1, double lon1, double lat2, double lon2, bool bInitialBearing=false );
	double GCDistance(double lat1, double lon1, double lat2, double lon2);
	double EllipsoidDistance(double lat1, double lon1, double lat2, double lon2);
};

/*!
 * @def				ST_IMA
 * @brief			�ν��Ͻ� ��ü�� ���´�.
 */
#define ST_IMA	CInverseMethod::GetInstance()