/**
 * @brief �� ���� -> �Ÿ�, ������ ����
 */

#pragma once

#define _USE_MATH_DEFINES

#include <math.h>

#include "VincentyParam.h" 	 


using namespace std;


#define DE2RA	((double) 0.01745329252 )
#define RA2DE	((double) 57.2957795129 )
#define ERAD	((double) 6378.137 )
#define ERADM	((double) 6378137.0 )
#define AVG_ERAD	((double) 6371.0 )
#define FLATTENING	((double) (1.000000/298.257223563))			// Earth flattening (WGS84)
#define	EPS			((double) 0.000000000005 )
#define	KM2MI		((double) 0.621371 )
#define	GEOSTATIONARY_ALT	((double) 35786.0 )
	
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
	//static CInverseMethod *m_pInstance;				///< ��ü ������
	double m_dDistance;												///< �� �������� �Ÿ� [km]
	double m_dFwdAz;													///< ������
	double m_dRevAz;													///< ������
	
public:
	bool VincentyInverse( sEllipsoid *e, double lat1, double lon1, double lat2, double lon2 );
	int VincentyInverse( SELDISTLOB *pResult, double lat1, double lon1, double lat2, double lon2 );

	//�� ��ǥ���� �Ÿ�
	double GetDistance( ){ return m_dDistance; }

	//�� ��ǥ���� ���� �ٶ󺸴� ������
	//�ŷڼ�. �Ľ� ������ ���� �ּ�ó����. ���� �ּ� Ǯ��. ����ö.
	double GetFwdAz() { return m_dFwdAz * (double)(180.0 / M_PI); }
	double GetRevAz() { return m_dRevAz * (double)(180.0 / M_PI); }

	double GCAzimuth(double lat1, double lon1, double lat2, double lon2, bool bInitialBearing=false );
	double GCDistance(double lat1, double lon1, double lat2, double lon2);
	double EllipsoidDistance(double lat1, double lon1, double lat2, double lon2);
};

/*!
 * @def				ST_IMA
 * @brief			�ν��Ͻ� ��ü�� ���´�.
 */
//#define ST_IMA	CInverseMethod::GetInstance()