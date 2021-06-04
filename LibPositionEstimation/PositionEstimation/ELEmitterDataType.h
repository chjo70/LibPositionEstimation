#ifndef _H_EL_EMITTER_DATATYPE
#define _H_EL_EMITTER_DATATYPE

//#include "pch.h"
#include <CString>
#include <string>
#include <list>
using namespace std;

#include "ELMsgDefn.h"

//#include "../Identify/ELCEDLibDataType2.h"


/**
* [�ĺ��� : D-GR-SDD-XXX]
* [������ : R-GR-SRS-XXX]
*
* @file ELEmitterDataType.h
* @brief ���ü ���â���� View/Logic �� �߻��ϴ� �������̽� data ����
* @author ������
* @date 2013.5.25
*   
* [�����̷�]
*  2013.5.25. �ű��ۼ�
*  2013.6.18. std::string strIsFisintTask ������ �߰�
*/

// ��ġ �����ϱ� ���� �ִ� LOB ���� ����
#define MAX_OF_LOBS_PE								    (30)					// LOB �ִ� ����


/************************************************************************************
*   ELINT Logic -> View�� ���Ǵ� �ڷ��� ����ü
*************************************************************************************/

/* 1.  E_EL_LV_SHOW_EMITTER_DATA �۽Ž� data_in�� ���Ǵ� ����ü. ���Ͼ���*/
#ifndef ST_EMIITERINFO_DISP
#define ST_EMIITERINFO_DISP
// stEmitterInfoDisp �� �ߺ� ���Ǹ� �������� #ifndef.. ���

#define SIZE_OF_CHAR_EMITTER_INFO	15
#define SIZE_OF_CHAR_MID_INFO		32

//#define MAX_IDCANDIDATE					10

// ������.
typedef struct stEmitterInfoDisp
{
	string strEmitterNum;
	string strDirection;
	string strTaskId;
	string strIsFisintTask; // 2013.6.18. �߰�
	string strFrqTypeRange;
	string strPRITypeRange;
	string strSJ;
	string strPwRange;
	string strPaRange;
	string strScanTypePeriod;
	string strElintNoteAir;
	string strPrimFuncAir;
	string strPlatformAir;
	string strIdNumberAir;
	string strConsistRateAir;
	string strElintNoteGnd;
	string strPrimFuncGnd;
	string strPlatformGnd;
	string strIdNumberGnd;
	string strConsisitRateGnd;
	string strEstimated;
	string strRptStat;
	string strSOI;
	string strPosEstimated;
	string strActivatedOrNot;
}_EMITTER_INFO;

// MAX_COLUMN ���� ���� 70�� �̸� �̻��� �ɶ� �� ���� ���� �������Ѿ� �Ѵ�. , ��ö��
// 68���� �ϸ� S/W�� �Ҿ���.
// UpdateHeaderCol()�� ���� ��ó�� ��. �� �̷��� ��������� �𸣰ڳ�...
#define NUM_OF_EMITTER_FIELD_NAME 90
static const char* strEmitterFieldName[] =
{
	"����",
	//"SDF ID",
	"AET��ȣ",
	"��ȣ�Ϸù�ȣ",
	"���� ����",

	"[G]PIN NR",
	"[G]��ɺ�ȣ",
	"[G]����",
	"[A]E-NOT",
	"[G]E-NOT",
	"[G]�ĺ�����",
	"LOB#",
	"����/�ð�",
	"�ĺ��ڵ�",

	"��ȣ����",
	"����[��]",
	"�ּ� ����[��]",
	"�ִ� ����[��]",
	"FISINT����",
	"���ļ�����",
	"���ļ�[MHz]",
	"�ּ� ���ļ�[MHz]",
	"�ִ� ���ļ�[MHz]",
	"���ļ�����/�ֱ�[ms]",
	"PRI����",
	"PRI[us]",
	"PRI �ּ�[us]",
	"PRI �ִ�[us]",

	"PRI����/�ֱ�[ms]",
	"S/J",

	"PW[us]",
	"PW �ּ�[us]",
	"PW �ִ�[us]",
	"PA[dBm]",
	"PA �ּ�[dBm]",
	"PA �ִ�[dBm]",
	"��ĵ����/�ֱ�[ms]",
	"�ؼ�/�ŷڵ�",
	"�켱����",
	"��������",
	"�޽���/�׷�",

	"����ID",
	"������",
	"Ž���뿪��ȣ",
	"��������[��]",
	"���ļ�����[MHz]",
	"PRF [PPS]",
	"�ּ� PRF[PPS]",
	"�ִ� PRF[PPS]",

	//"[A]PIN NR",
	"[A]��ɺ�ȣ",
	"[A]����",
	"[A]����",
	"[A]����ȣ",
	"[A]��ġ��[%]",
	"[A]��ġ��[%]",
	"[A]������ǥ",
	"[A]��[km]",
	"[A]�ŷڱ���[m]",
	"[A]CEP[km]",
	"[A]����[km]",
	"[A]����[km]",
	"[A]�������[��]",

	"[G]����",
	"[G]����ȣ",
	"[G]��ġ��[%]",
	"[G]��ġ��[%]",
	"[G]������ǥ(����)",
	"[G]������ǥ(�浵)",
	"[G]�ŷڱ���[m]",
	"[G]CEP[km]",
	"[G]����[km]",
	"[G]����[km]",
	"[G]�������[��]",
	"[G]������",
	"[G]�ĺ��ĺ�",
	"[G]�Ÿ� ����[km]",
	"[G]�Ÿ� ����[nm]",
	"[G]EOB ID",
	"��������",
	"�����溸",
	"���⿩��",
	"���ɽ�ȣ",
	"PDW ���忩��",
	"��������[PDW/IQ]",
	"��ȣ��",
	"PDW ��Ʈ",
	"Ȱ������",

	"��ȣID",

	"�ӹ���",
	"����ȣ",
	"LINK��ȣ", // 2015.1.5. SDF ID ��ſ�  LINK ��ȣ �ʵ�� �ٲ�. �ʵ� �̸��� �׳� �ξ��� !!
	
};
//
typedef struct stEmitterInfo
{
	CString strLinkNum;	 //#FA_Q_2502_T1						// "SDF ID";
	CString strEmitterNum;					// "���ü��ȣ",
	CString strSignalNum;						// "��ȣ�Ϸù�ȣ",
	CString strMissionId;						// "�ӹ���",
	CString strTaskId;							// "����ID",
	CString strTaskName;						// "������",
	CString strSrchBandId;					// "Ž���뿪 id"
	CString strDirection;						// "����[��]",
	CString strMinDirection;				// "�ּ� ����[��]",
	CString strMaxDirection;				// "�ִ� ����[��]",
	CString strDirectionDev;				// "��������[��]",
	CString strFisintTask;					// "FISINT����",
	CString strFrqType;							// "���ļ�����",
	CString strFrq;									// "���ļ�[Hz]",
	CString strMinFrq;							// "�ּ� ���ļ�[Hz]",
	CString strMaxFrq;							// "�ִ� ���ļ�[Hz]",
	CString strFrqPatternRange;			// "���ļ���������/����[ms]",
	CString strFrqRangeDev;					// "���ļ�����[Hz]",
	CString strPriority;						// "�켱����",
	CString strTime;								// "����/�ð�",
	CString strModType;							// "��������",
	CString strSigType;		//#FA_Q_2502_T1					// "��ȣ����",
	CString strPw;									// "PW[ns]",
	CString strMinPw;								// "�ּ� PW[ns]",
	CString strMaxPw;								// "�ִ� PW[ns]",
	CString strPa;									// "PA[dB]",
	CString strMinPa;								// "�ּ� PA[dBm]",
	CString strMaxPa;								// "�ִ� PA[dBm]",
	CString strPriType;							// "PRI����",
	CString strPri;									// "PRI[ns]",
	CString strMinPri;							// "�ּ� PRI",
	CString strMaxPri;							// "�ִ� PRI",
	CString strPRIPatternRange;			// "PRI��������/����[ns]",
	CString strPrf;									// "PRF����[PPS]",
	CString strMinPrf;							// "PRF����[PPS]",
	CString strMaxPrf;							// "PRF����[PPS]",
	CString strPulsePerGrp;					// "�޽���/�׷�",
	CString strSJ;									// "S/J",
	CString strScanTypePeroid;			// "��ĵ����/�ֱ�",
	CString strAirElintNotation;		// "[A]E-NOT",
	CString strAirFunc;							// "[A]��ɺ�ȣ",
	CString strAirSiteFunc;					// "[A]����Ʈ���",
	CString strAirThreat;						// "[A]����",
	CString strAirEqupNo;						// "[A]����ȣ",
	CString strAirConsistancyRatio;	// "[A]��ġ��[%]",
	CString strAirConsistancy;			// "[A]��ġ��",
	CString strAirGeoCoord;					// "[A]������ǥ",
	CString strAirAlt;							// "[A]��[m]",
	CString strAirConfIterval;			// "[A]�ŷڱ���[m]",
	CString strAirCep;							// "[A]CEP[m]",
	CString strAirMajorAxis;				// "[A]����[m]",
	CString strAirMinorAxis;				// "[A]����(m)",
	CString strAirMajorAxisAzimuth; // "[A]�������[��]",
	CString strGndElintNotation;		// "[G]E-NOT",
	CString strGndPinNr;						// "[G]PIN NR",
	CString strGndFunc;							// "[G]��ɺ�ȣ",
	CString strGndSiteFunc;					// "[G]����Ʈ���",
	CString strGndThreat;						// "[G]����",
	CString strGndEqupNo;						// "[G]����ȣ",
	CString strGndConsistancyRatio;	// "[G]��ġ��[%]",
	CString strGndConsistancy;			// "[G]��ġ��",
	CString strGndGeoLongitude;			// "[G]������ǥ(�浵)",
	CString strGndGeoLatitude;			// "[G]������ǥ(����)",
	CString strGndAlt;							// "[G]��[m]",
	CString strGndConfIterval;			// "[G]�ŷڱ���[m]",
	CString strGndCep;							// "[G]CEP[m]",
	CString strGndMajorAxis;				// "[G]����[m]",
	CString strGndMinorAxis;				// "[G]����(m)",
	CString strGndMajorAxisAzimuth; // "[G]�������[��]",												
	CString strGndBaseName;					// "[G]������",
	CString strGndIdInfo;						// "[G]�ĺ�����",
	CString strGndCandidate;				// "[G]�ĺ� �ĺ� ���",
	CString strGndDistError;				// "[G]�Ÿ� ����[km]",
	CString strGndDistErrorToNM;		// "[G]�Ÿ� ����[nm]",
	CString strGndEOBId;						// "[G] EOB ID",
	CString strFinalReport;					// "��������",							
	CString strFinalAlarm;					// "�����溸",							",
	CString strPosEstimation;				// "���⿩��",							
	CString strSOI;									// "���ɽ�ȣ"							
	CString strIsStorePDW;					// "PDW ���忩��"
	CString strNumOfCol;						// "��������[PDW/IQ]"
	CString strNumOfAmbigiousBeam;	// "��ȣ��",
	CString strNumOfPDWSet;					// "PDW ��Ʈ"
	CString strActivatedOrNot;			// "Ȱ������"
	CString strLobId;								// "LOB #" . 2014.10.8. �߰�. ������
	CString strDtctId;							// "Dtct ID" . 2015.3.3. �߰�. ������
	CString strIdCode;							// "�ĺ� �ڵ�" . 2015.3.31. �߰�. ��ö��
	CString strPolization;					// "�ؼ�" . 2015.4.1. �߰�. ��ö��
	CString strTaskOperator;				// "���� ����", 2015.4.7. �߰�. ��ö��
	CString strGndModeSymbol;					// ��� ��ȣ

	stEmitterInfo()
	{
		strLinkNum="";						
		strEmitterNum="";					
		strSignalNum="";					
		strMissionId="";						
		strTaskId="";							
		strTaskName="";					
		strSrchBandId="";					
		strDirection="";						
		strMinDirection="";				// 
		strMaxDirection="";				
		strDirectionDev="";				
		strFisintTask="";					// 
		strFrqType="";						
		strFrq="";								
		strMinFrq="";							
		strMaxFrq="";						
		strFrqPatternRange="";			
		strFrqRangeDev="";				
		strPriority="";						// 
		strTime="";							
		strModType="";						
		strSigType="";						
		strPw="";								
		strMinPw="";							
		strMaxPw="";							
		strPa="";								
		strMinPa="";							
		strMaxPa="";							
		strPriType="";						
		strPri="";								
		strMinPri="";							
		strMaxPri="";							
		strPRIPatternRange="";			
		strPrf="";								
		strMinPrf="";							
		strMaxPrf="";							
		strPulsePerGrp="";					
		strSJ="";								
		strScanTypePeroid="";			
		strAirElintNotation="";		// "[A
		strAirFunc="";						
		strAirSiteFunc="";					
		strAirThreat="";						
		strAirEqupNo="";					
		strAirConsistancyRatio="";	// 
		strAirConsistancy="";			// 
		strAirGeoCoord="";				
		strAirAlt="";							
		strAirConfIterval="";			// 
		strAirCep="";							
		strAirMajorAxis="";				// 
		strAirMinorAxis="";				// 
		strAirMajorAxisAzimuth=""; //
		strGndElintNotation="";		// 
		strGndPinNr="";						
		strGndFunc="";						
		strGndSiteFunc="";					
		strGndThreat="";					
		strGndEqupNo="";					
		strGndConsistancyRatio="";	
		strGndConsistancy="";			
		strGndGeoLongitude="";			
		strGndGeoLatitude="";			
		strGndAlt="";							
		strGndConfIterval="";			// 
		strGndCep="";						
		strGndMajorAxis="";				
		strGndMinorAxis="";				
		strGndMajorAxisAzimuth=""; 
			strGndBaseName="";				
		strGndIdInfo="";						
		strGndCandidate="";				
		strGndDistError="";				
		strGndDistErrorToNM="";		
		strGndEOBId="";					
		strFinalReport="";					
		strFinalAlarm="";					
		strPosEstimation="";				
		strSOI="";								
		strIsStorePDW="";					
		strNumOfCol="";					
		strNumOfAmbigiousBeam="";	
		strNumOfPDWSet="";				
		strActivatedOrNot="";			
		strLobId="";							
		strDtctId="";							
		strIdCode="";							
		strPolization="";					// 
		strTaskOperator="";				
		strGndModeSymbol="";			
	};

//	void stEmitterInfo::SetInterCommInfo(BYTE *i_pbyInData)
//	{
//		int nCount=0;
//		USHORT nSize=0;
//		
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strLinkNum.SetString((const char*)&i_pbyInData[nCount], nSize) ;								// "SDF ID";
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strEmitterNum.SetString((const char*)&i_pbyInData[nCount], nSize);					// "���ü��ȣ",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strSignalNum.SetString((const char*)&i_pbyInData[nCount], nSize);						// "��ȣ�Ϸù�ȣ",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strMissionId.SetString((const char*)&i_pbyInData[nCount], nSize);						// "�ӹ���",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strTaskId.SetString((const char*)&i_pbyInData[nCount], nSize);							// "����ID",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strSrchBandId.SetString((const char*)&i_pbyInData[nCount], nSize);					// "Ž���뿪 id"
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strDirection.SetString((const char*)&i_pbyInData[nCount], nSize);						// "����[��]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strDirectionDev.SetString((const char*)&i_pbyInData[nCount], nSize);				// "��������[��]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strFisintTask.SetString((const char*)&i_pbyInData[nCount], nSize);					// "FISINT����",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strFrqType.SetString((const char*)&i_pbyInData[nCount], nSize);				// "���ļ�����",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strFrq.SetString((const char*)&i_pbyInData[nCount], nSize);				// "���ļ�[Hz]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strMinFrq.SetString((const char*)&i_pbyInData[nCount], nSize);				// "�ּ� ���ļ�[Hz]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strMaxFrq.SetString((const char*)&i_pbyInData[nCount], nSize);				// "�ִ� ���ļ�[Hz]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strFrqPatternRange.SetString((const char*)&i_pbyInData[nCount], nSize);			// "���ļ���������/����[ms]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strFrqRangeDev.SetString((const char*)&i_pbyInData[nCount], nSize);					// "���ļ�����[Hz]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strPriority.SetString((const char*)&i_pbyInData[nCount], nSize);						// "�켱����",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strTime.SetString((const char*)&i_pbyInData[nCount], nSize);								// "����/�ð�",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strModType.SetString((const char*)&i_pbyInData[nCount], nSize);							// "��������",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strSigType.SetString((const char*)&i_pbyInData[nCount], nSize);							// "��ȣ����",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strPw.SetString((const char*)&i_pbyInData[nCount], nSize);							// "PW����[ns]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strMinPw.SetString((const char*)&i_pbyInData[nCount], nSize);							// "PW����[ns]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strMaxPw.SetString((const char*)&i_pbyInData[nCount], nSize);							// "PW����[ns]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strPa.SetString((const char*)&i_pbyInData[nCount], nSize);							// "PA����[dB]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strMinPa.SetString((const char*)&i_pbyInData[nCount], nSize);							// "PA����[dB]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strMaxPa.SetString((const char*)&i_pbyInData[nCount], nSize);							// "PA����[dB]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strPriType.SetString((const char*)&i_pbyInData[nCount], nSize);							// "PRI����",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strPri.SetString((const char*)&i_pbyInData[nCount], nSize);									// "PRI[ns]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strMinPri.SetString((const char*)&i_pbyInData[nCount], nSize);							// "�ּ� PRI[ns]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strMaxPri.SetString((const char*)&i_pbyInData[nCount], nSize);							// "�ִ� PRI[ns]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strPRIPatternRange.SetString((const char*)&i_pbyInData[nCount], nSize);			// "PRI��������/����[ns]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strPrf.SetString((const char*)&i_pbyInData[nCount], nSize);						// "PRF����[PPS]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strMinPrf.SetString((const char*)&i_pbyInData[nCount], nSize);						// "PRF����[PPS]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strMaxPrf.SetString((const char*)&i_pbyInData[nCount], nSize);						// "PRF����[PPS]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strPulsePerGrp.SetString((const char*)&i_pbyInData[nCount], nSize);					// "�޽���/�׷�",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strSJ.SetString((const char*)&i_pbyInData[nCount], nSize);									// "S/J",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strScanTypePeroid.SetString((const char*)&i_pbyInData[nCount], nSize);			// "��ĵ����/�ֱ�",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strAirElintNotation.SetString((const char*)&i_pbyInData[nCount], nSize);		// "[A]E-NOT",
//		nCount+=nSize;
//
//		//memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		//nCount+=sizeof(USHORT);
//		//strAirPinNr.SetString((const char*)&i_pbyInData[nCount], nSize);						// "[A]PIN NR",
//		//nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strAirFunc.SetString((const char*)&i_pbyInData[nCount], nSize);							// "[A]��ɺ�ȣ",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strAirSiteFunc.SetString((const char*)&i_pbyInData[nCount], nSize);					// "[A]����Ʈ���",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strAirThreat.SetString((const char*)&i_pbyInData[nCount], nSize);						// "[A]����",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strAirEqupNo.SetString((const char*)&i_pbyInData[nCount], nSize);						// "[A]����ȣ",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strAirConsistancy.SetString((const char*)&i_pbyInData[nCount], nSize);			// "[A]��ġ��[%]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strAirGeoCoord.SetString((const char*)&i_pbyInData[nCount], nSize);					// "[A]������ǥ",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strAirAlt.SetString((const char*)&i_pbyInData[nCount], nSize);							// "[A]��[m]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strAirConfIterval.SetString((const char*)&i_pbyInData[nCount], nSize);			// "[A]�ŷڱ���[m]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strAirCep.SetString((const char*)&i_pbyInData[nCount], nSize);							// "[A]CEP[m]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strAirMajorAxis.SetString((const char*)&i_pbyInData[nCount], nSize);				// "[A]����[m]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strAirMinorAxis.SetString((const char*)&i_pbyInData[nCount], nSize);				//"[A]����(m)",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strAirMajorAxisAzimuth.SetString((const char*)&i_pbyInData[nCount], nSize); //"[A]�������[��]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strGndElintNotation.SetString((const char*)&i_pbyInData[nCount], nSize);		// "[G]E-NOT",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strGndPinNr.SetString((const char*)&i_pbyInData[nCount], nSize);						// "[G]PIN NR",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strGndFunc.SetString((const char*)&i_pbyInData[nCount], nSize);							// "[G]��ɺ�ȣ",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strGndSiteFunc.SetString((const char*)&i_pbyInData[nCount], nSize);					// "[G]����Ʈ���",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strGndThreat.SetString((const char*)&i_pbyInData[nCount], nSize);						// "[G]����",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strGndEqupNo.SetString((const char*)&i_pbyInData[nCount], nSize);						// "[G]����ȣ",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strGndConsistancy.SetString((const char*)&i_pbyInData[nCount], nSize);			// "[G]��ġ��[%]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strGndGeoLongitude.SetString((const char*)&i_pbyInData[nCount], nSize);			//"[G]������ǥ(�浵)",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strGndGeoLatitude.SetString((const char*)&i_pbyInData[nCount], nSize);			//"[G]������ǥ(����)",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strGndAlt.SetString((const char*)&i_pbyInData[nCount], nSize);							// "[G]��[m]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strGndConfIterval.SetString((const char*)&i_pbyInData[nCount], nSize);			// "[G]�ŷڱ���[m]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strGndCep.SetString((const char*)&i_pbyInData[nCount], nSize);							// "[G]CEP[m]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strGndMajorAxis.SetString((const char*)&i_pbyInData[nCount], nSize);				// "[G]����[m]",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strGndMinorAxis.SetString((const char*)&i_pbyInData[nCount], nSize);				//"[G]����(m)",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strGndMajorAxisAzimuth.SetString((const char*)&i_pbyInData[nCount], nSize); //"[G]�������[��]",			
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strFinalReport.SetString((const char*)&i_pbyInData[nCount], nSize);					// "��������",			
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strFinalAlarm.SetString((const char*)&i_pbyInData[nCount], nSize);					// "�����溸",							",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strPosEstimation.SetString((const char*)&i_pbyInData[nCount], nSize);				// "���⿩��",							
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strSOI.SetString((const char*)&i_pbyInData[nCount], nSize);									// "���ɽ�ȣ"	
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strIsStorePDW.SetString((const char*)&i_pbyInData[nCount], nSize);					// "PDW ���忩��"
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strNumOfCol.SetString((const char*)&i_pbyInData[nCount], nSize);						// "��������[PDW/IQ]"
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strNumOfAmbigiousBeam.SetString((const char*)&i_pbyInData[nCount], nSize);	// "��ȣ��",
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strNumOfPDWSet.SetString((const char*)&i_pbyInData[nCount], nSize);					// "PDW ��Ʈ"
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strIdCode.SetString((const char*)&i_pbyInData[nCount], nSize);							// "�ĺ� �ڵ�"
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strPolization.SetString((const char*)&i_pbyInData[nCount], nSize);					// "�ؼ�"
//		nCount+=nSize;
//
//		memcpy(&nSize, &i_pbyInData[nCount], sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		strTaskOperator.SetString((const char*)&i_pbyInData[nCount], nSize);		// "���� ����"
//		nCount+=nSize;
//
//	}
//	int stEmitterInfo::GetInterCommInfo(BYTE * i_pbyOutData)
//	{
//		int nCount=0;
//		USHORT nSize=0;
//		
//		nSize=strLinkNum.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strLinkNum.GetBuffer(), nSize);								// "SDF ID";
//		nCount+=nSize;
//		
//		nSize=strEmitterNum.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strEmitterNum.GetBuffer(), nSize);					// "���ü��ȣ",
//		nCount+=nSize;
//
//		nSize=strSignalNum.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strSignalNum.GetBuffer(), nSize);						// "��ȣ�Ϸù�ȣ",
//		nCount+=nSize;
//
//		nSize=strMissionId.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strMissionId.GetBuffer(), nSize);						// "�ӹ���",
//		nCount+=nSize;
//
//		nSize=strTaskId.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strTaskId.GetBuffer(), nSize);							// "����ID",
//		nCount+=nSize;
//
//		nSize=strSrchBandId.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strSrchBandId.GetBuffer(), nSize);					// "Ž���뿪 id"
//		nCount+=nSize;
//
//		nSize=strDirection.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strDirection.GetBuffer(), nSize);						// "����[��]",
//		nCount+=nSize;
//
//		nSize=strDirectionDev.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strDirectionDev.GetBuffer(), nSize);				// "��������[��]",
//		nCount+=nSize;
//
//		nSize=strFisintTask.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strFisintTask.GetBuffer(), nSize);					// "FISINT����",
//		nCount+=nSize;
//
//		nSize=strFrqType.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strFrqType.GetBuffer(), nSize);				// "���ļ�����",
//		nCount+=nSize;
//
//		nSize=strFrq.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strFrq.GetBuffer(), nSize);				// "���ļ�[Hz]",
//		nCount+=nSize;
//
//		nSize=strMinFrq.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strMinFrq.GetBuffer(), nSize);				// "�ּ� ���ļ�",
//		nCount+=nSize;
//
//		nSize=strMaxFrq.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strMaxFrq.GetBuffer(), nSize);				// "�ִ� ���ļ�",
//		nCount+=nSize;
//
//		nSize=strFrqPatternRange.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strFrqPatternRange.GetBuffer(), nSize);			// "���ļ���������/����[ms]",
//		nCount+=nSize;
//
//		nSize=strFrqRangeDev.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strFrqRangeDev.GetBuffer(), nSize);					// "���ļ�����[Hz]",
//		nCount+=nSize;
//
//		nSize=strPriority.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strPriority.GetBuffer(), nSize);						// "�켱����",
//		nCount+=nSize;
//
//		nSize=strTime.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strTime.GetBuffer(), nSize);								// "����/�ð�",
//		nCount+=nSize;
//
//		nSize=strModType.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strModType.GetBuffer(), nSize);							// "��������",
//		nCount+=nSize;
//
//		nSize=strSigType.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strSigType.GetBuffer(), nSize);							// "��ȣ����",
//		nCount+=nSize;
//
//		nSize=strPw.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strPw.GetBuffer(), nSize);							// "PW����[ns]",
//		nCount+=nSize;
//
//		nSize=strMinPw.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strMinPw.GetBuffer(), nSize);							// "PW����[ns]",
//		nCount+=nSize;
//
//		nSize=strMaxPw.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strMaxPw.GetBuffer(), nSize);							// "PW����[ns]",
//		nCount+=nSize;
//
//		nSize=strPa.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strPa.GetBuffer(), nSize);							// "PA����[dB]",
//		nCount+=nSize;
//
//		nSize=strMinPa.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strMinPa.GetBuffer(), nSize);							// "PA����[dB]",
//		nCount+=nSize;
//
//		nSize=strMaxPa.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strMaxPa.GetBuffer(), nSize);							// "PA����[dB]",
//		nCount+=nSize;
//
//		nSize=strLinkNum.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strPriType.GetBuffer(), nSize);							// "PRI����",
//		nCount+=nSize;
//
//		nSize=strLinkNum.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strPri.GetBuffer(), nSize);									// "PRI[ns]",
//		nCount+=nSize;
//
//		nSize=strLinkNum.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strMinPri.GetBuffer(), nSize);							// "�ּ� PRI[ns]",
//		nCount+=nSize;
//
//		nSize=strLinkNum.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strMaxPri.GetBuffer(), nSize);							// "�ִ� PRI[ns]",
//		nCount+=nSize;
//
//		nSize=strPRIPatternRange.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strPRIPatternRange.GetBuffer(), nSize);			// "PRI��������/����[ns]",
//		nCount+=nSize;
//
//		nSize=strPrf.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strPrf.GetBuffer(), nSize);						// "PRF����[PPS]",
//		nCount+=nSize;
//
//		nSize=strMinPrf.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strMinPrf.GetBuffer(), nSize);						// "PRF����[PPS]",
//		nCount+=nSize;
//
//		nSize=strMaxPrf.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strMaxPrf.GetBuffer(), nSize);						// "PRF����[PPS]",
//		nCount+=nSize;
//
//		nSize=strPulsePerGrp.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strPulsePerGrp.GetBuffer(), nSize);					// "�޽���/�׷�",
//		nCount+=nSize;
//
//		nSize=strSJ.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strSJ.GetBuffer(), nSize);									// "S/J",
//		nCount+=nSize;
//
//		nSize=strScanTypePeroid.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strScanTypePeroid.GetBuffer(), nSize);			// "��ĵ����/�ֱ�",
//		nCount+=nSize;
//
//		nSize=strAirElintNotation.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strAirElintNotation.GetBuffer(), nSize);		// "[A]E-NOT",
//		nCount+=nSize;
//
//// 		nSize=strAirPinNr.GetLength();
//// 		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//// 		nCount+=sizeof(USHORT);
//// 		memcpy(&i_pbyOutData[nCount], strAirPinNr.GetBuffer(), nSize);						// "[A]PIN NR",
//// 		nCount+=nSize;
//
//		nSize=strAirFunc.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strAirFunc.GetBuffer(), nSize);							// "[A]��ɺ�ȣ",
//		nCount+=nSize;
//
//		nSize=strAirSiteFunc.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strAirSiteFunc.GetBuffer(), nSize);					// "[A]����Ʈ���",
//		nCount+=nSize;
//
//		nSize=strAirThreat.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strAirThreat.GetBuffer(), nSize);						// "[A]����",
//		nCount+=nSize;
//
//		nSize=strAirEqupNo.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strAirEqupNo.GetBuffer(), nSize);						// "[A]����ȣ",
//		nCount+=nSize;
//
//		nSize=strAirConsistancy.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strAirConsistancy.GetBuffer(), nSize);			// "[A]��ġ��[%]",
//		nCount+=nSize;
//
//		nSize=strAirGeoCoord.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strAirGeoCoord.GetBuffer(), nSize);					// "[A]������ǥ",
//		nCount+=nSize;
//
//		nSize=strAirAlt.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strAirAlt.GetBuffer(), nSize);							// "[A]��[m]",
//		nCount+=nSize;
//
//		nSize=strAirConfIterval.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strAirConfIterval.GetBuffer(), nSize);			// "[A]�ŷڱ���[m]",
//		nCount+=nSize;
//
//		nSize=strAirCep.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strAirCep.GetBuffer(), nSize);							// "[A]CEP[m]",
//		nCount+=nSize;
//
//		nSize=strAirMajorAxis.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strAirMajorAxis.GetBuffer(), nSize);				// "[A]����[m]",
//		nCount+=nSize;
//
//		nSize=strAirMinorAxis.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strAirMinorAxis.GetBuffer(), nSize);				//"[A]����(m)",
//		nCount+=nSize;
//
//		nSize=strAirMajorAxisAzimuth.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strAirMajorAxisAzimuth.GetBuffer(), nSize); //"[A]�������[��]",
//		nCount+=nSize;
//
//		nSize=strGndElintNotation.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strGndElintNotation.GetBuffer(), nSize);		// "[G]E-NOT",
//		nCount+=nSize;
//
//		nSize=strGndPinNr.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strGndPinNr.GetBuffer(), nSize);						// "[G]PIN NR",
//		nCount+=nSize;
//
//		nSize=strGndFunc.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strGndFunc.GetBuffer(), nSize);							// "[G]��ɺ�ȣ",
//		nCount+=nSize;
//
//		nSize=strGndSiteFunc.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strGndSiteFunc.GetBuffer(), nSize);					// "[G]����Ʈ���",
//		nCount+=nSize;
//
//		nSize=strGndThreat.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strGndThreat.GetBuffer(), nSize);						// "[G]����",
//		nCount+=nSize;
//
//		nSize=strLinkNum.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strGndEqupNo.GetBuffer(), nSize);						// "[G]����ȣ",
//		nCount+=nSize;
//
//		nSize=strGndConsistancy.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strGndConsistancy.GetBuffer(), nSize);			// "[G]��ġ��[%]",
//		nCount+=nSize;
//
//		nSize=strGndGeoLatitude.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strGndGeoLatitude.GetBuffer(), nSize);					//"[G]������ǥ",
//		nCount+=nSize;
//
//		nSize=strGndGeoLongitude.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strGndGeoLongitude.GetBuffer(), nSize);					//"[G]������ǥ",
//		nCount+=nSize;
//
//		nSize=strGndAlt.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strGndAlt.GetBuffer(), nSize);							// "[G]��[m]",
//		nCount+=nSize;
//
//		nSize=strGndConfIterval.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strGndConfIterval.GetBuffer(), nSize);			// "[G]�ŷڱ���[m]",
//		nCount+=nSize;
//
//		nSize=strGndCep.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strGndCep.GetBuffer(), nSize);							// "[G]CEP[m]",
//		nCount+=nSize;
//
//		nSize=strGndMajorAxis.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strGndMajorAxis.GetBuffer(), nSize);				// "[G]����[m]",
//		nCount+=nSize;
//
//		nSize=strGndMinorAxis.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strGndMinorAxis.GetBuffer(), nSize);				//"[G]����(m)",
//		nCount+=nSize;
//
//		nSize=strGndMajorAxisAzimuth.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strGndMajorAxisAzimuth.GetBuffer(), nSize); //"[G]�������[��]",			
//		nCount+=nSize;
//
//		nSize=strFinalReport.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strFinalReport.GetBuffer(), nSize);					// "��������",			
//		nCount+=nSize;
//		
//		nSize=strFinalAlarm.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strFinalAlarm.GetBuffer(), nSize);					// "�����溸",							",
//		nCount+=nSize;
//
//		nSize=strPosEstimation.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strPosEstimation.GetBuffer(), nSize);				// "���⿩��",							
//		nCount+=nSize;
//
//		nSize=strSOI.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strSOI.GetBuffer(), nSize);									// "���ɽ�ȣ"	
//		nCount+=nSize;
//						
//		nSize=strIsStorePDW.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strIsStorePDW.GetBuffer(), nSize);					// "PDW ���忩��"
//		nCount+=nSize;
//
//		nSize=strNumOfCol.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strNumOfCol.GetBuffer(), nSize);						// "��������[PDW/IQ]"
//		nCount+=nSize;
//
//		nSize=strNumOfAmbigiousBeam.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strNumOfAmbigiousBeam.GetBuffer(), nSize);	// "��ȣ��",
//		nCount+=nSize;
//
//		nSize=strNumOfPDWSet.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strNumOfPDWSet.GetBuffer(), nSize);					// "PDW ��Ʈ"
//		nCount+=nSize;
//
//		nSize=strActivatedOrNot.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strActivatedOrNot.GetBuffer(), nSize);					// "Ȱ������"
//		nCount+=nSize;
//
//		nSize=strLobId.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strLobId.GetBuffer(), nSize);					// "LOB ID"
//		nCount+=nSize;
//
//		nSize=strDtctId.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strDtctId.GetBuffer(), nSize);					// "DTCT ID"
//		nCount+=nSize;
//
//		nSize=strIdCode.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strIdCode.GetBuffer(), nSize);					// "�ĺ��ڵ�"
//		nCount+=nSize;
//
//		nSize=strPolization.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strPolization.GetBuffer(), nSize);			// "�ؼ�"
//		nCount+=nSize;
//
//		nSize=strTaskOperator.GetLength();
//		memcpy(&i_pbyOutData[nCount], &nSize, sizeof(USHORT));
//		nCount+=sizeof(USHORT);
//		memcpy(&i_pbyOutData[nCount], strPolization.GetBuffer(), nSize);			// "���� ����"
//		nCount+=nSize;
//
//		return nCount;
//	}
} I_EMITTER_STR;		

#define NUM_OF_LOB_LIST_FIELD_NAME 45
static const char* strLobListFieldName [] =
{
	"����",
	"LINK��ȣ", // 2015.1.5. ������ �߰�
	"E-NOT",
	"�ĺ��ڵ�",
	"���� ����",
	"����/�ð�",
	"LOB#",
	"AET��ȣ",
	"��Ʈ���޽� ��������",
	"���ļ�����",
	"���ļ�[MHz]",
	"�ּ� ���ļ�[MHz]",
	"�ִ� ���ļ�[MHz]",
	"PW[us]",
	"PRI[us]",
	"PRF[KHz]",
	"�׷���޽���",
	"PIT",
	"�����Ǽ�",
	"AMP",
	"RF Dev[MHz]",
	"DOA",
	"AOA",
	"AzStdDev",
	"AMB",
	"Sid",
	"BL",
	"FOV",
	"AC#",
	"AC��ġ",
	"AC Alt",
	"AC FOM",
	"Hdg",
	"Pitch",
	"Roll",
	"ǰ��",
	"PDW ����",
	"��ĵ����/�ֱ�[ms]",
	"PRI��������",
	"PRI������",
	"��������",
	"����ID",
	"��ȣID",
	"�б�",
	"�޽��� �׷�"
};

typedef struct stLOBInfoDisp
{
	bool bIsSelected;
	CString strElintNotation;	// "ELINT Notation"
	CString strLinkNum;	   //#FA_Q_2502_T1		// Link ��ȣ, 2015.1.5. ������ �߰�
	CString strSrchTime;			// 	"����/�ð�"
	CString strLobNum;			// "LOB#"
	CString strEmitterNum;		// "Emitter#"
	CString strModType;			// "��������"
	CString strFrqType;			// "���ļ�����",
	CString strFrq;					// "���ļ�[Hz]",
	CString strMinFrq;				// "�ּ� ���ļ�[Hz]",
	CString strMaxFrq;			// "�ִ� ���ļ�[Hz]",
	CString strPW;					// "PW",
	CString strPRI;					// "PRI",
	CString strPRF;					// "PRF",
	CString strMinPRF;			// "�ּ� PRF[PPS]",
	CString strMaxPRF;			// "�ִ� PRF[PPS]",
	CString strPulsePerGrp;		// "�׷���޽���",
	CString strPIT;					// "PIT",
	CString strNumOfPos;		// "�����Ǽ�",
	CString strAmp;				// "AMP",
	CString strRfDev;				// "RF Dev.",
	CString strDoa;				// "DOA",
	CString strAoa;				// "AOA",
	CString strAzStdDev;		// "AzStdDev",
	CString strAmb	;				// "AMB",
	CString strSid;					// "Sid",
	CString strBL;					// "BL",
	CString strFov;				// "FOV",
	CString strACNum;			// "AC#",
	CString strACPos;				// "AC��ġ",
	CString strACAlt;				//"AC Alt",
	CString strAcFom;			// "AC FOM",
	CString strHdg;				// "Hdg",
	CString strPitch;				// "Pitch",
	CString strRoll;					// "Roll",
	CString strQuality;			// "ǰ��",
	CString strNumOfPdw;		// "PDW ����",
	CString strScanTypePeriod;// "��ĵ����/�ֱ�",
	CString strPriModType;		// "PRI��������",
	CString strPriJitterRate;		//"PRI������"
	CString strFinalReport;		// "��������"
	CString strTaskId;				// "����ID" 2014.10.08 �߰�. ������.
	CString strTaskName;			// "������" 2015.06.03 �߰�. ��ö��
	CString strDtctId;				// DTCT ID . 2015.3.3. ������ �߰�
	CString strIdCode;				// �ĺ� �ڵ�, 2015.4.2. ��ö�� �߰�
	CString strTaskOperator;	// ���� �����, 2015.4.7. ��ö�� �߰�
	CString strPolization;	// �ؼ�. 2015.7.15. ������ �߰�
	CString strPPG;	// �޽��� �׷�. 2015.7.15. ������ �߰�

	stLOBInfoDisp()
		:bIsSelected(false)
	{
		strElintNotation="";
		strLinkNum="";			
		strSrchTime="";			
		strLobNum="";			
		strEmitterNum="";		
		strModType="";			
		strFrqType="";					
		strFrq="";							
		strMinFrq="";						
		strMaxFrq="";					
		strPW="";					
		strPRI="";					
		strPRF="";				
		strMinPRF="";			
		strMaxPRF="";			
		strPulsePerGrp="";		
		strPIT="";					
		strNumOfPos="";		
		strAmp="";				
		strRfDev="";				
		strDoa="";				
		strAoa="";				
		strAzStdDev="";		
		strAmb="";				
		strSid="";					
		strBL="";					
		strFov="";				
		strACNum="";			
		strACPos="";				
		strACAlt="";				
		strAcFom="";			
		strHdg="";				
		strPitch="";				
		strRoll="";					
		strQuality="";			
		strNumOfPdw="";		
		strScanTypePeriod="";
		strPriModType="";		
		strPriJitterRate="";		
		strFinalReport="";		
		strTaskId="";				
		strTaskName="";		
		strDtctId="";				
		strIdCode="";			
		strTaskOperator="";	
		strPolization="";
		strPPG="";
	};
}I_LOB_STR;	

// RAW DATA : PDW ������ Ÿ��
#define NUM_OF_PDW_HEAD_FIELD_NAME 9
static const char *strPdwHeadFieldName[] = {
	"����",
	"��ȣ",
	/*"PDW Set ID",*/
	"LINK��ȣ",
	"���� ID",
	"Ž���뿪��ȣ",
	/*"���ü ��ȣ",*/
	/*"�� ��ȣ",*/
	"LOB ��ȣ",
	"���ð���",
	"���� �޽� ��ȣ",
	"����"
};

#define NUM_OF_PDW_BODY_FIELD_NAME 19
static const char *strPdwBodyFieldName[] = {
	"�������۽ð�",
	"��� DTOA[us]",
	"��ȣ����",
	"�ؼ� ��ȿ��[%]",
	/*"FMOP[%]",
	"PMOP[%]",*/
	"Blanking[%]",
	//"ä�κ��濩��",
	/*"BLK[%]",*/
	"DV[%]",
	"FOV(IN)[%]",
	"Ch#",
	"��� ����[dBm]",
	"��� ���ļ�[MHz]",
	"��� ����[��]",
	"�ؼ�[%]",
	"PPF",
	"��� �޽���[ns]",
	//"I ������",
	//"Q ������",
	"PRF ID",
	"Spectrum ID",
	//"DtctId",
	//"���� ����",
	"���ϸ�",
	"������",
	"�� ���ϸ�"
};

struct I_PDWIQIF_STR
{
	CString strPdwSetId; // "PDW Set ID",
	CString strLinkNum; //#FA_Q_2502_T1     // "Link ��ȣ",
	CString strTaskId; // "���� ID",
	CString strSrchBandId; // "Ž���뿪��ȣ",
	CString strAetNum; // "AET ��ȣ"
	CString strAbtNum; //"ABT ��ȣ". 2016.6.8. ������ �߰�
	CString strLobNum; // "LOB��ȣ",
	CString strPdwCount; // "PDW����/IQ����",
	CString strPdwId; // PDW ID, PDW ���� ��ȣ
	CString strDataType; // "����", : PDW, IQ, IF
	CString strTime; // "�������۽ð�",
	CString strToa; // "TOA",
	CString strSignalType; // "��ȣ����",
	//CString strBitFlag; //"BIT Flag",
	CString strPolFlag; //"�ؼ� ��ȿ��",
	CString strFmopFlag; // "FMOP Flag",
	CString strPmopFlag; // "PMOP Flag",
	CString strBlankingTag; // "Blanking Tag",
	CString strChannelChangePOP; // "ä�κ��濩��",
	CString strBLK; // "ä�κ������",
	CString strDI; //"DI",
	CString strFovFlag; // "FOV Flag",
	CString strChannelNum; // "Ch#",
	CString strPa; // "��ȣ����",
	CString strFrq; // "�������ļ�",
	CString strSigDirection; // "��ȣ����",
	CString strPolarization; // "��ȣ�ؼ�",
	CString strPPFTag; // "PPF Tag",
	CString strPw; // "�޽���",
	CString strIData; // "I������"
	CString strQData; // "Q������",	
	CString strFilename;		// ���ϸ�CString strDtctId; // DTCT ID
	CString strTaskName;
	CString strFatherFilename;
	CString strDtctId; // "DTCT ID", 2015.3.3. ������ �߰�
		
	UINT nPRFDBSetID;	//ȭ��ǥ�þȵ�, PRF ���-DB�� �ĺ� ID�θ� ���
	UINT nSpectrumDBSetID;	//ȭ��ǥ�þȵ�, Spectrum ���-DB�� �ĺ� ID�θ� ���
	
	//@start_WJH
	I_PDWIQIF_STR()
	{	
	//	strPdwSetId= _T("");
	//	strLinkNum=_T("");// "Link ��ȣ",
	//	strTaskId=_T("");		 // "���� ID",
	//	strSrchBandId=_T(""); // "Ž���뿪��ȣ",
	//	strAetNum=_T(""); // "AET ��ȣ"
	//	strAbtNum=_T(""); // "ABT ��ȣ". 2016.6.8. ������ �߰�
	//	strLobNum=_T(""); // "LOB��ȣ",
	//	strPdwCount=_T(""); // "PDW����/IQ����",
	//	strPdwId=_T(""); // PDW ID, PDW ���� ��ȣ
	//	strDataType=_T(""); // "����", : PDW, IQ, IF
	//	strTime=_T(""); // "�������۽ð�",
	//	strToa=_T(""); // "TOA",
	//	strSignalType=_T(""); // "��ȣ����",
	//	strBitFlag=_T(""); //"BIT Flag",
	//	strFmopFlag=_T(""); // "FMOP Flag",
	//	strPmopFlag=_T(""); // "PMOP Flag",
	//	strBlankingTag=_T(""); // "Blanking Tag",
	//	strChannelChangePOP=_T(""); // "ä�κ��濩��",
	//	strBLK=_T(""); // "ä�κ������",
	//	strDI=_T(""); //"DI",
	//	strFovFlag=_T(""); // "FOV Flag",
	//	strChannelNum=_T(""); // "Ch#",
	//	strPa=_T(""); // "��ȣ����",
	//	strFrq=_T(""); // "�������ļ�",
	//	strSigDirection=_T(""); // "��ȣ����",
	//	strPolarization=_T(""); // "��ȣ�ؼ�",
	//	strPPFTag=_T(""); // "PPF Tag",
	//	strPw=_T(""); // "�޽���",
	//	strIData=_T(""); // "I������"
	//	strQData=_T(""); // "Q������",	
	//	strFilename=_T("");		// ���ϸ�		
	//	strTaskName=_T(""); // DTCT ID				
	//	strFatherFilename=_T("");
		strDtctId=_T(""); // DTCT ID		

		nPRFDBSetID = 0;
		nSpectrumDBSetID = 0;	
	}	  
	//@end_WJH
} ;




//=======> ���� ������ ����
						  
											
typedef struct STempThreatData
{
	unsigned short usSdfId;
	unsigned char szTaskId[LENGTH_OF_TASK_ID];
	unsigned short usSearchBandId;
	unsigned short usSerialNum;
	unsigned char ucAnalysisStatus;
	STempThreatData(){
		usSdfId=0;
		memset(&szTaskId, 0, LENGTH_OF_TASK_ID);
		usSearchBandId=0;
		usSerialNum=0;
		ucAnalysisStatus=0;
	}
}I_THREAT_DATA;
static const I_THREAT_DATA stInitThreatData;

typedef struct 
{
	unsigned short usAetId;
 	unsigned short usThreatId;
 	unsigned short usEmitterId;
 	unsigned short usBeamId;

	unsigned int uiOpInitID;
	unsigned char	aucTaskID[LENGTH_OF_TASK_ID];

}I_AET_DATA;

//static const I_AET_DATA g_stInitAETData;


static const char* strIdResult[] = 
{
	"-",
	"�⼺",
	"����",
	"�Ҹ�"
};
enum EnumIdResult
{
	E_EL_PE_UNK_ID=0,
	E_EL_OLD_ID=1,						// �⼺
	E_EL_NEW_ID,							// ����
	E_EL_UNK_ID								// �Ҹ�
};

struct STR_CEDEOBID_INFO {
	EnumIdResult eIdResult;					// 0=�⼺, 1=����, 2=�Ҹ�

	// CED�� ���̴� ��� �ε���
	int nCoRadarModeIndex;
	int nRadarModeIndex[MAX_IDCANDIDATE];

	// CED�� ���̴� �ε���
	//int nCoRadarIndex;
	//int nRadarIndex[MAX_IDCANDIDATE];

	// EOB�� ����/��ġ �ε���
	int nThreatIndex;
	int nDeviceIndex;

// 	STR_CEDEOBID_INFO::STR_CEDEOBID_INFO() : eIdResult(E_EL_OLD_ID), nCoRadarModeIndex(0), nCoRadarIndex(0), nThreatIndex(0), nDeviceIndex(0)
// 	{
// 		memset( nRadarModeIndex, 0, sizeof(nRadarModeIndex) );
// 		memset( nRadarIndex, 0, sizeof(nRadarIndex) );
// 	}

} ;

// ��ġ ���� �ϱ� ���� ����ü ����
#ifndef _STR_LOBS
#define _STR_LOBS
typedef struct {
	float fDoa;
	float fLatitude;
	float fLongitude;
	float fAltitude;

	int iCollectorID;

	unsigned int uiLOBID;				// Ŭ�����͸� �ϸ鼭 ������ �ϱ� ���� ���� �߰�

} STR_LOBS;
#endif

/**
 * @enum      E_BEAM_CODE
 * @brief     �� ���� �ڵ� ����
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 */
enum E_BEAM_CODE
{
	E_UNKNOWN_CODE=0,
	E_CREATE_FREQ_TYPE_CODE=1,
	E_CREATE_PRI_TYPE_CODE,
	E_CREATE_INTRA_TYPE_CODE,
	E_CREATE_SIGNAL_TYPE_CODE,

	E_CREATE_DIST_CODE,

};

// LOB�� ������Ʈ�ϱ� ���� ����ü ����
typedef struct {
	unsigned int uiAETID;
	unsigned int uiABTID;
	unsigned int uiLOBID;
	string strAcqTime;

	unsigned int uiSeqNum;

} STR_UPDATE_LOB_ABT_AET ;

enum ENUM_PE_STAT {
	E_EL_PESTAT_FAIL=0,					// ��ġ ���� �� ������ �� ����
	E_EL_PESTAT_SUCCESS,				// ��ġ ���� ����� ������ ����
	E_EL_PESTAT_NOT_YET,				// LOB�� ��� ��ġ ���� ���� ���� ����
	E_EL_PESTAT_IMPOSSIBILITY		// �װ��Ⱑ ���� ��ġ�� �����ϱ� ����� ����
};

struct STR_POSITION_ESTIMATION {
	ENUM_PE_STAT enValid;
	float fLongitude;							// [deg]
	float fLatitude;							// [deg]
	//int iAltidude;							// [m]
	float fCEP;										// [m]
	float fMajorAxis;							// [m]
	float fMinorAxis;							// [m]
	float fTheta;									// [0.1��]
	float fDistanceErrorOfThreat;	// [m]

	//int iManualLongitude;				// [deg], MANUAL_POS_EST_LONG
	//int iManualLatitude;				// [deg], MANUAL_POS_EST_LAT

	bool bApplayOfLOBClustering;// TRUE: ����, FALSE: ������

// 	STR_POSITION_ESTIMATION::STR_POSITION_ESTIMATION() : enValid(E_EL_PESTAT_FAIL), fLongitude(0.0), fLatitude(0.0), fCEP(0), fMajorAxis(0), fMinorAxis(0), fTheta(0), fDistanceErrorOfThreat(0), bApplayOfLOBClustering(false)
// 	{
// 
// 	}

} ;


typedef struct STempAmbiBeam
{
	unsigned short usAmbiguousBeamId;
	STempAmbiBeam(){
		usAmbiguousBeamId=0;
	}
}I_AmbiBeam;


/*--------------------------------------------------------------------------------
*						GMI LOB ó���� ������ Ÿ�� �ű� ���� 2016.2.13. ������
*--------------------------------------------------------------------------------*/
/*! Beam �Ǵ� Emitter ���¸� ������ enum
 * @enum      E_BEAM_EMITTER_STAT
 * @brief 
 * @author    ������ (jeongnam.lee@lignex1.com)
 * @date      2016-02-13 ���� 8:27
 */
enum E_BEAM_EMITTER_STAT
{
	E_ES_NOT_AVAILABLE = 0,			// ��ü�� �������� ���� ��� (�ʿ伺�� ���ؼ��� ���� �ʿ�)

	// ���������� ���� ���� ���� �� : ������ �ʿ䰡 ����.
	E_ES_NEW,										// �ű� ó�� (���������� Ȱ����, ����), Ȱ��
	E_ES_UPDATE,								// Ȱ��
	E_ES_DELETE,								// ��Ȱ�� ó��
	E_ES_DEACTIVATED,						// Ȱ�� ���� 
	E_ES_REACTIVATED,						// Ȱ�� �簳

};

static const char* strBeamEmitterStat[] = 
{
	"-",

	"Ȱ��(�ű�)",
	"Ȱ��",
	"��Ȱ��",
	"����",
	"Ȱ��(�簳)",

};

// ������� ���� POSN���� OPCODE ����
enum E_EMITTER_OPCODE
{
	E_EO_NOT_AVAILABLE = 0,					// ��ü�� �������� ���� ��� (�ʿ伺�� ���ؼ��� ���� �ʿ�)

	// �Ʒ� �ڵ���� �ش� ���ü/�� ���â�� �����ؾ� ��.
	E_MR_UPDATE_INFO,								// ���� ó�� (ELNOT, ��ġ ���� �� ����), ���ü/�� ���̺� ������(ELNOT, ��ġ ���� ����) �߰�
	E_MR_REMOVE_AETABT,							// ���ü/�� ���â���� ����
	E_MR_REMOVE_LOB,								// LOB ���â���� ����

	// �Ʒ� �ڵ���� �ش� LOB ���â�� �����ؾ� ��.
	E_MR_CHANGE,										// (���ü/)��/LOB ��ȣ ����

	//////////////////////////////////////////////////////////////////////////
	E_MR_UPDATE_STAT,								// ���� ���� ����, (Ȱ��, ��Ȱ�� ��)

	E_MR_ALERT_UPDATE,							// �溸 ó��
	E_MR_REPORT_UPDATE,							// ���� ó��

	E_MR_REMOVE_AET									// ���ü���� �ش� �׸� ����

};

static const char* strBeamEmitterOpcode[] = 
{
 	"���������ȵ�",

	"���� ",
	"����(E/B)",
	"����(L) ",

	"���� ",
	//////////////////////////////////////////////////////////////////////////
	"���� ",

	"�溸 ",
	"����",

	"����(E)",

	"���� ��������(NoChange)"
};


typedef struct {
	int uiSeqNum;

	//////////////////////////////////////////////////////////////////////////
	// LOB ������ �ʿ��� ����
	int iLinkNum;																								// ������ ���� �������� ��ũ ��ȣ

	bool isFiltered;																						// ���� ���� ����
	bool isManualEdited;																				// ���� ���� ����

	unsigned int uiAETID;
	unsigned int uiABTID;
	unsigned int uiLOBID;

	__time32_t tiAcqTime;																						// �װ����� �м� �ð�, LOB �޽����� time_t �� �Ǿ� ������ ������.
	int tiContactTimems;																				// �װ����� �м� �ð�

	E_BEAM_EMITTER_STAT enEmitterStat;													// LOB ����

	int iTaskType;																							// ���� ���� ����

	unsigned int uiCoLOB;																				// LOB ����

	int iBeamValidity;

	//////////////////////////////////////////////////////////////////////////
	// �ĺ� �߰� ����
	//
	// �ĺ� �ĺ� ����
	bool bOverCount;
	//int usCoCandidate;
	int usThreatId[MAX_IDCANDIDATE];

	// ���� ����
	//unsigned short usPriorityLevel;

	// ��ġ��
	int iIdRatio[MAX_IDCANDIDATE];

	// �ű� CED �� EOB�� �ĺ��� ����
	STR_CEDEOBID_INFO idInfo;

	// ��ġ ���� ����
	STR_POSITION_ESTIMATION peInfo;

	//char chELNOT[_MAX_ELNOT_STRING_SIZE_];					// MAX_SIZE_OF_ELNOT

	// EOB �Ÿ� ���� [km]
	double dEOBErrorDistance;

	// �ĺ� ���
	//EnumLibType sIdObject;

	// ���� �ڵ�
	//char modCode[MAX_SIZE_OF_MODULATIONCODE];

	E_BEAM_CODE eBeamCode;

	bool bCompFreq;
	bool bCompPRI;

} I_AET_ANAL;

//static I_AET_ANAL g_stInitAETAnal;





typedef struct STempAET
{
	I_THREAT_DATA stThreat;
	I_AET_DATA stAet;
	I_AET_ANAL stAnal;
	STempAET(){
		//memset(&stThreat, 0, sizeof(I_THREAT_DATA));
		//memset(&stAet, 0, sizeof(I_AET_DATA));
		//memset(&stAnal, 0, sizeof(I_AET_ANAL));

	/*	g_stInitThreatData;
		g_stInitAETData;
		g_stInitAETAnal;*/
	}
}I_AET;
static const I_AET g_stInitAET;

/**
 * @typedef   SELLOBDATA_EXT
 * @brief     ���� ���� ������ ����ü �̿ܿ� �߰��� ���� ����ü ����
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 */
typedef struct {
	I_AET_ANAL aetAnal;
	I_AET_DATA aetData;

} SELLOBDATA_EXT;


typedef struct {
	SRxLOBData stLOBData;
	SELLOBDATA_EXT stLOBDataExt;

} STR_LOBDATAEXT ;

typedef struct {
	int iSignalType;
	int iFrqType;
	int iPRIType;
	//int iScanType;
	//int iMOPType;

} STR_ID_TYPE;

/**
 * @typedef   SELABTDATA_EXT
 * @brief     ABT ����ü ����
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 */
typedef struct {
	// ������ ��ü PDW ����
	//int nCoTotalPdw;
	//int nCoTotalIQ;

	// ��Ʈ�� ���� ����
	//bool bIntraMop;

	// ��ġ ���� ����-Covariace ��
	double dPECoVar[4];

	float fLastAOA;
	double dRadarCollectionLatitude;
	double dRadarCollectionLongitude;

	// ���� ���� ����
	//bool bIsManualEdited;

	// ���ü�� Ȱ�� ���¸� ����
	enum E_BEAM_EMITTER_STAT enBeamEmitterStat;

	// �� ��ȿ�� ����
	int nCoBeamValidity;

	// ���� ���� ���� �÷���
	//UELMANUALVAL xManualEdited;

	// FISINT�� ����
	//bool bIsFISINTTask;
	UINT uiOpInitID;

	// ����/�ڵ� ��ġ ��� ���� ����
	//bool bManualPosEstPreferred;
	unsigned int uiSeqNum;

	//SLOBOtherInfo stLOBOtherInfo;

	int iLOBPoolIndex;

	STR_ID_TYPE stIDType;

	//enTHREAT_PLATFORM enPlatform;

	int nCoIdEOB;
	//STR_EOB_RESULT stEOBResult[MAX_CANDIDATE_EOB];

	// �� ���� ����
	bool bCompFreq;
	bool bCompPRI;

	// �ű� CED �� EOB�� �ĺ��� ����
	STR_CEDEOBID_INFO idInfo;

	ENUM_PE_STAT enValid;
	bool bApplayOfLOBClustering;// TRUE: ����, FALSE: ������

	UINT uiPE;
	bool bFullOfPE;
	double dBLatitude[MAX_OF_LOBS_PE];
	double dBLongitude[MAX_OF_LOBS_PE];
	double dEasting[MAX_OF_LOBS_PE];
	double dNorthing[MAX_OF_LOBS_PE];

} SELABTDATA_EXT;

typedef struct {
	// ������ ��ü PDW ����
	int nCoTotalPdw;
	//int nCoTotalIQ;

	// ��Ʈ�� ���� ����
	// bool bIntraMop;
	
	// ���ü ���� �ĺ��� ���� ����
	//char szEOBELNOT[_MAX_ELNOT_STRING_SIZE_];					// EOB ���� ELNOT

	// ���ü�� Ȱ�� ���¸� ����
	enum E_BEAM_EMITTER_STAT enBeamEmitterStat;

	// ���� ���� ���� �÷���
	bool bIsManualEdited;

	bool bIsFISINTTask;

	// ���� ���� ����
	//UELMANUALVAL xMannualEdited;

	int iPinNumber;
	//char szELNOT[_MAX_ELNOT_STRING_SIZE_];

	//int iTaskType;		// ���� ����

	// �������� �̽ĺ� �� �߿��� ��ǥ ABT ��ȣ
	unsigned int nUnIDABTID;

	// ����/�ڵ� ��ġ ��� ���� ����
	bool bManualPosEstPreferred;
	unsigned int uiSeqNum;


} SELAETDATA_EXT;

typedef struct {
	SRxLOBData stLOBData;
	SRxABTData stABTData;
	SELABTDATA_EXT stABTDataExt;

} STR_ABTDATAEXT ;


/*!
 * @typedef   SELIDENTIFICATION_INFO
 * @brief     ��ȣ �ĺ� �ɼ� â�� ���� ����ü ���� (GMI ���� ó���� ����ü)
 * @author    ������(jeongnam.lee@lignex1.com)
 * @date      2016-02-13 
 */
// typedef struct {
// 	SRxThreatData				stMsgData;		// �װ����� ���ŵ� �޽����� Data (Header) ����
// 	SRxThreatDataGroup		stMsgDataGrp;	// �װ����� ���ŵ� �޽����� DataGrp ����
// 	SPosEstData stPosEst;							// ���󿡼� ��ġ ������ ����
// 	I_AET stCDFAet;									// ��ȣ �ĺ� �Է� �� ��� ����
// 
// } SELIDENTIFICATION_INFO_MR;

/*! Beam (ABT) ������ ��� �ִ� ����ü
 * @struct     SBeamInfoFamily
 * @brief 
 * @author    ������ (jeongnam.lee@lignex1.com)
 * @date      2016-02-13 ���� 8:28 
 */
#define LENGTH_OF_GMI_CHAR 64
struct SThreatFamilyInfo
{
	E_EMITTER_OPCODE enOpcode;						// ������� POSN���� ��� ��

	unsigned int nSeqNum;									// DB ���̺��� SEQ_NUM

	unsigned int iAETID;									// ���ü #
	unsigned int iABTID;									// ABT #
	unsigned int iLOBID;									// LOB #

	bool bApplySearchFilter;							// ApplySearchFilterToAlarmAndMapDisplay() �Լ� ���� ����

	bool bIsFISINTTask;										// FISINT�� ����

	// �溸 �� ���� ����/���� ���� �ð�
	time_t ti_FirstTime;									// ���� �ð�
	time_t ti_FinalTime;								// ���� �ð�(�溸)

	E_BEAM_CODE eBeamCode;								// �� ���� ����, ������� ���� �ڵ�

	E_BEAM_EMITTER_STAT enEmitterStat;					// ���� ���� ����

	unsigned int iChangedAETID;							// AET #
	unsigned int iChangedABTID;							// ABT #

	SThreatFamilyInfo() : enOpcode(E_EO_NOT_AVAILABLE), nSeqNum(0), iAETID(0), iABTID(0), iLOBID(0), bApplySearchFilter(false), bIsFISINTTask(false), eBeamCode(E_UNKNOWN_CODE), enEmitterStat(E_ES_NOT_AVAILABLE), iChangedABTID(0), iChangedAETID(0)
	{

	}

} ;



typedef struct SEmitterFilter
{
	bool bFrqUse;
	float fFrqMin;
	float fFrqMax;
	bool bPriUse;
	float fPriMin;
	float fPriMax;
	bool bPwUse;
	float fPwMin;
	float fPwMax;
	SEmitterFilter()
		:bFrqUse(false)
		,bPriUse(false)
		,bPwUse(false)
		,fFrqMin(0.0)
		,fFrqMax(0.0)
		,fPriMin(0.0)
		,fPriMax(0.0)
		,fPwMin(0.0)
		,fPwMax(0.0)
	{
	}
}I_EMITTER_FILTER;

// typedef struct STempPDW
// {
// 	SRxPDWData stData;
// 	SAvgPDWDataRGroup stAvgData;
// 	//std::list<SRXPDWDataRGroup> listGroup; // 2014.07.15.������. ������� �ʾƼ� �ּ�ó��
// 	STempPDW()
// 	{
// 		memset(&stData, 0, sizeof(SRxPDWData));
// 		memset(&stAvgData, 0, sizeof(SRXPDWDataRGroup));
// 		//stData.Init();
// 		//stAvgData.Init();
// 		//listGroup.clear(); // 2014.07.15.������. ������� �ʾƼ� �ּ�ó��
// 	}
// }I_PDW;

// ���ü/LOB ���â���� ���� �ۼ� �޴� Ŭ�� �� �߻�
// struct SEmitterInfoForReport
// {
// 	std::string strDtctId;
// 	std::string strTaskId;
// 	int nAetId;
// 	int nLobId;
// 	bool bFromEmitterWnd;
// 	SEmitterInfoForReport()
// 	: strTaskId("")
// 	, strDtctId("")
// 	, nAetId(-1)
// 	, nLobId(-1)
// 	, bFromEmitterWnd(true){};
// };

struct SEmitterActivationInfo
{
	int nDetectYear; //���� �⵵
	int nDetectMon; // ���� ��
	int nDetectDay; // ���� ��
	int nDetectTimeHour;	// ���� �ð� (��)
	int nDetectTimeMin;		// ���� �ð� (��)
	int nDetectTimeSec;		// ���� �ð� (��)
	int nRecentUpdateTimeHour;  // ���� �ð� (��)
	int nRecentUpdateTimeMin;	 // ���� �ð� (��)
	int nRecentUpdateTimeSec;	 // ���� �ð� (��)
	int nTotalActTimeHour;			// �� Ȱ�� �ð� (��)
	int nTotalActTimeMin;			 // �� Ȱ�� �ð� (��)
	int nTotalActTimeSec;			// �� Ȱ�� �ð� (��)
	bool bDeactivatedNow;			// Ȱ������ ����
	int nDeactivatedTimeHour;			// Ȱ������ �ð� (��)
	int nDeactivatedTimeMin;			 // Ȱ������ �ð� (��)
	int nDeactivatedTimeSec;			// Ȱ������ �ð� (��)
	bool bTimeOfConvergenceIsValid; // 
	int nTimeOfConvergenceHour;
	int nTimeOfConvergenceMin;
	int nTimeOfConvergenceSec;
	int nTimeOfConvergenceLobNum;
	int nNumOfDetectionToday;
	int nNumOfTotalLob;
	int nNumOfLobUsedToPosEst;
	int nNumOfTotalPdw;
	bool bIsSoi; // ���ɽ�ȣ����. 2015.3.30. ������. ���ɽ�ȣ�� �溸/�������� ���� ��ȣ������ �ǹ���.
	int nFinalAircraftLatitude;
	int nFinalAircraftLongitude;
	int nFinalAircraftDoa;
	SEmitterActivationInfo()
		:nDetectYear(0),
		 nDetectMon(0),
	     nDetectDay(0),
		 nDetectTimeHour(0),
		 nDetectTimeMin(0),
		 nDetectTimeSec(0),
		 nRecentUpdateTimeHour(0),
		 nRecentUpdateTimeMin(0),
		 nRecentUpdateTimeSec(0),
		 nTotalActTimeHour(0),
		 nTotalActTimeMin(0),
		 nTotalActTimeSec(0),
		 bDeactivatedNow(false),
		 nDeactivatedTimeHour(0),
		 nDeactivatedTimeMin(0),
		 nDeactivatedTimeSec(0),
		 bTimeOfConvergenceIsValid(false),
		 nTimeOfConvergenceHour(0),
		 nTimeOfConvergenceMin(0),
		 nTimeOfConvergenceSec(0),
		 nTimeOfConvergenceLobNum(0),
		 nNumOfDetectionToday(0),
		 nNumOfTotalLob(0),
		 nNumOfLobUsedToPosEst(0),
		 nNumOfTotalPdw(0),
		 bIsSoi(false),
		 nFinalAircraftLatitude(0),
		 nFinalAircraftLongitude(0),
         nFinalAircraftDoa(0)
	{};
};

typedef struct stIdntfyReslt
{
	string			strIdnfyRsltID;				// "0" �� �ڵ����� DB���� Unique ���� �Ҵ�
	string			strDtctSigID;						// ���� : ����.
	string			strIdnfyTime;					// ���� : ��:��:��
	int	nRadarModeIndex;							// nRadarModeIndex
	int	nThreatIndex;	
	int	nDeviceIndex;
	unsigned int uiBeamID;
	int uiLobID;
	int uiAetID;
	unsigned int	uiCandidateNum;
	string			strMsgAetID;	
	string			strElnotAir;
	string			strElnotGnd;
	stIdntfyReslt()
		:uiBeamID(0)
		,uiLobID(0)
		,uiAetID(0)
		,nRadarModeIndex(0)
		,nThreatIndex(0)
		,nDeviceIndex(0)
		,uiCandidateNum(0)
		,strIdnfyRsltID("")		
		,strDtctSigID("")			
		,strIdnfyTime("")
		,strMsgAetID("")	
		,strElnotAir("")
		,strElnotGnd("")
	{		
	};
}SELIdentify;

/**
 * @typedef   I_PDW_IQ_DATA
 * @brief     PDW �� IQ ������ ���� ���θ� �����ϴ� ����ü
 * @author    ��ö�� (churlhee.jo@lignex1.com)
 */
typedef struct SPDWIQData
{
	bool bPDW;
	bool bIQ;
	SPDWIQData(){
		bIQ = false;
		bPDW = false;
	}
} I_PDW_IQ_DATA;

// ���ü/��/LOB ���� ����
#define _ACTIVE_STAT_						"Ȱ��"
#define _DEACTIVE_STAT_				"����"
#define _REACTIVE_STAT_					"�簳"
#define _DELETE_STAT_					"��Ȱ��"

/*! Emitter ������ ��� �ִ� ����ü
 * @struct     SEmitterInfoFamily
 * @brief 
 * @author    ������ (jeongnam.lee@lignex1.com)
 * @date      2016-02-13 ���� 8:28 
 */
// struct SEmitterInfoFamily
// {
// 	E_BEAM_EMITTER_STAT eStat;					// Emitter Status 
// 	char szEmitterIdGnd[LENGTH_OF_GMI_CHAR];							// ���󿡼� ������ Emitter Unique ID
// 	char szPastEmitterIdGnd[LENGTH_OF_GMI_CHAR];					// Emitter Status�� E_MR_NEW_BY_DEMERGE �� ���, ���� Emitter ID�� ����.
// 	SELIDENTIFICATION_INFO_MR stEmitter;	// �ĺ�, ��ġ����, ���� �� Emitter ����	
// 	
// 	SEmitterInfoFamily(){
// 		memset(szEmitterIdGnd, 0, LENGTH_OF_GMI_CHAR);
// 		memset(szPastEmitterIdGnd, 0, LENGTH_OF_GMI_CHAR);
// 	}
// };

/*! Emitter �Ǵ� Beam�� Merge�� ���Ͽ� Delete ����� ���� ���, ID�� ����ϱ� ���� ����ü
 * @struct     SDeleteID
 * @brief 
 * @author    ������ (jeongnam.lee@lignex1.com)
 * @date      2016-02-13 ���� 8:40 
 */
struct SDeleteID
{
	char szIdGnd[LENGTH_OF_GMI_CHAR];					// ���� ����� ID. (Emitter  Unique ID �Ǵ� Beam Unique ID)
	
	SDeleteID(){
		memset(szIdGnd, 0, LENGTH_OF_GMI_CHAR);	
	}
} ;

// VER.3 ���â�� ����ü && ENUM Ÿ��
enum E_GMI_PROC_MSG_INFO_TYPE 
{	
	E_GMI_LOB = 0,
	E_GMI_ABT,
	E_GMI_AET,
};

// #define LENGTH_OF_SHORT_CHAR 64
// #define LENGTH_OF_FRQ_OR_PRI_SEQ_CHAR 1024
// struct SLobInfoToWnd
// {
// 	unsigned long long ullSeqNum;
// 	char szMissionId[LENGTH_OF_MISSION_ID+1];
// 	char szMissionName[SIZE_OF_TASK_NAME];
// 	char szTaskId[LENGTH_OF_TASK_ID];
// 	char szTaskName[SIZE_OF_TASK_NAME];
// 	int nSearchBandId;
// 	char szTaskType[LENGTH_OF_SHORT_CHAR];
// 	char szDetAntDirection[LENGTH_OF_SHORT_CHAR];
// 	char szLinkNum[LENGTH_OF_SHORT_CHAR];
// 	char szAircraftNum[LENGTH_OF_SHORT_CHAR];
// 	char szRxPath[LENGTH_OF_SHORT_CHAR];
// 	unsigned int uiLobId;
// 	unsigned int uiAbtId;
// 	unsigned int uiAetId;
// 	char szElnotPri[LENGTH_OF_SHORT_CHAR];
// 	char szModeCodePri[LENGTH_OF_SHORT_CHAR];
// 	char szModulationCode[LENGTH_OF_SHORT_CHAR];
// 	char szNickName[LENGTH_OF_SHORT_CHAR];
// 	char szPriFuncCode[LENGTH_OF_SHORT_CHAR];
// 	char szElnotSec[LENGTH_OF_SHORT_CHAR];
// 	char szModeCodeSec[LENGTH_OF_SHORT_CHAR];
// 	char szElnotTert[LENGTH_OF_SHORT_CHAR];
// 	char szModeCodeTert[LENGTH_OF_SHORT_CHAR];
// 	char szSigType[LENGTH_OF_SHORT_CHAR];
// 	int iPolarization;
// 	char szPolarization[LENGTH_OF_SHORT_CHAR];
// 	int iRatioOfPOL;
// 	//int iIsFISINTTask;
// 	char szIsFISINTTask[LENGTH_OF_SHORT_CHAR];
// 	
// 	int nNumOfPPG;
// 	int nNumOfPulse;
// 	char szTimeInfo[LENGTH_OF_SHORT_CHAR];
// 	char szValidity[LENGTH_OF_SHORT_CHAR];
// 	char szBL[LENGTH_OF_SHORT_CHAR];
// 	int iFov;
// 	int iDoaMean; // �ػ� 0.1 deg
// 	int iDoaMax;
// 	int iDoaMin;
// 	int iDoaDev;
// 	int iDoaStd;
// 	int iDIRatio;
// 	char szFrqType[LENGTH_OF_SHORT_CHAR];
// 	char szFrqPatternType[LENGTH_OF_SHORT_CHAR];
// 	char szHasFrqPeriod[LENGTH_OF_SHORT_CHAR];
// 	int nFrqPeriod; // �ػ� 10khz
// 	int nFrqPositionCount;
// 	int nFrqElementCount;
// 	int nFrqMean; // �ػ� 10khz
// 	int nFrqMax;
// 	int nFrqMin;
// 	int nFrqDev;
// 	int arrFrqSeq[MAX_FREQ_PRI_STEP];
// 	char szFrqSeq[LENGTH_OF_FRQ_OR_PRI_SEQ_CHAR];	
// 	char szPriType[LENGTH_OF_SHORT_CHAR];
// 	char szPriPatternType[LENGTH_OF_SHORT_CHAR];
// 	char szHasPriPeriod[LENGTH_OF_SHORT_CHAR];
// 	int nPriPeriod; // us
// 	int nPriPositionCount;
// 	int nPriElementCount;
// 	int nPriMean;
// 	int nPriMax;
// 	int nPriMin;
// 	int nPriDev;
// 	int iPriJitterRatio;
// 	int arrPRISeq[MAX_FREQ_PRI_STEP];
// 	char szPRISeq[LENGTH_OF_FRQ_OR_PRI_SEQ_CHAR];
// 	int nPrfMean;
// 	int nPrfMax;
// 	int nPrfMin;
// 	int nPrfDev;
// 	int iPaMean;
// 	int iPaMax;
// 	int iPaMin;
// 	int iPaDev;
// 	int iPwMean;
// 	int iPwMax;
// 	int iPwMin;
// 	int iPwDev;
// 	char szScanType[LENGTH_OF_SHORT_CHAR];
// 	char szScanTypeDetail[LENGTH_OF_SHORT_CHAR];
// 	int nScanPeriodMicroSec;
// 	float fScanPeriodHz;
// 	char szIntraType[LENGTH_OF_SHORT_CHAR];
// 	char szIntraDetailType[LENGTH_OF_SHORT_CHAR];
// 	int nIntraFrqMean; // 10khz�ػ�
// 	int nIntraFrqMax;
// 	int nIntraFrqMin;
// 	int nIntraFrqChangeWidth;
// 	int arrPriPerGroup[MAX_PRI_PER_GROUP];
// 	int arrPaDiffPerGroup[MAX_PADIFF_PER_GROUP];
// 	int nAircraftLatitude;
// 	int nAircraftLongitude;
// 	int nAircraftPitch; // 0.1�� �ػ�
// 	int nAircraftRoll;
// 	int nAircraftHeading;
// 	int nAircraftAltitude;
// 	int nAircraftFOM;
// 	char szIsPDWRestored[LENGTH_OF_SHORT_CHAR];
// 	char szIsIQRestored[LENGTH_OF_SHORT_CHAR];
// 	char szFirstReportTime[LENGTH_OF_SHORT_CHAR];
// 	char szFinalReportTime[LENGTH_OF_SHORT_CHAR];
// 	char szStat[LENGTH_OF_SHORT_CHAR];
// 	int nHour;
// 	int nMin;
// 	int nSec;
// 	int nIsSelect;
// 	unsigned int nAirLobId;
// 	int nIsFiltered;
// 	unsigned int uiNewAbtId;
// 	unsigned int uiNewAetId;
// 	E_EMITTER_OPCODE eOpCodeForUpdate; // ������Ʈ �Ǵ� ������ �پ��� ����, Update Key�� �ּ� �ĺ��ϱ�� ��. ������. 2018.4.10.
// 
// 	// ��ö��, �ð� ������ �����ϱ� ���� �߰��� �Լ�
// 	bool operator<(const SLobInfoToWnd & other ) 
// 	{
// 		return strcmp( this->szTimeInfo, other.szTimeInfo ) < 0;
// 	}
// 
// 	SLobInfoToWnd()
// 		:ullSeqNum(0),
// 		nSearchBandId(0),
// 	uiLobId(0),
// 	uiAbtId(0),
// 	uiAetId(0),
// 	iPolarization(0),
// 	iRatioOfPOL(0),
// 	nNumOfPPG(0),
// 	nNumOfPulse(0),
// 	iDoaMean(0), // �ػ� 0.1 deg
// 	iDoaMax(0),
// 	iDoaMin(0),
// 	iDoaDev(0),
// 	iDoaStd(0),
// 	iDIRatio(0),
// 	iFov(0),
// 	nFrqPeriod(0), // �ػ� 10khz
// 	nFrqPositionCount(0),
// 	nFrqElementCount(0),
// 	nFrqMean(0), // �ػ� 10khz
// 	nFrqMax(0),
// 	nFrqMin(0),
// 	nFrqDev(0),
// 	nPriPeriod(0), // us
// 	nPriPositionCount(0),
// 	nPriElementCount(0),
// 	nPriMean(0),
// 	nPriMax(0),
// 	nPriMin(0),
// 	nPriDev(0),
// 	iPriJitterRatio(0),
// 	iPwMean(0),
// 	iPwMax(0),
// 	iPwMin(0),
// 	iPwDev(0),
// 	nPrfMean(0),
// 	nPrfMax(0),
// 	nPrfMin(0),
// 	nPrfDev(0),
// 	iPaMean(0),
// 	iPaMax(0),
// 	iPaMin(0),
// 	iPaDev(0),
// 	nScanPeriodMicroSec(0),
// 	fScanPeriodHz(0.0),
// 	nIntraFrqMean(0), // 10khz�ػ�
// 	nIntraFrqMax(0),
// 	nIntraFrqMin(0),
// 	nIntraFrqChangeWidth(0),
// 	nAircraftLatitude(0),
// 	nAircraftLongitude(0),
// 	nAircraftPitch(0), // 0.1�� �ػ�
// 	nAircraftRoll(0),
// 	nAircraftHeading(0),
// 	nAircraftAltitude(0),
// 	nAircraftFOM(0),
// 	nHour(0),
// 	nMin(0),
// 	nSec(0),
// 	nIsSelect(0),
// 	nAirLobId(0),
// 	nIsFiltered(0),
// 	uiNewAbtId(0),
// 	uiNewAetId(0),
// 	eOpCodeForUpdate(E_EO_NOT_AVAILABLE)
// 	{
// 		memset(szMissionId, 0, LENGTH_OF_MISSION_ID+1);
// 		memset(szMissionName, 0,SIZE_OF_TASK_NAME);
// 		memset(szTaskId, 0,LENGTH_OF_TASK_ID);
// 		memset(szTaskName, 0,SIZE_OF_TASK_NAME);
// 		memset(szTaskType, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szDetAntDirection, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szLinkNum, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szAircraftNum, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szElnotPri, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szModeCodePri, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szModulationCode, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szNickName, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szPriFuncCode, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szElnotSec, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szModeCodeSec, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szElnotTert, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szModeCodeTert, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szSigType, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szTimeInfo, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szValidity, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szBL, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szFrqType, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szFrqPatternType, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szHasFrqPeriod, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szPriType, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szPriPatternType, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szHasPriPeriod, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szScanType, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szScanTypeDetail, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szIntraType, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szIntraDetailType, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szIsPDWRestored, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szIsIQRestored, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szFirstReportTime, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szFinalReportTime, 0,LENGTH_OF_SHORT_CHAR);	
// 		memset(szStat, 0, LENGTH_OF_SHORT_CHAR);
// 		memset(arrFrqSeq, 0, MAX_FREQ_PRI_STEP*sizeof(int));
// 		memset(arrPRISeq, 0, MAX_FREQ_PRI_STEP*sizeof(int));
// 		memset(arrPriPerGroup, 0, MAX_PRI_PER_GROUP*sizeof(int));
// 		memset(arrPaDiffPerGroup, 0, MAX_PADIFF_PER_GROUP*sizeof(int));
// 		memset(szIsFISINTTask, 0, LENGTH_OF_SHORT_CHAR);
// 		memset(szPRISeq, 0, LENGTH_OF_FRQ_OR_PRI_SEQ_CHAR);
// 		memset(szFrqSeq, 0, LENGTH_OF_FRQ_OR_PRI_SEQ_CHAR);
// 		memset(szPolarization, 0, LENGTH_OF_SHORT_CHAR);
// 		memset(szRxPath, 0, LENGTH_OF_SHORT_CHAR);//LDRA�߰�
// 
// 	}
// };
// 
// struct SAbtInfoToWnd
// {
// 	unsigned long long ullSeqNum;
// 	char szMissionId[LENGTH_OF_MISSION_ID+1];
// 	char szMissionName[SIZE_OF_TASK_NAME];
// 	unsigned int uiAbtId;
// 	unsigned int uiAetId;
// 	int iRadarModeIndex;
// 	int iThreatIndex;
// 	char szTaskType[LENGTH_OF_SHORT_CHAR];
// 	char szElnotPri[LENGTH_OF_SHORT_CHAR];
// 	char szModeCodePri[LENGTH_OF_SHORT_CHAR];
// 	char szModulationCode[LENGTH_OF_SHORT_CHAR];
// 	char szNickName[LENGTH_OF_SHORT_CHAR];
// 	char szPriFuncCode[LENGTH_OF_SHORT_CHAR];
// 	char szPlatform[LENGTH_OF_SHORT_CHAR];
// 	char szPlaceNameKor[LENGTH_OF_SHORT_CHAR]; // 2018-02-22 �ѱ����� �߰�
// 	int nRadarModePriority;
// 	float fDistanceErrFromThreat;
// 	int nRadarPriority;
// 	char szSigType[LENGTH_OF_SHORT_CHAR];
// 	int nNumOfPPG;
// 	int nNumOfLOB;
// 	int iDoaMean; // �ػ� 0.1 deg
// 	int iDoaMax;
// 	int iDoaMin;
// 	int iDoaDev;
// 	int iDoaStd;
// 	char szFirstDetectTimeInfo[LENGTH_OF_SHORT_CHAR];
// 	char szFinalDetectTimeInfo[LENGTH_OF_SHORT_CHAR];
// 	char szValidity[LENGTH_OF_SHORT_CHAR];
// 	char szFrqType[LENGTH_OF_SHORT_CHAR];
// 	char szFrqPatternType[LENGTH_OF_SHORT_CHAR];
// 	char szHasFrqPeriod[LENGTH_OF_SHORT_CHAR];
// 	int nFrqPeriodMean; // �ػ� 10khz
// 	int nFrqPeriodMax; // �ػ� 10khz
// 	int nFrqPeriodMin; // �ػ� 10khz
// 	int nFrqPositionCount;
// 	int nFrqElementCount;
// 	int nFrqMean; // �ػ� 10khz
// 	int nFrqMax;
// 	int nFrqMin;
// 	int nFrqDev;
// 	int arrFrqElement[64];
// 	char szFrqElement[LENGTH_OF_FRQ_OR_PRI_SEQ_CHAR];		
// 	char szPriType[LENGTH_OF_SHORT_CHAR];
// 	char szPriPatternType[LENGTH_OF_SHORT_CHAR];
// 	char szHasPriPeriod[LENGTH_OF_SHORT_CHAR];
// 	int nPriPeriodMean; // us
// 	int nPriPeriodMax; // us
// 	int nPriPeriodMin; // us
// 	int nPriPositionCount;
// 	int nPriElementCount;
// 	int nPriMean;
// 	int nPriMax;
// 	int nPriMin;
// 	int nPriDev;
// 	int iPriJitterRatio;
// 	int arrPriElement[64];
// 	char szPriElement[LENGTH_OF_FRQ_OR_PRI_SEQ_CHAR];
// 	int nPwMean;
// 	int nPwMax;
// 	int nPwMin;
// 	int nPwDev;
// 	int nPaMean;
// 	int nPaMax;
// 	int nPaMin;
// 	int nPaDev;
// 	char szScanType[LENGTH_OF_SHORT_CHAR];
// 	int nScanPeriodMinMicroSec;
// 	int nScanPeriodMaxMicroSec;
// 	float fScanPeriodMinHz;
// 	float fScanPeriodMaxHz;
// 	char szHasIntraModulation[LENGTH_OF_SHORT_CHAR];
// 	int nIntraFrqChangeWidthMin;
// 	int nIntraFrqChangeWidthMax;
// 	int arrPriPerGroup[MAX_PRI_PER_GROUP];
// 	int arrPaDiffPerGroup[MAX_PADIFF_PER_GROUP];
// 	// ��ġ���� ����
// 	int nPosEstLat;
// 	int nPosEstLong;
// 	int nManualPosEstLat;
// 	int nManualPosEstLong;
// 	int nRepresentPosEstLat;
// 	int nRepresentPosEstLong;
// 	int nAltitude; // meter
// 	int nCEP;// meter
// 	int nLengthOfMajorAxis; // meter
// 	int nLengthOfMinorAxis;
// 	int nEEPTiltAngle; // 0.1�� �ػ�
// 	int nManualPosEstPreferred; // ������ġ������ �켱�Ѵٴ� indicator.
// 	unsigned int uiNewAetId;
// 	E_EMITTER_OPCODE eOpCodeForUpdate; // ������Ʈ �Ǵ� ������ �پ��� ����, Update Key�� �ּ� �ĺ��ϱ�� ��. ������. 2018.4.10.
// 
// 	char szFirstReportTime[LENGTH_OF_SHORT_CHAR];
// 	char szFinalReportTime[LENGTH_OF_SHORT_CHAR];
// 	char szFinalAlarmTime[LENGTH_OF_SHORT_CHAR];
// 	char szPolarization[LENGTH_OF_SHORT_CHAR];
// 	char szStat[LENGTH_OF_SHORT_CHAR];
// 
// 	SAbtInfoToWnd()
// 		:ullSeqNum(0),
// 	uiAbtId(0),
// 	uiAetId(0),
// 	iRadarModeIndex(0), //LDRA�߰�
// 	iThreatIndex(0),//LDRA�߰�
// 	nRadarModePriority(0),
// 	fDistanceErrFromThreat(0.0),
// 	nRadarPriority(0),	
// 	nNumOfPPG(0),
// 	nNumOfLOB(0),
// 	iDoaMean(0), // �ػ� 0.1 deg
// 	iDoaMax(0),
// 	iDoaMin(0),
// 	iDoaDev(0),
// 	iDoaStd(0),
// 	nFrqPeriodMean(0), // �ػ� 10khz
// 	nFrqPeriodMax(0), // �ػ� 10khz
// 	nFrqPeriodMin(0), // �ػ� 10khz
// 	nFrqPositionCount(0),
// 	nFrqElementCount(0),
// 	nFrqMean(0), // �ػ� 10khz
// 	nFrqMax(0),
// 	nFrqMin(0),
// 	nFrqDev(0),	
// 	nPriPeriodMean(0), // us
// 	nPriPeriodMax(0), // us
// 	nPriPeriodMin(0), // us
// 	nPriPositionCount(0),
// 	nPriElementCount(0),
// 	nPriMean(0),
// 	nPriMax(0),
// 	nPriMin(0),
// 	nPriDev(0),
// 	iPriJitterRatio(0),	
// 	nPwMean(0),
// 	nPwMax(0),
// 	nPwMin(0),
// 	nPwDev(0),
// 	nPaMean(0),
// 	nPaMax(0),
// 	nPaMin(0),
// 	nPaDev(0),
// 	nScanPeriodMinMicroSec(0),
// 	nScanPeriodMaxMicroSec(0),
// 	fScanPeriodMinHz(0),
// 	fScanPeriodMaxHz(0),
// 	nIntraFrqChangeWidthMin(0),
// 	nIntraFrqChangeWidthMax(0),
// 	nPosEstLat(0),
// 	nPosEstLong(0),
// 	nManualPosEstLat(0),
// 	nManualPosEstLong(0),
// 	nRepresentPosEstLat(0),
// 	nRepresentPosEstLong(0),
// 	nAltitude(0), // meter
// 	nCEP(0),// meter
// 	nLengthOfMajorAxis(0), // meter
// 	nLengthOfMinorAxis(0),
// 	nEEPTiltAngle(0),  // 0.1�� �ػ�
// 	nManualPosEstPreferred(0),
// 	uiNewAetId(0),
// 	eOpCodeForUpdate(E_EO_NOT_AVAILABLE)
// 	{
// 		//@start_�̽ÿ�
// 		memset(szMissionId, 0,LENGTH_OF_MISSION_ID+1);
// 		memset(szMissionName, 0,SIZE_OF_TASK_NAME);
// 		memset(szTaskType, 0, LENGTH_OF_SHORT_CHAR);
// 		memset(szElnotPri, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szModeCodePri, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szModulationCode, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szNickName, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szPriFuncCode, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szPlatform, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szSigType, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szFirstDetectTimeInfo, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szFinalDetectTimeInfo, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szValidity, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szFrqType, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szFrqPatternType, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szHasFrqPeriod, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(arrFrqElement, 0, 64);
// 		memset(szFrqElement,0,LENGTH_OF_FRQ_OR_PRI_SEQ_CHAR);
// 		memset(szPriType, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szPriPatternType, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szHasPriPeriod, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(arrPriElement,0, 64);
// 		memset(szPriElement, 0, LENGTH_OF_FRQ_OR_PRI_SEQ_CHAR);
// 		memset(szScanType, 0,LENGTH_OF_SHORT_CHAR);
// 		memset(szHasIntraModulation, 0,LENGTH_OF_SHORT_CHAR);		
// 		memset(arrPriPerGroup, 0, sizeof(int)*MAX_PRI_PER_GROUP);
// 		memset(arrPaDiffPerGroup, 0, sizeof(int)*MAX_PADIFF_PER_GROUP);
// 		memset(szFirstReportTime,0, LENGTH_OF_SHORT_CHAR);
// 		memset(szFinalReportTime,0, LENGTH_OF_SHORT_CHAR);
// 		memset(szFinalAlarmTime,0, LENGTH_OF_SHORT_CHAR);
// 		memset(szPolarization, 0, LENGTH_OF_SHORT_CHAR);		
// 		memset(szStat, 0, LENGTH_OF_SHORT_CHAR);
// 		memset(szPlaceNameKor, 0, LENGTH_OF_SHORT_CHAR);
// 		//@end_�̽ÿ�
// 	}
// };
// 
// struct SAetInfoToWnd
// {
// 	unsigned long long ullSeqNum;
// 	char szMissionId[LENGTH_OF_MISSION_ID+1];
// 	char szMissionName[SIZE_OF_TASK_NAME];
// 	char szTaskType[LENGTH_OF_SHORT_CHAR];
// 	unsigned int uiAetId;
// 
// 	//CED ���� ����
// 	int iRadarIndex;
// 	int iThreatIndex;
// 	char szElnotPri[LENGTH_OF_SHORT_CHAR];
// 	char szIdResult[LENGTH_OF_SHORT_CHAR];
// 	char szNickName[LENGTH_OF_SHORT_CHAR];
// 	char szPriFuncCode[LENGTH_OF_SHORT_CHAR];
// 	int nRadarPriority;
// 	
// 	//EOB ����
// 	unsigned int uiPinNum;
// 	char szThreatName[LENGTH_OF_SHORT_CHAR];
// 	char szBENumber[LENGTH_OF_SHORT_CHAR];
// 	char szThreatFuncCode[LENGTH_OF_SHORT_CHAR];
// 	int nThreatPriority;
// 	float fDistanceErrFromThreat; // Nautical Mile
// 	int iEquipNumber;
// 	//char szEquipNumber[LENGTH_OF_SHORT_CHAR];
// 
// 	int nNumOfLOB;
// 	int nNumOfBeam;
// 	char szFirstDetectTimeInfo[LENGTH_OF_SHORT_CHAR];
// 	char szFinalDetectTimeInfo[LENGTH_OF_SHORT_CHAR];
// 	char szValidity[LENGTH_OF_SHORT_CHAR];
// 	int nFrqMean; // �ػ� 10khz
// 	int nFrqMax;
// 	int nFrqMin;
// 	int nFrqDev;
// 	int nPriMean;
// 	int nPriMax;
// 	int nPriMin;
// 	int nPriDev;
// 	int nPrfMean;
// 	int nPrfMax;
// 	int nPrfMin;
// 	int nPrfDev;
// 	int nPrfMinPPS;
// 	int nPrfMaxPPS;
// 	int nPwMean;
// 	int nPwMax;
// 	int nPwMin;
// 	int nPwDev;
// 	int nPaMean;
// 	int nPaMax;
// 	int nPaMin;
// 	int nPaDev;
// 	// ��ġ���� ����
// 	char szPosEstValidity[LENGTH_OF_SHORT_CHAR];
// 	int nPosEstLat; // �ڵ����� ���ŵǴ� ��ġ ���Ⱚ
// 	int nPosEstLong; // �ڵ����� ���ŵǴ� ��ġ ���Ⱚ
// 	int nManualPosEstLat; // ������ ������ġ ���Ⱚ
// 	int nManualPosEstLong; // ������ ������ġ ���Ⱚ	
// 	int nRepresentPosEstLat; // ��ǥ ��ġ ����
// 	int nRepresentPosEstLong; // ��ǥ ��ġ �浵
// 	int nAltitude; // meter
// 	int nCEP;// meter
// 	int nLengthOfMajorAxis; // meter
// 	int nLengthOfMinorAxis;
// 	int nEEPTiltAngle; // 0.1�� �ػ�
// 	int nManualPosEstPreferred; // ������ġ������ �켱�Ѵٴ� indicator.
// 
// 	E_EMITTER_OPCODE eOpCodeForUpdate; // ������Ʈ �Ǵ� ������ �پ��� ����, Update Key�� �ּ� �ĺ��ϱ�� ��. ������. 2018.4.10.
// 	char szManualPosEstPreferred[LENGTH_OF_SHORT_CHAR];
// 	char szFirstReportTime[LENGTH_OF_SHORT_CHAR];
// 	char szFinalReportTime[LENGTH_OF_SHORT_CHAR];
// 	char szFinalAlarmTime[LENGTH_OF_SHORT_CHAR];
// 
// 	char szStat[LENGTH_OF_SHORT_CHAR];
// 
// 	SAetInfoToWnd()
// 	:ullSeqNum(0)
// 	 ,uiAetId(0)
// 	 ,iRadarIndex(0)
// 	 ,iThreatIndex(0)
// 	 ,nRadarPriority(0)
// 	 ,uiPinNum(0)
// 	 ,nThreatPriority(0)
// 	 ,fDistanceErrFromThreat(0.0) // Nautical Mile
// 	 ,iEquipNumber(0)
// 	 ,nNumOfLOB(0)
// 	 ,nNumOfBeam(0)
// 	 ,nFrqMean(0) // �ػ� 10khz
// 	 ,nFrqMax(0)
// 	 ,nFrqMin(0)
// 	 ,nFrqDev(0)
// 	 ,nPriMean(0)
// 	 ,nPriMax(0)
// 	 ,nPriMin(0)
// 	 ,nPriDev(0)
// 	 ,nPrfMean(0)
// 	 ,nPrfMax(0)
// 	 ,nPrfMin(0)
// 	 ,nPrfDev(0)
// 	 ,nPrfMinPPS(0)
// 	 ,nPrfMaxPPS(0)
// 	 ,nPwMean(0)
// 	 ,nPwMax(0)
// 	 ,nPwMin(0)
// 	 ,nPwDev(0)
// 	 ,nPaMean(0)
// 	 ,nPaMax(0)
// 	 ,nPaMin(0)
// 	 ,nPaDev(0)
// 	 ,nPosEstLat(0)
// 	 ,nPosEstLong(0)
// 	 ,nManualPosEstLat(0)
// 	 ,nManualPosEstLong(0)
// 	 ,nRepresentPosEstLat(0)
// 	 ,nRepresentPosEstLong(0)
// 	 ,nAltitude(0) // meter
// 	 ,nCEP(0)// meter
// 	 ,nLengthOfMajorAxis(0) // meter
// 	 ,nLengthOfMinorAxis(0)
// 	 ,nEEPTiltAngle(0) // 0.1�� �ػ�
// 	 ,nManualPosEstPreferred(0)
// 	 ,eOpCodeForUpdate(E_EO_NOT_AVAILABLE)
// 	 //,szManualPosEstPreferred(0)
// 	{
// 		//@start_�̽ÿ�
// 		memset(szMissionId, 0, LENGTH_OF_MISSION_ID+1);
// 		memset(szMissionName, 0, SIZE_OF_TASK_NAME);
// 		memset(szTaskType, 0, SIZE_OF_TASK_NAME);
// 		memset(szElnotPri, 0, LENGTH_OF_SHORT_CHAR);
// 		memset(szIdResult, 0, LENGTH_OF_SHORT_CHAR);
// 		memset(szNickName, 0, LENGTH_OF_SHORT_CHAR);
// 		memset(szPriFuncCode, 0, LENGTH_OF_SHORT_CHAR);		
// 		memset(szThreatName, 0, LENGTH_OF_SHORT_CHAR);
// 		memset(szBENumber, 0, LENGTH_OF_SHORT_CHAR);
// 		memset(szThreatFuncCode, 0, LENGTH_OF_SHORT_CHAR);	
// 		//memset(szEquipNumber, 0, LENGTH_OF_SHORT_CHAR);	
// 		memset(szFirstDetectTimeInfo, 0, LENGTH_OF_SHORT_CHAR);
// 		memset(szFinalDetectTimeInfo, 0, LENGTH_OF_SHORT_CHAR);
// 		memset(szValidity, 0, LENGTH_OF_SHORT_CHAR);	
// 		memset(szPosEstValidity, 0, LENGTH_OF_SHORT_CHAR);
// 		memset(szManualPosEstPreferred, 0, LENGTH_OF_SHORT_CHAR);
// 		memset(szFirstReportTime, 0, LENGTH_OF_SHORT_CHAR);
// 		memset(szFinalReportTime, 0, LENGTH_OF_SHORT_CHAR);
// 		memset(szFinalAlarmTime, 0, LENGTH_OF_SHORT_CHAR);
// 		memset(szStat, 0, LENGTH_OF_SHORT_CHAR);
// 
// 		//@end_�̽ÿ�
// 		
// 	}
// };

// GUI->LOGIC
// struct SCommonListWndSelect
// {
// 	int nAetId;
//     int nAbtId;
// 	int nLobId;
// 	std::string strRecentTime;
// 	SCommonListWndSelect()
// 	:nAetId(0)
// 	,nAbtId(0)
// 	,nLobId(0)
// 	,strRecentTime("")
//     {}
// };


#define NUM_OF_MANUAL_INPUT_COUNT	_	36

#define DF_MANUAL_SIG_SIGTYPE				(0x0000000000000001)		/*��ȣ����-----*/
#define DF_MANUAL_SIG_MODTYPE				(0x0000000000000002)		/*��������-----*/
#define DF_MANUAL_SIG_POLARIZATION		(0x0000000000000004)		/*�ؼ�*/
#define DF_MANUAL_DOA_MEAN					(0x0000000000000008)		/*������ ���-----*/
#define DF_MANUAL_DOA_MAX					(0x0000000000000010)		/*������ �ִ�-----*/
#define DF_MANUAL_DOA_MIN						(0x0000000000000020)		/*������ �ּ�-----*/
#define DF_MANUAL_POSEST_NS					(0x0000000000000040)		/*NS*/
#define DF_MANUAL_POSEST_LAT				(0x0000000000000080)		/*����*/
#define DF_MANUAL_POSEST_EW					(0x0000000000000100)		/*EW*/
#define DF_MANUAL_POSEST_LONG				(0x0000000000000200)		/*�浵*/
#define DF_MANUAL_FRQ_FRQTYPE				(0x0000000000000400)		/*���ļ� ����-----*/
#define DF_MANUAL_FRQ_MEAN					(0x0000000000000800)		/*���ļ� ���-----*/
#define DF_MANUAL_FRQ_MAX						(0x0000000000001000)		/*���ļ� �ִ�-----*/
#define DF_MANUAL_FRQ_MIN						(0x0000000000002000)		/*���ļ� �ּ�-----*/
#define DF_MANUAL_FRQ_CHNG_WIDTH		(0x0000000000004000)		/*���ļ� ������*/
#define DF_MANUAL_FRQ_PERIOD				(0x0000000000008000)		/*�ֱ�-----*/
#define DF_MANUAL_FRQ_STEP					(0x0000000000010000)		/*�ܼ�*/
#define DF_MANUAL_FRQ_LEVEL					(0x0000000000020000)		/*����*/
#define DF_MANUAL_PRI_PRITYPE					(0x0000000000040000)		/*PRI ����-----*/
#define DF_MANUAL_PRI_MEAN					(0x0000000000080000)		/*PRI ���-----*/
#define DF_MANUAL_PRI_MAX						(0x0000000000100000)		/*PRI �ִ�-----*/
#define DF_MANUAL_PRI_MIN						(0x0000000000200000)		/*PRI �ּ�-----*/
#define DF_MANUAL_PRI_CHNG_WIDTH		(0x0000000000400000)		/*PRI ������*/
#define DF_MANUAL_PRI_PERIOD					(0x0000000000800000)		/*�ֱ�-----*/
#define DF_MANUAL_PRI_STEP						(0x0000000001000000)		/*�ܼ�*/
#define DF_MANUAL_PRI_CHNG_RATIO			(0x0000000002000000)		/*PRI �����*/
#define DF_MANUAL_PRI_LEVEL					(0x0000000004000000)		/*����*/
#define DF_MANUAL_ELNOT							(0x0000000008000000)		/*ELNOT-----*/
#define DF_MANUAL_PA_MEAN						(0x0000000010000000)		/*PA ���-----*/
#define DF_MANUAL_PA_MAX						(0x0000000020000000)		/*PA �ִ�-----*/
#define DF_MANUAL_PA_MIN						(0x0000000040000000)		/*PA �ּ�-----*/
#define DF_MANUAL_PA_SCAN_TYPE				(0x0000000080000000)		/*��ĵŸ��-----*/
#define DF_MANUAL_PA_SCAN_PERIOD			(0x0000000100000000)		/*��ĵ�ֱ�*/
#define DF_MANUAL_PW_MEAN					(0x0000000200000000)		/*PW ���*/
#define DF_MANUAL_PW_MAX					(0x0000000400000000)		/*PW �ִ�*/
#define DF_MANUAL_PW_MIN						(0x0000000800000000)		/*PW �ּ�*/
#define DF_MANUAL_PIN_NUM					(0x000000100000000)			/*PIN NUMBER _ 2015.10.31. �߰�*/

struct SFLTMData
{
	int nFrqType;
	int nPriType;
	int nMeanFrq;
	int nMeanPri;
	int nMeanPw;
	int nMeanPa;
	int nMeanDoa;
	int nAntDir;
	SFLTMData()
		:nFrqType(0)
		,nPriType(0)
		,nMeanFrq(0)
		,nMeanPri(0)
		,nMeanPw(0)
		,nMeanPa(0)
		,nMeanDoa(0)
		,nAntDir(0)
	{};
};


//////////////////////////////////////////////////////////////////////////
// ���ü/��/LOB â���� CED/EOB â���� �����͸� �����ϱ� ���� ����ü ����
// typedef struct {
// 	SRadar stSRadar;
// 	EnumLibType enLibType;
// 
// } SELCALLCEDLIB;

// CED ������ �����Ǵ� ����
// typedef struct {
// 	int nAETId;
// 	int nABTId;
// 
// 	char szMissionId[LENGTH_OF_MISSION_ID+2];
// 
// } SELCEDCREATE ;


typedef struct {
	bool bCheckBox;
	float low;
	float high;
} SELCHECKRANGE;

typedef enum {
	enumIncludeAllInput=0,						// �Է¹����� �������
	enumOverlapWithInnerInput,				// �Է� ���� ���ο� ��ø
	enumOverlapWithInput							// �Է� ������ ��ø

} EnumRangeSearchRef;


// typedef struct {
// 	//EnumLibType enLibType;	// �⺻�� �Ǵ� �ǹ���
// 
// 	CString strELNOT;
// 	CString strNickName;
// 	int iPin;
// 	EnumFunctionCodes eFunctionCodes_ForGUI;	// �ֱ�� �ڵ�
// 	EnumRadarStatus eStatus;									// ����
// 	CString strAssocWeapSys;									// ���� ���� ü��
// 	CString strAssocPlatform;									// ���� �÷���
// 
// 	EnumRangeSearchRef eRangeSearchRef;			// ���� �˻� ����
// 
// 	SELCHECKRANGE stFreq;											// ���ļ� ����
// 	SELCHECKRANGE stPRI;											// PRI ����
// } SELRADARLIST_SEARCH_FILTER;


typedef struct {
	int low;
	int high;
} SELRANGE;

// typedef struct stSELRADARMODELIST_SEARCH_FILTER{
// 	EnumLibType enLibType;	// �⺻�� �Ǵ� �ǹ���
// 
// 	// ���ü
// 	CString strELNOT;
// 	CString strNickName;
// 	int iPin;
// 	EnumRadarStatus eStatus;									// ����
// 	CString strAssocWeapSys;									// ���� ���� ü��
// 	CString strAssocPlatform;									// ���� �÷���
// 
// 	// �������
// 	PlatformCode::EnumPlatformCode ePlatform;
// 	EnumValidationCode eValidation;
// 	EnumFunctionCodes	eFunctionCode;	// �ֱ�� �ڵ�
// 	SignalType::EnumSignalType eSignalType;									//��ȣ���� (Pulsed, CW, EA) enum����
// 	ScanType::EnumScanType	eScanPrimaryType;								//�� ��ĵŸ�� �ڵ�(SCAN_TYPE_CODE ����)
// 	CString strModulationCode;												//* �����ڵ�(2) [����ʿ�]
// 	PolizationCode::EnumPolizationCode ePolarization;
// 	EnumRangeSearchRef eRangeSearchRef;
// 	float fRF_TypicalMin;													//�ְ�������: ���ļ� Typical �ּ�(MHz)
// 	float fRF_TypicalMax;													//�ְ�������: ���ļ� Typical �ִ�(MHz)
// 	float fPRI_TypicalMin;													//PRI TYPICAL (USEC) �ּ�
// 	float fPRI_TypicalMax;													//PRI TYPICAL (USEC) �ִ�
// 	float fScanPrimaryTypicalMin;											//�� ��ĵ �ֱⰪ�� TYPICAL (SEC) �ּ�
// 	float fScanPrimaryTypicalMax;											//�� ��ĵ �ֱⰪ�� TYPICAL (SEC) �ִ�
// 	float fPD_TypicalMin;													//PD TYPICAL �� (USEC) �ּ�
// 	float fPD_TypicalMax;													//PD TYPICAL �� (USEC) �ִ�
// 
// 	// ���ļ�
// 	ContinuityCode::EnumContinuityCode eRF_Continuity;						//���ļ� ���������Ӽ�: RF ��ȭ�� ���Ӽ� (CONTINUITY_CODE ����)
// 	PatternCode::EnumPatternCode eRF_Pattern;								//RF ��ȭ�� ���� ���� (PATTERN_CODE ����)
// 	EnumRF_LagacyTypeCode eRF_LagacyType;									//RF_LEGACY_TYPE_CODE ����
// 	SELRANGE stRF_NumElements;													//NUMBER	PRI ELEMENT ��
// 	SELRANGE stRF_NumPositions;													//NUMBER	PRI POSITION ��
// 
// 	// PD
// 	ContinuityCode::EnumContinuityCode ePD_Continuity;						//PD���� ���Ӽ� (CONTINUITY_CODE ����)
// 	PatternCode::EnumPatternCode ePD_Pattern;								//PD���� ���Ͽ��� (PATTERN_CODE ����)
// 	SELRANGE stPD_NumElements;													//PD ELEMENT ��
// 	SELRANGE stPD_NumPositions;													//PD Positon ��
// 	float fPD_PatternPeriodMin;												//PD �����ֱ� (USEC) �ּ�
// 	float fPD_PatternPeriodMax;												//PD �����ֱ� (USEC) �ִ�
// 
// 	// PRI
// 	ContinuityCode::EnumContinuityCode ePRI_Continuity;						//PRI ��ȭ�� ���Ӽ� (CONTINUITY_CODE ����)
// 	PatternCode::EnumPatternCode ePRI_Pattern;								//NUMBER PRI��ȭ�� ���Ͽ���	(PATTERN_CODE ����)
// 	EnumPRI_LegacyTypeCode ePRI_LagacyType;									//PRI_LAGACY_TYPE_CODE ����
// 	SELRANGE stPRI_NumElements;													//NUMBER	PRI ELEMENT ��
// 	SELRANGE stPRI_NumPositions;													//NUMBER	PRI POSITION ��
// 	float fPRI_FramePeriodMin;												//PRI FRAME �ֱ� (USEC) �ּ�
// 	float fPRI_FramePeriodMax;												//PRI FRAME �ֱ� (USEC) �ִ�
// 	float fPRI_SubframePeriodMin;											//PRI SUBFRAME �ֱ� (USEC) �ּ�
// 	float fPRI_SubframePeriodMax;											//PRI SUBFRAME �ֱ� (USEC) �ִ�
// 	float fPRI_PPG_Min;														//Pulse Per Group �ּ�
// 	float fPRI_PPG_Max;														//Pulse Per Group �ִ�
// 
// 	// �޽��� ����
// 	MOP_CW_ModulationType::EnumMOP_CW_ModulationType eMOP_CW_ModulationType;	//�޽� �� �Ǵ� CW �������� (MOP_CW_MOD_TYPE_CODE ����)
// 	EnumMOP_CW_LegacyType eMOP_CW_LegacyType;					//���������� LEGACY TERM (MOP_CW_LEGACY_TYPE_CODE)
// 
// 	stSELRADARMODELIST_SEARCH_FILTER()
// 	{
// 		enLibType = E_EL_LIB_TYPE_NORMAL;	// �⺻�� �Ǵ� �ǹ���
// 
// 		// ���ü
// 		strELNOT="";
// 		strNickName="";
// 		iPin = 0;
// 		eStatus = enumUndefinedRadarStatus;									// ����
// 		strAssocWeapSys="";									// ���� ���� ü��
// 		strAssocPlatform="";									// ���� �÷���
// 
// 		// �������
// 		ePlatform = PlatformCode::enumUndefinedPlatformCode;
// 		eValidation=enumUndefinedValidationCode;
// 		eFunctionCode=enumUndefinedFunctionCode;	// �ֱ�� �ڵ�
// 		eSignalType=SignalType::enumSignalUndefined;									//��ȣ���� (Pulsed, CW, EA) enum����
// 		eScanPrimaryType=ScanType::enumUndefinedScanType;								//�� ��ĵŸ�� �ڵ�(SCAN_TYPE_CODE ����)
// 		strModulationCode="";												//* �����ڵ�(2) [����ʿ�]
// 		ePolarization=PolizationCode::enumUndefinedPolization;
// 		eRangeSearchRef=enumIncludeAllInput;
// 		fRF_TypicalMin=0.0f;													//�ְ�������: ���ļ� Typical �ּ�(MHz)
// 		fRF_TypicalMax=0.0f;													//�ְ�������: ���ļ� Typical �ִ�(MHz)
// 		fPRI_TypicalMin=0.0f;													//PRI TYPICAL (USEC) �ּ�
// 		fPRI_TypicalMax=0.0f;													//PRI TYPICAL (USEC) �ִ�
// 		fScanPrimaryTypicalMin=0.0f;											//�� ��ĵ �ֱⰪ�� TYPICAL (SEC) �ּ�
// 		fScanPrimaryTypicalMax=0.0f;											//�� ��ĵ �ֱⰪ�� TYPICAL (SEC) �ִ�
// 		fPD_TypicalMin=0.0f;													//PD TYPICAL �� (USEC) �ּ�
// 		fPD_TypicalMax=0.0f;													//PD TYPICAL �� (USEC) �ִ�
// 
// 		// ���ļ�
// 		eRF_Continuity=ContinuityCode::enumUndefinedContinuityCode;						//���ļ� ���������Ӽ�: RF ��ȭ�� ���Ӽ� (CONTINUITY_CODE ����)
// 		eRF_Pattern=PatternCode::enumUndefinedPatternCode;								//RF ��ȭ�� ���� ���� (PATTERN_CODE ����)
// 		eRF_LagacyType=enumUndefinedRF_LagacyType;									//RF_LEGACY_TYPE_CODE ����
// 		stRF_NumElements = SELRANGE();													//NUMBER	PRI ELEMENT ��
// 		stRF_NumPositions = SELRANGE();													//NUMBER	PRI POSITION ��
// 
// 		// PD
// 		ePD_Continuity=ContinuityCode::enumUndefinedContinuityCode;						//PD���� ���Ӽ� (CONTINUITY_CODE ����)
// 		ePD_Pattern=PatternCode::enumUndefinedPatternCode;								//PD���� ���Ͽ��� (PATTERN_CODE ����)
// 		stPD_NumElements = SELRANGE();													//PD ELEMENT ��
// 		stPD_NumPositions = SELRANGE();													//PD Positon ��
// 		fPD_PatternPeriodMin=0.0f;												//PD �����ֱ� (USEC) �ּ�
// 		fPD_PatternPeriodMax=0.0f;												//PD �����ֱ� (USEC) �ִ�
// 
// 		// PRI
// 		ePRI_Continuity=ContinuityCode::enumUndefinedContinuityCode;						//PRI ��ȭ�� ���Ӽ� (CONTINUITY_CODE ����)
// 		ePRI_Pattern=PatternCode::enumUndefinedPatternCode;								//NUMBER PRI��ȭ�� ���Ͽ���	(PATTERN_CODE ����)
// 		ePRI_LagacyType=enumUndefinedPRI_LegacyType;									//PRI_LAGACY_TYPE_CODE ����
// 		stPRI_NumElements = SELRANGE();													//NUMBER	PRI ELEMENT ��
// 		stPRI_NumPositions = SELRANGE();													//NUMBER	PRI POSITION ��
// 		fPRI_FramePeriodMin=0.0f;												//PRI FRAME �ֱ� (USEC) �ּ�
// 		fPRI_FramePeriodMax=0.0f;												//PRI FRAME �ֱ� (USEC) �ִ�
// 		fPRI_SubframePeriodMin=0.0f;											//PRI SUBFRAME �ֱ� (USEC) �ּ�
// 		fPRI_SubframePeriodMax=0.0f;											//PRI SUBFRAME �ֱ� (USEC) �ִ�
// 		fPRI_PPG_Min=0.0f;														//Pulse Per Group �ּ�
// 		fPRI_PPG_Max=0.0f;														//Pulse Per Group �ִ�
// 
// 		// �޽��� ����
// 		eMOP_CW_ModulationType=MOP_CW_ModulationType::enumUndefinedMOP_CW_ModulationType;	//�޽� �� �Ǵ� CW �������� (MOP_CW_MOD_TYPE_CODE ����)
// 		eMOP_CW_LegacyType=enumUndefinedMOP_CW_LegacyType;					//���������� LEGACY TERM (MOP_CW_LEGACY_TYPE_CODE)
// 	}
// 
// } SELRADARMODELIST_SEARCH_FILTER;

// typedef struct {
// 	EnumLibType enLibType;	// �⺻�� �Ǵ� �ǹ���
// 
// 	int iPin;
// 	CString strBE_Number;
// 	CString strFacilityName;											//������Ī (������)
// 
// 	EnumFunctionCodes ePrimaryFunction_ForGUI;		//�ֱ���ڵ�
// 	FriendOrFOE::EnumFriendOrFOE eFriendOrFOE;		//���Ʊ��� enum��
// 
// 	CString strPlaceNameKor;										//�ѱ����� (50)
// 	CString strPlaceNameEng;										//�������� (50)
// 
// 	CString strCategory;												//��ü ���� ī�װ� (99999 ����)(5)
// 	CString strADA;														//���������� (ADA) (AA123 ����) (5)
// 
// 	CountryCode::EnumCountryCode eUserCountry;	//��뱹���� enum���� ǥ��
// 	PlatformCode::EnumPlatformCode ePlatform;			//ž�� �÷����� ���� (PLATFORM_CODE ����)
// 
// } SELTHREATLIST_SEARCH_FILTER;

#endif




/************************************************************************************
*   ELINT View -> Logic�� ���Ǵ� �ڷ��� ����ü
*************************************************************************************/

#endif