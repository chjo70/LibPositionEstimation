
// CallDllPositionEstimationDlg.cpp : 구현 파일
//

#include "stdafx.h"
#include "CallDllPositionEstimation.h"
#include "CallDllPositionEstimationDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#include "../DllPositionEstimation/Dll/PositionEstimation.h"


// 응용 프로그램 정보에 사용되는 CAboutDlg 대화 상자입니다.

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// 대화 상자 데이터입니다.
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 지원입니다.

// 구현입니다.
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// CCallDllPositionEstimationDlg 대화 상자




CCallDllPositionEstimationDlg::CCallDllPositionEstimationDlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(CCallDllPositionEstimationDlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CCallDllPositionEstimationDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CCallDllPositionEstimationDlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_BUTTON_FILE, &CCallDllPositionEstimationDlg::OnBnClickedButtonFile)
END_MESSAGE_MAP()


// CCallDllPositionEstimationDlg 메시지 처리기

BOOL CCallDllPositionEstimationDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// 시스템 메뉴에 "정보..." 메뉴 항목을 추가합니다.

	// IDM_ABOUTBOX는 시스템 명령 범위에 있어야 합니다.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// 이 대화 상자의 아이콘을 설정합니다. 응용 프로그램의 주 창이 대화 상자가 아닐 경우에는
	//  프레임워크가 이 작업을 자동으로 수행합니다.
	SetIcon(m_hIcon, TRUE);			// 큰 아이콘을 설정합니다.
	SetIcon(m_hIcon, FALSE);		// 작은 아이콘을 설정합니다.

	// TODO: 여기에 추가 초기화 작업을 추가합니다.

	return TRUE;  // 포커스를 컨트롤에 설정하지 않으면 TRUE를 반환합니다.
}

void CCallDllPositionEstimationDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

// 대화 상자에 최소화 단추를 추가할 경우 아이콘을 그리려면
//  아래 코드가 필요합니다. 문서/뷰 모델을 사용하는 MFC 응용 프로그램의 경우에는
//  프레임워크에서 이 작업을 자동으로 수행합니다.

void CCallDllPositionEstimationDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // 그리기를 위한 디바이스 컨텍스트입니다.

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// 클라이언트 사각형에서 아이콘을 가운데에 맞춥니다.
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// 아이콘을 그립니다.
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

// 사용자가 최소화된 창을 끄는 동안에 커서가 표시되도록 시스템에서
//  이 함수를 호출합니다.
HCURSOR CCallDllPositionEstimationDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}



void CCallDllPositionEstimationDlg::OnBnClickedButtonFile()
{
	// TODO: 여기에 컨트롤 알림 처리기 코드를 추가합니다.
	int i=0;
	FILE *fi;

	SELPE_RESULT stSELPE_RESULT;

	STR_LOBS strLobs[1000];

	CString strPathName, strTitle, strValue;

	GetDlgItem(IDC_STATIC_LATITUDE)->SetWindowText( "-" );
	GetDlgItem(IDC_STATIC_LONGITUDE)->SetWindowText( "-" );
	GetDlgItem(IDC_STATIC_ALTITUDE)->SetWindowText( "-" );

	GetDlgItem(IDC_STATIC_CEP)->SetWindowText( "-" );
	GetDlgItem(IDC_STATIC_MAJOR_EEP)->SetWindowText( "-" );
	GetDlgItem(IDC_STATIC_MINOR_EEP)->SetWindowText( "-" );
	GetDlgItem(IDC_STATIC_THETA)->SetWindowText( "-" );

#if 1
	if( true == OpenFile( strPathName, "방탐 및 센서 좌표 파일 읽어오기..." ) ) {
		fi = fopen( (LPCSTR)(LPCTSTR) strPathName, "rt" );
		if( fi != NULL ) {
			memset( & strLobs[0], 0, sizeof(strLobs) );
			while( ! feof( fi ) ) {
				fscanf( fi, "%f,%f,%f,%f" , & strLobs[i].fDoa, & strLobs[i].fLatitude, & strLobs[i].fLongitude, & strLobs[i].fAltitude );

				if( i >= 1000 || ( strLobs[i].fLatitude == 0. && strLobs[i].fLongitude ==0. ) ) {
					break;
				}
				++i;
			}
			fclose( fi );

		}

		strTitle.Format( "%s: %d 개" , strPathName, i );
		SetWindowText( strTitle );

	}

#else
	strPathName = "D:/LibPositionEstimation/Debug/data.csv";

	SetWindowText( strPathName );

	fi = fopen( (LPCSTR)(LPCTSTR) strPathName, "rt" );
	if( fi != NULL ) {
		memset( & strLobs[0], 0, sizeof(strLobs) );
		while( ! feof( fi ) ) {
			fscanf( fi, "%f,%f,%f,%f" , & strLobs[i].fDoa, & strLobs[i].fLatitude, & strLobs[i].fLongitude, & strLobs[i].fAltitude );

			if( i >= 1000 || strLobs[i].fLatitude == 0. && strLobs[i].fLongitude ==0. ) {
				break;
			}
			++i;
		}
		fclose( fi );

	}

	strTitle.Format( "%s: %d 개" , strPathName, i );
	SetWindowText( strTitle );

#endif

	PE::CPositionEstimation::RunPositionEstimation( & stSELPE_RESULT, i, & strLobs[0] );

	if( stSELPE_RESULT.bResult == true ) {
		strValue.Format( "%.6f", stSELPE_RESULT.dLatitude );
		GetDlgItem(IDC_STATIC_LATITUDE)->SetWindowText( strValue );

		strValue.Format( "%.6f", stSELPE_RESULT.dLongitude );
		GetDlgItem(IDC_STATIC_LONGITUDE)->SetWindowText( strValue );

		strValue.Format( "%.1f", stSELPE_RESULT.dAltitude );
		GetDlgItem(IDC_STATIC_ALTITUDE)->SetWindowText( strValue );

		strValue.Format( "%.2f", stSELPE_RESULT.dCEP_error/1. );
		GetDlgItem(IDC_STATIC_CEP)->SetWindowText( strValue );

		strValue.Format( "%.2f", stSELPE_RESULT.dEEP_major_axis/1. );
		GetDlgItem(IDC_STATIC_MAJOR_EEP)->SetWindowText( strValue );
		strValue.Format( "%.2f", stSELPE_RESULT.dEEP_minor_axis/1. );
		GetDlgItem(IDC_STATIC_MINOR_EEP)->SetWindowText( strValue );

		strValue.Format( "%.2f", stSELPE_RESULT.dEEP_theta );
		GetDlgItem(IDC_STATIC_THETA)->SetWindowText( strValue );
	}
}

/**
 * @brief		OpenFile
 * @param		CString & strPathname
 * @param		char * pTitle
 * @return		bool
 * @author		조철희 (churlhee.jo@lignex1.com)
 * @version		0.0.1
 * @date		2021/02/19 18:28:25
 * @warning		
 */
bool CCallDllPositionEstimationDlg::OpenFile( CString &strPathname, char *pTitle )
{
    bool bRet = true;
    CFileDialog *pWndFile;
    TCHAR szinitDir[MAX_PATH];

    CString strFilepath;

    strFilepath = GetFilePath();

    pWndFile = new CFileDialog(TRUE, NULL, NULL, OFN_ENABLESIZING | OFN_NONETWORKBUTTON | OFN_SHOWHELP | OFN_HIDEREADONLY, _T("수집 목록 파일들 (*.csv)|*.csv*|All Files (*.*)|*.*||") );

    _tcscpy_s( szinitDir, MAX_PATH, strFilepath.GetBuffer(0) );
    strFilepath.ReleaseBuffer();

    // Initializes m_ofn structure
    //pWndFile->m_ofn.lpstrTitle = ( LPCSTR ) pTitle;			// 타이틀명
    pWndFile->m_ofn.lpstrInitialDir = szinitDir;			// 타이틀명

    // Call DoModal
    if (pWndFile->DoModal() == IDOK) {
        //m_strWindowtitle = pWndFile->GetFileTitle() + '.' + pWndFile->GetFileExt();
        strPathname = pWndFile->GetPathName();

    }
    else {
        bRet = false;
    }

    delete pWndFile;

    return bRet;

}