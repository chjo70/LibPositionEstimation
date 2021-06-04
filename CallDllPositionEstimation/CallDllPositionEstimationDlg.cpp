
// CallDllPositionEstimationDlg.cpp : ���� ����
//

#include "stdafx.h"
#include "CallDllPositionEstimation.h"
#include "CallDllPositionEstimationDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#include "../DllPositionEstimation/Dll/PositionEstimation.h"


// ���� ���α׷� ������ ���Ǵ� CAboutDlg ��ȭ �����Դϴ�.

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// ��ȭ ���� �������Դϴ�.
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV �����Դϴ�.

// �����Դϴ�.
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


// CCallDllPositionEstimationDlg ��ȭ ����




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


// CCallDllPositionEstimationDlg �޽��� ó����

BOOL CCallDllPositionEstimationDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// �ý��� �޴��� "����..." �޴� �׸��� �߰��մϴ�.

	// IDM_ABOUTBOX�� �ý��� ��� ������ �־�� �մϴ�.
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

	// �� ��ȭ ������ �������� �����մϴ�. ���� ���α׷��� �� â�� ��ȭ ���ڰ� �ƴ� ��쿡��
	//  �����ӿ�ũ�� �� �۾��� �ڵ����� �����մϴ�.
	SetIcon(m_hIcon, TRUE);			// ū �������� �����մϴ�.
	SetIcon(m_hIcon, FALSE);		// ���� �������� �����մϴ�.

	// TODO: ���⿡ �߰� �ʱ�ȭ �۾��� �߰��մϴ�.

	return TRUE;  // ��Ŀ���� ��Ʈ�ѿ� �������� ������ TRUE�� ��ȯ�մϴ�.
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

// ��ȭ ���ڿ� �ּ�ȭ ���߸� �߰��� ��� �������� �׸�����
//  �Ʒ� �ڵ尡 �ʿ��մϴ�. ����/�� ���� ����ϴ� MFC ���� ���α׷��� ��쿡��
//  �����ӿ�ũ���� �� �۾��� �ڵ����� �����մϴ�.

void CCallDllPositionEstimationDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // �׸��⸦ ���� ����̽� ���ؽ�Ʈ�Դϴ�.

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Ŭ���̾�Ʈ �簢������ �������� ����� ����ϴ�.
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// �������� �׸��ϴ�.
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

// ����ڰ� �ּ�ȭ�� â�� ���� ���ȿ� Ŀ���� ǥ�õǵ��� �ý��ۿ���
//  �� �Լ��� ȣ���մϴ�.
HCURSOR CCallDllPositionEstimationDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}



void CCallDllPositionEstimationDlg::OnBnClickedButtonFile()
{
	// TODO: ���⿡ ��Ʈ�� �˸� ó���� �ڵ带 �߰��մϴ�.
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
	if( true == OpenFile( strPathName, "��Ž �� ���� ��ǥ ���� �о����..." ) ) {
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

		strTitle.Format( "%s: %d ��" , strPathName, i );
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

	strTitle.Format( "%s: %d ��" , strPathName, i );
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
 * @author		��ö�� (churlhee.jo@lignex1.com)
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

    pWndFile = new CFileDialog(TRUE, NULL, NULL, OFN_ENABLESIZING | OFN_NONETWORKBUTTON | OFN_SHOWHELP | OFN_HIDEREADONLY, _T("���� ��� ���ϵ� (*.csv)|*.csv*|All Files (*.*)|*.*||") );

    _tcscpy_s( szinitDir, MAX_PATH, strFilepath.GetBuffer(0) );
    strFilepath.ReleaseBuffer();

    // Initializes m_ofn structure
    //pWndFile->m_ofn.lpstrTitle = ( LPCSTR ) pTitle;			// Ÿ��Ʋ��
    pWndFile->m_ofn.lpstrInitialDir = szinitDir;			// Ÿ��Ʋ��

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