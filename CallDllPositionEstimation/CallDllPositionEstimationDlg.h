
// CallDllPositionEstimationDlg.h : ��� ����
//

#pragma once


// CCallDllPositionEstimationDlg ��ȭ ����
class CCallDllPositionEstimationDlg : public CDialogEx
{
// �����Դϴ�.
public:
	CCallDllPositionEstimationDlg(CWnd* pParent = NULL);	// ǥ�� �������Դϴ�.

// ��ȭ ���� �������Դϴ�.
	enum { IDD = IDD_CALLDLLPOSITIONESTIMATION_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV �����Դϴ�.

public:
	bool OpenFile( CString &strPathname, char *pTitle );


// �����Դϴ�.
protected:
	HICON m_hIcon;

	// ������ �޽��� �� �Լ�
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedButtonFile();
};
