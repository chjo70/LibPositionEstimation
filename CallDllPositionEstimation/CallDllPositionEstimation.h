
// CallDllPositionEstimation.h : PROJECT_NAME ���� ���α׷��� ���� �� ��� �����Դϴ�.
//

#pragma once

#ifndef __AFXWIN_H__
	#error "PCH�� ���� �� ������ �����ϱ� ���� 'stdafx.h'�� �����մϴ�."
#endif

#include "resource.h"		// �� ��ȣ�Դϴ�.


// CCallDllPositionEstimationApp:
// �� Ŭ������ ������ ���ؼ��� CallDllPositionEstimation.cpp�� �����Ͻʽÿ�.
//

class CCallDllPositionEstimationApp : public CWinApp
{
public:
	CCallDllPositionEstimationApp();

// �������Դϴ�.
public:
	virtual BOOL InitInstance();

// �����Դϴ�.

	DECLARE_MESSAGE_MAP()
};

extern CCallDllPositionEstimationApp theApp;