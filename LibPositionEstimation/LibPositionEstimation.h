
// LibPositionEstimation.h : PROJECT_NAME ���� ���α׷��� ���� �� ��� �����Դϴ�.
//

#pragma once

#ifndef __AFXWIN_H__
	#error "PCH�� ���� �� ������ �����ϱ� ���� 'stdafx.h'�� �����մϴ�."
#endif

#include "resource.h"		// �� ��ȣ�Դϴ�.


// CLibPositionEstimationApp:
// �� Ŭ������ ������ ���ؼ��� LibPositionEstimation.cpp�� �����Ͻʽÿ�.
//

class CLibPositionEstimationApp : public CWinApp
{
public:
	CLibPositionEstimationApp();

// �������Դϴ�.
public:
	virtual BOOL InitInstance();

// �����Դϴ�.

	DECLARE_MESSAGE_MAP()
};

extern CLibPositionEstimationApp theApp;