// DllPositionEstimation.h : DllPositionEstimation DLL�� �⺻ ��� �����Դϴ�.
//

#pragma once

#ifndef __AFXWIN_H__
	#error "PCH�� ���� �� ������ �����ϱ� ���� 'stdafx.h'�� �����մϴ�."
#endif

#include "resource.h"		// �� ��ȣ�Դϴ�.


// CDllPositionEstimationApp
// �� Ŭ������ ������ ������ DllPositionEstimation.cpp�� �����Ͻʽÿ�.
//

class CDllPositionEstimationApp : public CWinApp
{
public:
	CDllPositionEstimationApp();

// �������Դϴ�.
public:
	virtual BOOL InitInstance();

	DECLARE_MESSAGE_MAP()
};
