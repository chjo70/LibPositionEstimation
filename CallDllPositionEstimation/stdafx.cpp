
// stdafx.cpp : 표준 포함 파일만 들어 있는 소스 파일입니다.
// CallDllPositionEstimation.pch는 미리 컴파일된 헤더가 됩니다.
// stdafx.obj에는 미리 컴파일된 형식 정보가 포함됩니다.

#include "stdafx.h"


CCriticalSection g_criticalExe;


CString GetFilePath()
{
	g_criticalExe.Lock();

	static TCHAR pBuf[256] = { 0, } ;

	memset( pBuf, NULL, sizeof(pBuf) );

	GetModuleFileName( NULL, pBuf, sizeof(pBuf) );

	CString strFilePath;

	strFilePath.Format( _T("%s"), pBuf );
	strFilePath = strFilePath.Left( strFilePath.ReverseFind( _T('\\') ) );

	g_criticalExe.Unlock();

	return strFilePath;

}