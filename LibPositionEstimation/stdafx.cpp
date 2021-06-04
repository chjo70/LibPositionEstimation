
// stdafx.cpp : 표준 포함 파일만 들어 있는 소스 파일입니다.
// LibPositionEstimation.pch는 미리 컴파일된 헤더가 됩니다.
// stdafx.obj에는 미리 컴파일된 형식 정보가 포함됩니다.

#include "pch.h"

#ifdef _DEBUG


#ifdef ATLTRACE 
#undef ATLTRACE
#undef ATLTRACE2

#define ATLTRACE CustomTrace
#define ATLTRACE2 ATLTRACE
#endif // ATLTRACE

#endif

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

void CustomTrace(const TCHAR * format, ...)
{
    const int TraceBufferSize = 1024;
    TCHAR buffer[TraceBufferSize];

    va_list argptr; va_start(argptr, format);
    sprintf_s(buffer, format, argptr);
    va_end(argptr);

    ::OutputDebugString(buffer);
}

void CustomTrace(int dwCategory, int line, const TCHAR* format, ...)
{
    va_list argptr; va_start(argptr, format);
    CustomTrace(format, argptr);
    va_end(argptr);
}