//////////////////////////////////////////////////////////////////////////
/*!
 * @file      ISNUMBER.h
 * @brief     
 * @author    ??ö?? (churlhee.jo@lignex1.com)
 * @date      2015-06-17 ???? 10:01 
 * @warning   
 */

#pragma once

#include <iostream>
#include <math.h>				
#include <float.h>
#include <limits>

using namespace std;


template<typename T>
bool is_infinite( const T &value )
{
	T max_value = (std::numeric_limits<T>::max)();
	T min_value = - max_value;

	return ! ( min_value <= value && value <= max_value );
}

template<typename T>
bool is_nan( const T &value )
{
	bool bRet=false;

	//return value != value;
	if( value > value || value < value ) {
		bRet = true;
	}

	return bRet;
}

template<typename T>
bool IsNumber( const T &value )
{
	return ! is_infinite( value ) && ! is_nan( value );
}







