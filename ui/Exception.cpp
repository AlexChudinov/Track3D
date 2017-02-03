/**
 * AC 23/06/2016
 */
#include <stdafx.h>
#include <exception>
#include "Exception.h"
#include "Base.h"

Exception::Exception(const char* String) throw()
	:
	m_String(String)
{
}

Exception::Exception(const std::string & str) throw() : m_String(str)
{
	;
}

Exception::Exception(const Exception& ex) throw()
	:
	m_String(ex.what())
{
}

Exception::Exception(const CBase & baseObj) throw() : m_String(baseObj.toString())
{
	;
}

Exception::~Exception() throw()
{
}

const char* Exception::what() const throw()
{
	return m_String.c_str();
}