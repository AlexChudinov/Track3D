#pragma once
#ifndef _EXCEPT_
#define _EXCEPT_

class CBase;

/**
@brief Exception is an exception handler class
AC 23/06/2016
*/
class Exception : public std::exception
{
	const std::string m_String;
public:
	Exception(const char* String) throw();
	Exception(const std::string& str) throw();
	Exception(const Exception& ex) throw();
	Exception(const CBase& baseObj) throw();
	~Exception() throw();
	const char* what() const throw();
};

#endif

