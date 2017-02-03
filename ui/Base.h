#pragma once
#ifndef _CBASE_
#define _CBASE_

#include <sstream>

/**
 * @ Base class to implement some basic io operations
 * AC 23/06/2016
 */

class CBase
{
public:
	CBase() { ; }
	virtual ~CBase() { ; }

	/**
	 * @brief Every class can be represented in a string form
	 */
	virtual std::string toString() const 
	{ return std::string("This is a simple base class!\n"); }

	/**
	 * @brief Every class can be printed on the screen
	 */
	virtual void print(std::ostream& out) const { out << toString(); }
};

#endif // !_BASE_

