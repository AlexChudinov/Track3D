#pragma once
#ifndef _Vector3D_
#define _Vector3D_

#include <float.h>
#include "math.h" 
#include "constant.hpp"

namespace EvaporatingParticle
{

//-----------------------------------------------------------------------------
//                       Vector3<class Element>
//-----------------------------------------------------------------------------
template<class Element> 
class Vector3
{
public:
  Element x, y, z;

  Vector3()
    : x(0), y(0), z(0)
  {
  }

  Vector3(const Element& X, const Element& Y, const Element& Z)
    : x(X), y(Y), z(Z)
  {
  }

  template<class Element1>
  Vector3(const Vector3<Element1>& v)
    : x(v.x), y(v.y), z(v.z)
  {
  }

  double      length() const;

// Square of length of vector.
  double      sqlength() const;

// Normalized vector. Returns V/||V|| or {0,0} if ||V||=0
  void        normalize();
  Vector3     normalized() const;

// Arithmetics operations
  Vector3     operator - () const;
  Vector3     operator + () const { return *this; }

  Vector3&    operator += (const Vector3& V2);
  Vector3&    operator -= (const Vector3& V2);

  Vector3&    operator /= (const Element& a);
  Vector3&    operator *= (const Element& a);

  Vector3     operator + (const Vector3& V2) const;
  Vector3     operator - (const Vector3& V2) const;
  Vector3     operator * (const Element& a) const;
  Vector3     operator / (const Element& a) const;

// Algebraic operations on vectors
  double      operator & (const Vector3& V2) const; // scalar product of two vectors
  Vector3     operator * (const Vector3& V2) const; // vector product of two vectors

// Angle between two vectors
  double      operator ^ (const Vector3& V2) const;

// Compare vectors
  bool        operator == (const Vector3& V2) const;
  bool        operator != (const Vector3& V2) const;
};

//-----------------------------------------------------------------------------
//                Additional functions
//-----------------------------------------------------------------------------
template<class Element>
Vector3<Element> operator * (const Element& a, const Vector3<Element>& V);

template<class Element>
double abs(const Vector3<Element>& v);

template<class Element>
void Vector3<Element>::normalize()
{
  double xnorm = length();
  if(xnorm > Const_Almost_Zero)
  {
    xnorm = 1./xnorm;
    *this *= (Element)xnorm;
  }
  else
  {
    x = y = z = 0;
  }
}

template<class Element>
double Vector3<Element>::length() const
{
  double n = sqrt( (*this & *this) );
  if( n > Const_Almost_Zero )
    return n;
  else
    return 0.;
}

template<class Element>
double Vector3<Element>::operator ^ ( const Vector3<Element>&  V2 ) const
{
  double nvv = sqrt( double( (*this & *this) * (V2 & V2) ) );
  if( -FLT_MIN <= nvv && nvv <= FLT_MIN )
    return 0.;
  double svv = (*this & V2);
  if(nvv < (svv * Const_Almost_Zero))
    return 0.;

  return acos(svv / nvv);
}

//-----------------------------------------------------------------------------
//                Inline functions
//-----------------------------------------------------------------------------
template<class Element>
inline Vector3<Element> Vector3<Element>::operator - () const
{ return Vector3<Element>( -x, -y, -z ); }

template<class Element>
inline Vector3<Element>& Vector3<Element>::operator += (const Vector3<Element>& V2)
{ x += V2.x;  y += V2.y;  z += V2.z; return *this; }

template<class Element>
inline Vector3<Element>& Vector3<Element>::operator -= (const Vector3<Element>& V2)
{ x -= V2.x;  y -= V2.y;  z -= V2.z; return *this; }

template<class Element>
inline Vector3<Element>& Vector3<Element>::operator /= (const Element& a)
{ x /= a;  y /= a;  z /= a; return *this; }

template<class Element>
inline Vector3<Element>& Vector3<Element>::operator *= (const Element& a)
{ x *= a;  y *= a;  z *= a; return *this; }

template<class Element>
inline Vector3<Element> Vector3<Element>::operator + (const Vector3<Element>&  V2) const
{ return  Vector3<Element>(x + V2.x, y + V2.y, z + V2.z); }

template<class Element>
inline Vector3<Element> Vector3<Element>::operator - (const Vector3<Element>&  V2) const
{ return  Vector3<Element>(x - V2.x, y - V2.y, z - V2.z); }

template<class Element>
inline Vector3<Element> Vector3<Element>::operator * (const Element& a) const
{ return  Vector3<Element>(x * a, y * a, z * a); }

template<class Element>
inline Vector3<Element> operator * (const Element& a, const Vector3<Element>& V)
{ return V * a; }

template<class Element>
inline Vector3<Element> Vector3<Element>::operator / (const Element& a) const
{ return Vector3<Element>(x / a, y / a, z / a); }

template<class Element>
inline double Vector3<Element>::operator & (const Vector3<Element>& V2) const
{ return x * V2.x + y * V2.y + z * V2.z; }

template<class Element>
inline Vector3<Element> Vector3<Element>::operator * (const Vector3<Element>& V2) const
{ return Vector3<Element>(y * V2.z - z * V2.y, z * V2.x - x * V2.z, x * V2.y - y * V2.x); }

template<class Element>
inline Vector3<Element> Vector3<Element>::normalized() const
{ Vector3<Element> v(*this); v.normalize(); return v; }

template<class Element>
inline double Vector3<Element>::sqlength() const
{ return (*this & *this); }

template<class Element>
inline bool Vector3<Element>::operator == ( const Vector3<Element>& V2 ) const
{ return (x == V2.x) && (y == V2.y) && (z == V2.z); }

template<class Element>
inline bool Vector3<Element>::operator != ( const Vector3<Element>& V2 ) const
{ return !(*this == V2); }

template<class Element>
inline double abs(const Vector3<Element>& v)
{ return v.length(); }

/**
 * Elementwise products of two vectors
 * AC 24/06/2016
 */
template<class Element>
inline Vector3<Element> operator && (const Vector3<Element>& v1, const Vector3<Element>& v2)
{
	return Vector3<Element>(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}

typedef Vector3<double> Vector3D;
typedef Vector3<float> Vector3F;

};  // namespace EvaporatingParticle

#endif // #ifndef _Vector3D_
