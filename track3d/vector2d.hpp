#ifndef _Vector2D_
#define _Vector2D_

namespace EvaporatingParticle
{

//-----------------------------------------------------------------------------
//                       Vector2D
//-----------------------------------------------------------------------------
class Vector2D
{
public:
  double x, y;

  Vector2D()
    : x(0), y(0)
  {
  }

  Vector2D(double X, double Y)
    : x(X), y(Y)
  {
  }

  Vector2D(const Vector2D& v)
    : x(v.x), y(v.y)
  {
  }

  double      length() const;

// Square of length of vector.
  double      sqlength() const;

// Normalized vector. Returns V/||V|| or {0,0} if ||V||=0
  void        normalize();
  Vector2D    normalized() const;

// Polar coordinates
  static Vector2D rphi(const Vector2D& V);

//Arithmetics operations
  Vector2D    operator - () const;
  Vector2D    operator + () const { return *this; }

  Vector2D&   operator += (const Vector2D& V2);
  Vector2D&   operator -= (const Vector2D& V2);

  Vector2D&   operator /= (const double a);
  Vector2D&   operator *= (const double a);

  Vector2D    operator + (const Vector2D& V2) const;
  Vector2D    operator - (const Vector2D& V2) const;
  Vector2D    operator * (const double  a) const;
  Vector2D    operator / (const double  a) const;

// Algebraic operations on vectors
  double      operator & (const Vector2D& V2) const; // scalar product of two vectors

// Angle between two vectors
  double      operator ^ (const Vector2D& V2) const;

//Compare vectors
  bool        operator == (const Vector2D& V2) const;
  bool        operator != (const Vector2D& V2) const;
};

//-----------------------------------------------------------------------------
//                Additional functions
//-----------------------------------------------------------------------------

Vector2D operator * ( const double  a, const Vector2D&  V );

double abs( const Vector2D& v );

//-----------------------------------------------------------------------------
//                Inline functions
//-----------------------------------------------------------------------------
inline Vector2D Vector2D::operator - () const
{ return Vector2D( -x, -y ); }

inline Vector2D& Vector2D::operator += (const Vector2D& V2)
{ x += V2.x;  y += V2.y;  return *this; }

inline Vector2D& Vector2D::operator -= (const Vector2D& V2)
{ x -= V2.x;  y -= V2.y;  return *this; }

inline Vector2D& Vector2D::operator /= (const double a)
{ x /= (double)a;  y /= (double)a;  return *this; }

inline Vector2D& Vector2D::operator *=( const double a )
{ x *= (double)a;  y *= (double)a;  return *this; }

inline Vector2D Vector2D::operator + (const Vector2D&  V2) const
{ return  Vector2D(x + V2.x, y + V2.y); }

inline Vector2D Vector2D::operator - (const Vector2D&  V2) const
{ return  Vector2D( x - V2.x, y - V2.y ); }

inline Vector2D Vector2D::operator * (const double  a) const
{ return  Vector2D( (double)(x * a), (double)(y * a) ); }

inline Vector2D operator * ( const double a, const Vector2D& V )
{ return V * a; }

inline Vector2D Vector2D::operator / ( const double  a ) const
{ return Vector2D( (double)(x / a), (double)(y / a) ); }

inline double Vector2D::operator & (const Vector2D& V2) const
{ return (double)x * V2.x + (double)y * V2.y; }

inline Vector2D Vector2D::normalized() const
{ Vector2D v(*this); v.normalize(); return v; }

inline double Vector2D::sqlength() const
{ return (*this & *this); }

inline bool Vector2D::operator == ( const Vector2D& V2 ) const
{ return x == V2.x && y == V2.y; }

inline bool Vector2D::operator != ( const Vector2D& V2 ) const
{ return !(*this == V2); }

inline double abs(const Vector2D& v)
{ return v.length(); }

};  // namespace EvaporatingParticle

#endif // #ifndef _Vector2D_
