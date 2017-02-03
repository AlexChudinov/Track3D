
#include "stdafx.h"

#include <float.h>
#include "math.h" 
#include "vector2d.hpp"
#include "constant.hpp"

namespace EvaporatingParticle
{

void Vector2D::normalize()
{
  double xnorm = length();
  if(xnorm > Const_Almost_Zero )
  {
    xnorm = 1./xnorm;
    *this *= xnorm;
  }
  else
    x = y = 0;
}

double Vector2D::length() const
{
  double n = sqrt( (*this & *this) );
  if( n > Const_Almost_Zero )
    return n;
  else
    return 0.;
}

Vector2D Vector2D::rphi(const Vector2D& V)
{
  Vector2D rp(0.,0.);

  if ( (rp.x = V.length()) < Const_Almost_Zero )
    return rp;

  if( fabs(V.x) < Const_Almost_Zero )
    rp.y = ( V.y > 0. ) ? Const_Half_PI : 3.* Const_Half_PI;
  else if( fabs(V.y) < Const_Almost_Zero )
    rp.y = ( V.x > 0. ) ? 0. : Const_PI ;
  else if ( V.x > 0. && V.y > 0. )
    rp.y = atan( double(V.y/V.x) );
  else if ( V.x > 0. && V.y < 0. )
    rp.y = 2.*Const_PI - atan( double(-V.y/V.x) ) ;
  else if ( V.x < 0. && V.y > 0. )
    rp.y = Const_PI - atan( double(-V.y/V.x) );
  else if ( V.x < 0. && V.y < 0. )
    rp.y = Const_PI + atan( double(V.y/V.x) );

  rp.y = max(rp.y, (double)0);
  rp.y = min(rp.y, (double)(2.*Const_PI));

  return rp;
}

double Vector2D::operator ^ ( const Vector2D&  V2 ) const
{
  double nvv = sqrt( double( (*this & *this) * (V2 & V2) ) );
  if( -FLT_MIN <= nvv && nvv <= FLT_MIN )
    return 0.;
  double svv = (*this & V2);
  if(nvv < (svv * Const_Almost_Zero))
    return 0.;

  return acos(svv / nvv);
}

}; // namespace EvaporatingParticle