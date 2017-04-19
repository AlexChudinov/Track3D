
#pragma once

#include "CObject.h"
#include "Elements.h"

namespace EvaporatingParticle
{

struct RGB_Color
{
  RGB_Color(unsigned char r = 0, unsigned char g = 0, unsigned char b = 0)
    : red(r), green(g), blue(b)
  {
  }

  RGB_Color(COLORREF clr)
    : red(GetRValue(clr)), green(GetGValue(clr)), blue(GetBValue(clr))
  {
  }

  unsigned char red,
                green,
                blue;
};

typedef std::vector<double> CValueVector;
typedef std::vector<RGB_Color> CColorVector;
//---------------------------------------------------------------------------------------
// CColorImage - a base class for countours, vector plots.
//---------------------------------------------------------------------------------------
class CColorImage : public CObject
{
public:
  CColorImage();
  virtual ~CColorImage();

  enum  // Scale type
  {
    scLinear = 0,
    scLog    = 1,
    scCount  = 2
  };

  enum  // Color map
  {
    cmRainbow     = 0,  // 4 intervals: red - yellow, yellow - green, green - cyan, cyan - blue.
    cmRainbowExt  = 1,  // 5 intervals: magenta - red, red - yellow, yellow - green, green - cyan, cyan - blue.
    cmCount       = 2
  };

  virtual void        draw() = 0;

  bool                get_enable_image() const;
  DWORD_PTR           get_enable_image_ptr() const;
  void                set_enable_image(bool bEnable);

  int                 get_scale_type() const;
  DWORD_PTR           get_scale_type_ptr() const;
  void                set_scale_type(int nType);

  UINT                get_levels_count() const;
  DWORD_PTR           get_levels_count_ptr() const;
  void                set_levels_count(UINT nCount);

  int                 get_color_map_type() const;
  DWORD_PTR           get_color_map_type_ptr() const;
  void                set_color_map_type(int nType);

  bool                get_enable_user_range() const;
  DWORD_PTR           get_enable_user_range_ptr() const;
  void                set_enable_user_range(bool bEnable);

  CString             get_drawn_reg_names() const;
  DWORD_PTR           get_drawn_reg_names_ptr() const;

  void                clear_reg_names();
  void                invalidate();

  virtual void        save(CArchive& archive);
  virtual void        load(CArchive& archive);

  virtual void        restore_user_range() = 0;
  virtual void        get_min_max() = 0;

  const char*         get_scale_name(int nScaleType) const;
  const char*         get_clr_map_name(int nClrMapType) const;

protected:
  virtual void        set_default();

  virtual bool        build_draw_array() = 0;

  CRegionsCollection  get_regions() const;

  void                get_values_array();

  RGB_Color           get_color(double fVal) const;
  void                get_colors_array();

  unsigned char       clamp(double fColorComp) const;

// User interface:
  bool                m_bEnable;      // enable to draw the image.

  bool                m_bUserDefRange;

  double              m_fMinVal,
                      m_fMaxVal;

  UINT                m_nLevelsCount; // inner levels count; the intervals count is m_nLevelsCount + 1.

  int                 m_nColorMapType,
                      m_nScaleType;   // the scale type can be either linear or logarithmic.

  CStringVector       m_vRegNames;

// Run-time:
  CValueVector        m_vValues;  // size of this is m_nLevelsCount, the count of inner values.
  CColorVector        m_vColors;  // size of this is m_nLevelsCount + 1, count of colored intervals.

  bool                m_bReady;
};

//---------------------------------------------------------------------------------------
// Inline implementation.
//---------------------------------------------------------------------------------------
inline bool CColorImage::get_enable_image() const
{
  return m_bEnable;
}

inline DWORD_PTR CColorImage::get_enable_image_ptr() const
{
  return (DWORD_PTR)&m_bEnable;
}

inline void CColorImage::set_enable_image(bool bEnable)
{
  m_bEnable = bEnable;
}

inline int CColorImage::get_scale_type() const
{
  return m_nScaleType;
}

inline DWORD_PTR CColorImage::get_scale_type_ptr() const
{
  return (DWORD_PTR)&m_nScaleType;
}

inline void CColorImage::set_scale_type(int nType)
{
  if(m_nScaleType != nType)
  {
    m_nScaleType = nType;
    m_bReady = false;
  }
}

inline UINT CColorImage::get_levels_count() const
{
  return m_nLevelsCount;
}

inline DWORD_PTR CColorImage::get_levels_count_ptr() const
{
  return (DWORD_PTR)&m_nLevelsCount;
}

inline void CColorImage::set_levels_count(UINT nCount)
{
  if(m_nLevelsCount != nCount)
  {
    m_nLevelsCount = nCount;
    m_bReady = false;
  }
}

inline int CColorImage::get_color_map_type() const
{
  return m_nColorMapType;
}

inline DWORD_PTR CColorImage::get_color_map_type_ptr() const
{
  return (DWORD_PTR)&m_nColorMapType;
}

inline void CColorImage::set_color_map_type(int nType)
{
  if(m_nColorMapType != nType)
  {
    m_nColorMapType = nType;
    m_bReady = false;
  }
}

inline bool CColorImage::get_enable_user_range() const
{
  return m_bUserDefRange;
}

inline DWORD_PTR CColorImage::get_enable_user_range_ptr() const
{
  return (DWORD_PTR)&m_bUserDefRange;
}

inline CString CColorImage::get_drawn_reg_names() const
{
  return CObject::compile_string(m_vRegNames);
}

inline DWORD_PTR CColorImage::get_drawn_reg_names_ptr() const
{
  return (DWORD_PTR)&m_vRegNames;
}

inline void CColorImage::clear_reg_names()
{
  if(m_vRegNames.size() != 0)
  {
    m_vRegNames.clear();
    m_bReady = false;
  }
}

inline void CColorImage::invalidate()
{
  m_bReady = false;
}

inline unsigned char CColorImage::clamp(double fColorComp) const
{
  return fColorComp < 0 ? 0 : (fColorComp > 255 ? 255 : (unsigned char)fColorComp);
}

inline const char* CColorImage::get_scale_name(int nScaleType) const
{
  switch(nScaleType)
  {
    case scLinear: return _T("Linear");
    case scLog:    return _T("Logarithmic");
  }

  return _T("");
}

inline const char* CColorImage::get_clr_map_name(int nClrMapType) const
{
  switch(nClrMapType)
  {
    case cmRainbow:     return _T("Rainbow");
    case cmRainbowExt:  return _T("Extended Rainbow");
  }

  return _T("");
}

}; // namespace EvaporatingParticle
