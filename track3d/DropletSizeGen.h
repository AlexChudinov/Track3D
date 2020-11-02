#pragma once

#include "intsafe.h"
#include <random>
#include <vector>

class CArchive;

namespace EvaporatingParticle
{

class CDropletSizeGen
{
public:
  CDropletSizeGen();

  enum
  {
    dstConst    = 0,
    dstUniform  = 1,
    dstGauss    = 2,
    dstCount    = 3
  };

  void                  generate();
  double                get_droplet_size(size_t i) const;

  int                   get_distr_type() const;
  DWORD_PTR             get_distr_type_ptr() const;
  void                  set_distr_type(int nType);

  size_t                get_droplets_count() const;
  DWORD_PTR             get_droplets_count_ptr() const;
  void                  set_droplets_count(size_t nCount);

  double                get_mean_size() const;
  DWORD_PTR             get_mean_size_ptr() const;
  void                  set_mean_size(double fSize);

  double                get_std_dev() const;
  DWORD_PTR             get_std_dev_ptr() const;
  void                  set_std_dev(double fStdDev);

  double                get_min_size() const;
  DWORD_PTR             get_min_size_ptr() const;
  void                  set_min_size(double fSize);

  double                get_max_size() const;
  DWORD_PTR             get_max_size_ptr() const;
  void                  set_max_size(double fSize);

  unsigned long long    get_rand_seed() const;
  DWORD_PTR             get_rand_seed_ptr() const;
  void                  set_rand_seed(unsigned long long nSeed);

  static const char*    distr_name(int nType);

// Streams support:
  void                  save(CArchive& ar);
  void                  load(CArchive& ar);

protected:
  void                  set_default();
  void                  debug_output();

private:
  std::mt19937_64       m_RndGen;
  unsigned long long    m_nRandSeed;

  std::vector<double>   m_vSizeArr;

  double                m_fMeanSize,
                        m_fStdDev,

                        m_fMinSize,
                        m_fMaxSize;

  int                   m_nDistrType;
};

// Inline implementation
inline int CDropletSizeGen::get_distr_type() const
{
  return m_nDistrType;
}

inline DWORD_PTR CDropletSizeGen::get_distr_type_ptr() const
{
  return (DWORD_PTR)&m_nDistrType;
}

inline void CDropletSizeGen::set_distr_type(int nType)
{
  m_nDistrType = nType;
}

inline size_t CDropletSizeGen::get_droplets_count() const
{
  return m_vSizeArr.size();
}

inline DWORD_PTR CDropletSizeGen::get_droplets_count_ptr() const
{
  return (DWORD_PTR)&m_vSizeArr;
}

inline void CDropletSizeGen::set_droplets_count(size_t nCount)
{
  m_vSizeArr.resize(nCount, 0.0);
}

inline double CDropletSizeGen::get_mean_size() const
{
  return m_fMeanSize;
}

inline DWORD_PTR CDropletSizeGen::get_mean_size_ptr() const
{
  return (DWORD_PTR)&m_fMeanSize;
}

inline void CDropletSizeGen::set_mean_size(double fSize)
{
  m_fMeanSize = fSize;
}

inline double CDropletSizeGen::get_std_dev() const
{
  return m_fStdDev;
}

inline DWORD_PTR CDropletSizeGen::get_std_dev_ptr() const
{
  return (DWORD_PTR)&m_fStdDev;
}

inline void CDropletSizeGen::set_std_dev(double fStdDev)
{
  m_fStdDev = fStdDev;
}

inline double CDropletSizeGen::get_min_size() const
{
  return m_fMinSize;
}

inline DWORD_PTR CDropletSizeGen::get_min_size_ptr() const
{
  return (DWORD_PTR)&m_fMinSize;
}

inline void CDropletSizeGen::set_min_size(double fSize)
{
  m_fMinSize = fSize;
}

inline double CDropletSizeGen::get_max_size() const
{
  return m_fMaxSize;
}

inline DWORD_PTR CDropletSizeGen::get_max_size_ptr() const
{
  return (DWORD_PTR)&m_fMaxSize;
}

inline void CDropletSizeGen::set_max_size(double fSize)
{
  m_fMaxSize = fSize;
}

inline unsigned long long CDropletSizeGen::get_rand_seed() const
{
  return m_nRandSeed;
}

inline DWORD_PTR CDropletSizeGen::get_rand_seed_ptr() const
{
  return (DWORD_PTR)&m_nRandSeed;
}

inline void CDropletSizeGen::set_rand_seed(unsigned long long nSeed)
{
  m_nRandSeed = nSeed;
}

};  // namespace EvaporatingParticle
