#pragma once

#include "CObject.h"
#include "Elements.h"

namespace EvaporatingParticle
{

//-------------------------------------------------------------------------------------------------
// CImportOpenFOAM a class for reading data from OpenFOAM file.
//-------------------------------------------------------------------------------------------------
class CImportOpenFOAM : public CObject
{
public:
  CImportOpenFOAM();
  ~CImportOpenFOAM();

  enum  // Note: Import must be performed in two steps, 1 and 2.
  {
// Preliminary step. Reading OpenFOAM file, searching for proper nodes in a short mesh used for OpenFOAM export and changing corresponding
// *.var file on disk. Input: *.geom file of a short mesh as the "Gas-Dynamic data file" and a *.csv OpenFOAM file as the "Import data file".
// At this step the proper nodes are searched by a simple comparison of the distance with a given tolerance eps ~ 0.01 mm.
    opFirstStep   = 0,

// Run-time step. Input: *.geom file of a full mesh sutable for tracking as the "Gas-Dynamic data file" and a modified *.geom file of the 
// short mesh modified at the first step as the "Import data file". At this step OpenFOAM values at the nodes of the full mesh are interpolated
// from the short mesh.
    opSecondStep  = 1,

// Run-time operation. The original ANSYS values are read from the *.var file of the full mesh (untouched by the import operations).
    opRestore     = 2
  };

  bool                do_import();

  int                 get_step() const;
  DWORD_PTR           get_step_ptr() const;
  void                set_step(int nStep);

  const char*         get_filename() const;
  DWORD_PTR           get_filename_ptr() const;
  void                set_filename(const char* pName);

  double              get_x_shift() const;
  DWORD_PTR           get_x_shift_ptr() const;
  void                set_x_shift(double fShift);

  void                save(CArchive& ar);
  void                load(CArchive& ar);

protected:
  void                set_default();

// First step (reading from the OpenFOAM file, identifying nodes and updating the file of variables of the "short" mesh on the disk):
  bool                first_step();

  CNode3D*            find_node(float fX, float fY, float fZ) const;

  void                convert_to_cgs(float& fN, float& fVx, float& fVy, float& fVz, float& fX, float& fY, float& fZ) const;
  void                convert_to_si(float& fPress, float& fDens, float& fVx, float& fVy, float& fVz) const;

  void                print_statistics(UINT nSuccess, UINT nFail);
  bool                update_var_file();

// Second step, run-time (reading the *.geom, *.var and *.rgn files of the "short" mesh, updated at the previous step, interpolating and
// substituting the pressure, temperature and velocity fields in the nodes of the full mesh):
  bool                second_step();

// Restore the original gas-dynamic data:
  bool                restore_ansys_data();

private:
// User-defined data:
  int                 m_nStep;    // can be either opFirstStep or opSecondStep.
  std::string         m_sFile;

  double              m_fShiftX;      // cm.
};

inline int CImportOpenFOAM::get_step() const
{
  return m_nStep;
}

inline DWORD_PTR CImportOpenFOAM::get_step_ptr() const
{
  return (DWORD_PTR)&m_nStep;
}

inline void CImportOpenFOAM::set_step(int nStep)
{
  m_nStep = nStep;
}

inline const char* CImportOpenFOAM::get_filename() const
{
  return m_sFile.c_str();
}

inline DWORD_PTR CImportOpenFOAM::get_filename_ptr() const
{
  return (DWORD_PTR)&m_sFile;
}

inline void CImportOpenFOAM::set_filename(const char* pName)
{
  m_sFile = pName;
}

inline double CImportOpenFOAM::get_x_shift() const
{
  return m_fShiftX;
}

inline DWORD_PTR CImportOpenFOAM::get_x_shift_ptr() const
{
  return (DWORD_PTR)&m_fShiftX;
}

inline void CImportOpenFOAM::set_x_shift(double fShift)
{
  m_fShiftX = fShift;
}

};  // namespace EvaporatingParticle