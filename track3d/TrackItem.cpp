
#include "stdafx.h"
#include "TrackItem.h"
//#include "Integrators.hpp"
//#include "BarnesHut.h"
//#include "StochasticGasDynamic.hpp"
//#include "Tracker.hpp"
//#include "ParticleTracking.h"
//#include <algorithm>


namespace EvaporatingParticle
{

CBaseTrackItem* CBaseTrackItem::create(int nType)
{
  switch(nType)
  {
    case CTrack::ptDroplet:
    {
      CDropletTrackItem* pDropletItem = new CDropletTrackItem();
      return (CBaseTrackItem*)pDropletItem;
    }
    case CTrack::ptIon:
    {
      CIonTrackItem* pIonItem = new CIonTrackItem();
      return (CBaseTrackItem*)pIonItem;
    }
  }

  return NULL;
}

//---------------------------------------------------------------------------------------
//
// CDropletTrackItem:
// pState is an array of size DROPLET_STATE_SIZE = 8.
// pItem->pos(pState[0], pState[1], pState[2]);
// pItem->vel(pState[3], pState[4], pState[5]);
// pItem->temp = pState[6];
// pItem->mass = pState[7];
//
// CIonTrackItem:
// pState is an array of size ION_STATE_SIZE = 8.
// pItem->pos(pState[0], pState[1], pState[2]);
// pItem->vel(pState[3], pState[4], pState[5]);
// pItem->temp = pState[6];
// pItem->unfragm = pState[7];
//
//---------------------------------------------------------------------------------------
CBaseTrackItem* CBaseTrackItem::create(int nType, int nElemId, double fTime, double* pState)
{
  if(pState == NULL)
    return NULL;

  Vector3D vPos(pState[0], pState[1], pState[2]);
  Vector3D vVel(pState[3], pState[4], pState[5]);
  switch(nType)
  {
    case CTrack::ptDroplet:
    {
      double fTemp = pState[6];
      double fMass = pState[7];
      CDropletTrackItem* pDropletItem = new CDropletTrackItem(nElemId, vPos, vVel, fTemp, fMass, fTime);

      return (CBaseTrackItem*)pDropletItem;
    }
    case CTrack::ptIon:
    {
      double fIonTemp = pState[6];
      double fUnfragm = pState[7];
      CIonTrackItem* pIonItem = new CIonTrackItem(nElemId, vPos, vVel, fIonTemp, fUnfragm, fTime);

      return (CBaseTrackItem*)pIonItem;
    }
  }
  
  return NULL;  
}

void CBaseTrackItem::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  ar << nElemId;
  ar << pos.x;
  ar << pos.y;
  ar << pos.z;
  ar << vel.x;
  ar << vel.y;
  ar << vel.z;
  ar << time;
}

void CBaseTrackItem::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  ar >> nElemId;
  ar >> pos.x;
  ar >> pos.y;
  ar >> pos.z;
  ar >> vel.x;
  ar >> vel.y;
  ar >> vel.z;
  ar >> time;
}

//---------------------------------------------------------------------------------------
//  CIonTrackItem
//---------------------------------------------------------------------------------------
void CIonTrackItem::state(double* pState) const
{
  pState[0] = pos.x;
  pState[1] = pos.y;
  pState[2] = pos.z;
  pState[3] = vel.x;
  pState[4] = vel.y;
  pState[5] = vel.z;
  pState[6] = temp;
  pState[7] = unfragm;
}

CBaseTrackItem* CIonTrackItem::copy() const
{
  CIonTrackItem* pItem = new CIonTrackItem(nElemId, pos, vel, temp, unfragm, time);
  pItem->tempinf = tempinf;
  return (CBaseTrackItem*)pItem;
}

double CIonTrackItem::get_mass() const
{
  return 0.;
}

double CIonTrackItem::get_temp() const
{
  return temp;
}

void CIonTrackItem::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  CBaseTrackItem::save(ar);

  ar << temp;
  ar << tempinf;
  ar << unfragm;
}

void CIonTrackItem::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CBaseTrackItem::load(ar);

  ar >> temp;
  ar >> tempinf;
  ar >> unfragm;
}

//---------------------------------------------------------------------------------------
//  CDropletTrackItem
//---------------------------------------------------------------------------------------
void CDropletTrackItem::state(double* pState) const
{
  pState[0] = pos.x;
  pState[1] = pos.y;
  pState[2] = pos.z;
  pState[3] = vel.x;
  pState[4] = vel.y;
  pState[5] = vel.z;
  pState[6] = temp;
  pState[7] = mass;
}

CBaseTrackItem* CDropletTrackItem::copy() const
{
  CDropletTrackItem* pItem = new CDropletTrackItem(nElemId, pos, vel, temp, mass, time);
  return (CBaseTrackItem*)pItem;
}

double CDropletTrackItem::get_mass() const
{
  return mass;
}

double CDropletTrackItem::get_temp() const
{
  return temp;
}

void CDropletTrackItem::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  CBaseTrackItem::save(ar);

  ar << temp;
  ar << mass;
}

void CDropletTrackItem::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CBaseTrackItem::load(ar);

  ar >> temp;
  ar >> mass;
}

//---------------------------------------------------------------------------------------
// CTrack
//---------------------------------------------------------------------------------------
void CTrack::get_track_item(size_t nInd, CTrackItem& item) const
{
  CBaseTrackItem* pItem = at(nInd);
  item.nElemId = pItem->nElemId;
  item.pos = pItem->pos;
  item.vel = pItem->vel;
  item.time = pItem->time;

  switch(m_nType)
  {
    case ptDroplet:
    {
      CDropletTrackItem* pDropletItem = (CDropletTrackItem*)pItem;
      item.temp = pDropletItem->temp;
      item.mass = pDropletItem->mass;
      item.tempinf = 0;
      item.unfragm = 1;
    }
    case ptIon:
    {
      CIonTrackItem* pIonItem = (CIonTrackItem*)pItem;
      item.temp = pIonItem->temp;
      item.tempinf = pIonItem->tempinf;
      item.unfragm = pIonItem->unfragm;
      item.mass = 0;
    }
  }
}

};  // namespace EvaporatingParticle

#if 0

/**
 * AC 28/06/2016
 * modified by AC 06/07/2016
 */

using namespace EvaporatingParticle;

BasePathIntegrator::operator BasePathIntegrator::_DY() const
{
	//Here is the start values set for an every integration step
	_DY result;

	result. vec3D<1>() = item().vel;
	result. vec3D<2>() = item().pos;

	result.scalar<1>() = item().temp;
	result.scalar<2>() = item().diffCoef; 
	result.scalar<3>() = item().stationaryTemperature; //[AC] 22/07/2016
	result.scalar<4>() = item().conc; //[AC] 11/08/2016

	return result;
}

BasePathIntegrator & BasePathIntegrator::operator=(const _DY & dy)
{//Here is we return the values back from the integrator
	item().vel = dy.vec3D<1>();
	item().pos = dy.vec3D<2>();

	item().temp = dy.scalar<1>();
	item().diffCoef = dy.scalar<2>();
	item().stationaryTemperature = dy.scalar<3>(); //[AC] 22/07/2016
	item().conc = dy.scalar<4>(); //[AC] 11/08/2016

	return *this;
}

std::size_t IonMobilityPathIntegrator::switchIntegrator() const
{
	double maxFreq = m_field->fieldHighestFrequency();

	//Set current mesh point and calculate field values
	if(!(m_isExist = m_field->setNode(item().pos)))
		return 0; //If not existing, just return any valuable
	Vd scalarFieldVals  = m_field->scalarFields(item().pos, item().time);
	
	//Calculate current pressure and temperature
	const double& T = scalarFieldVals[2];
	const double& P = scalarFieldVals[1];

	//current ion mobility
	double kappa =
		item().kappa * T / std::sqrt(Const_T0*item().temp)
		* Const_One_Atm_CGS / P;

	if (item().charge * m_h / (kappa*item().mass) < 1e-2)
	{
		return LOW_PRESS;
	}
	else
	{
		return HIGH_PRESS;
	}
}

BasePathIntegrator::FieldValues 
BasePathIntegrator::fieldValues(const V3D& pos, double time) const
{
	FieldValues fieldValues;

	if(!(m_isExist = m_field->setNode(item().pos)))
		return fieldValues;

	double rfTime = time - item().phase;
	Vd scalarFieldValues  = m_field->scalarFields(pos, rfTime);
// [MS] 20-07-2016 added the ion velocity and current to support the Coulomb field in the axially-symmetric case.
  VV3D vectorFieldValues = m_field->vectorFields(pos, item().vel, rfTime, item().curr);

	fieldValues.E = vectorFieldValues[0];
	fieldValues.U = vectorFieldValues[1];
	fieldValues.P = scalarFieldValues[1];
	fieldValues.T = scalarFieldValues[2];
	fieldValues.Cv= scalarFieldValues[3] * item().gasMolMass - Const_Boltzmann;

	return fieldValues;
}

BasePathIntegrator::BasePathIntegrator(const CTracker * tracker)
	:
	m_field(new AnsysField(tracker)),
	m_isExist(true)
{
	//Choose integrator type
	switch (tracker->integratorType())
	{
	case CTracker::RungeKutta2:
		m_integrator = 
			new math::IntegratorRK2 <BasePathIntegrator, _DY, double>;
		break;
	case CTracker::RungeKutta4:
		m_integrator =
			new math::IntegratorRK4 <BasePathIntegrator, _DY, double>;
		break;
	case CTracker::PredictorCorrector:
		m_integrator =
			new math::IntegratorPC <BasePathIntegrator, _DY, double>;
		break;
	default:
		throw Exception("Unsupported integrator type!");
		break;
	}
}

BasePathIntegrator::~BasePathIntegrator()
{
	delete m_field;
	delete m_integrator;
}

IonMobilityPathIntegrator::IonMobilityPathIntegrator
(
	const CTracker * tracker,
	double h,
	double maxTime
) :
	BasePathIntegrator(tracker),
	m_h(h),
	m_maxTime(maxTime)
{
	_DFuns dFuns(2);
	dFuns[0] = &IonMobilityPathIntegrator::diffLP<math::MobilityCalculatorUsingGasTemperature<double>>;
	dFuns[1] = &IonMobilityPathIntegrator::diffHP<math::MobilityCalculatorUsingGasTemperature<double>>;

	m_switchDiffScheme =
		new math::SwitchDiffScheme
		<
			IonMobilityPathIntegrator,
			_DY, double
		>
		(
			*this, dFuns, item().time,
			m_h, m_integrator,
			&IonMobilityPathIntegrator::switchIntegrator
		);
}

CTrack EvaporatingParticle::IonMobilityPathIntegrator::estimateTrack
(
	int timeStepsPerSample, 
	const CTrackItem & startConditions,
	const bool * terminateFlag
)
{
	CTrack track;
	track.reserve(100);

	//keep first value
	track.push_back(startConditions);

	item() = startConditions;

	int nSteps = 0;
	
	while (!(*terminateFlag) && m_isExist && maxTime() > item().time)
	{
		m_switchDiffScheme->doIntegratorStep();
		if ((++nSteps) == timeStepsPerSample && m_isExist)
		{
			nSteps = 0;
			track.push_back(item());
		}
	}

	return track;
}

AnsysMesh::AnsysMesh(const CTracker * mesh)
	:
	Mesh<CTracker, CElem3D, CNode3D>
	(
		mesh,
		&CElem3D::inside,
		&CElem3D::getNeighborElems,
		&CElem3D::fieldScalarVals,
		&CElem3D::fieldVectorVals,
		&CElem3D::idx,
		&CTracker::searchElemFromIdx,
		&CTracker::searchElemFromBegin,
		&CTracker::get_symmetry_type
	)
{
}

AnsysField::AnsysField(const CTracker * tracker)
	:
	m_mesh(new AnsysMesh(tracker)),
	m_tracker(tracker)
{
}

AnsysField::~AnsysField()
{
	delete m_mesh;
}

bool AnsysField::setNode(const V3D& pos)
{
	return m_mesh->setNode(pos);
}

// [MS] 20-07-2016 added the ion velocity to support the Coulomb field in the axially-symmetric case.
VV3D AnsysField::vectorFields(const V3D& pos, const V3D& vel, double time, double curr) const
{
	VV3D field = m_mesh->getVectorVals();
	VV3D fieldVals(2); //0-electric field, 1-gas flow

	fieldVals[0] = m_tracker->rfField(field[1], pos, time)
		+ m_tracker->dcField(field[0], pos) + m_tracker->spaceChargeField(pos, vel, curr);

	fieldVals[1] = field[2]; //Gas velocity field

	return fieldVals;
}

Vd AnsysField::scalarFields(const V3D & /*pos*/, double /*time*/) const
{
	return m_mesh->getScalarVals();
}

double AnsysField::fieldHighestFrequency() const
{
	return std::max<double>
		(
			std::max<double>
			(
				m_tracker->get_rf_frequency(), m_tracker->get_rf_Q00_freq()
			),
			m_tracker->get_rf_flatapole_freq()
		);
}

double EvaporatingParticle::AnsysField::fragmentationEnergy() const
{
	return m_tracker->get_act_energy() / Const_Erg_to_EV;
}

double EvaporatingParticle::AnsysField::ionCollisionCrossection() const
{
	return m_tracker->get_ion_cross_section();
}

IonMobilityPathIntegratorDiffusionJumps::
IonMobilityPathIntegratorDiffusionJumps
(
	const CTracker * tracker, 
	double h,
	double maxTime
)
	:IonMobilityPathIntegrator(tracker,h,maxTime)
{
}

CTrack IonMobilityPathIntegratorDiffusionJumps::
estimateTrack
(
	int timeStepsPerSample,
	const CTrackItem & startConditions,
	const bool * terminateFlag
)
{
	CTrack track;
	track.reserve(100);

	//keep first value
	track.push_back(startConditions);

	item() = startConditions;

	int nSteps = 0;
	while (!(*terminateFlag) && m_isExist && maxTime() > item().time)
	{
		differScheme()->doIntegratorStep();

		//Do diffusion jump
		item().pos += std::sqrt(2 * item().diffCoef*h())*math::randOnSphere();

		if ((++nSteps) == timeStepsPerSample && m_isExist)
		{
			nSteps = 0;
			track.push_back(item());
		}
	}

	return track;
}

#endif  // #if 0