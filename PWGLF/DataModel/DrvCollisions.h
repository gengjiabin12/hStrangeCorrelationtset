
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace drvcollision
{
DECLARE_SOA_INDEX_COLUMN(BC, bc);        //! Most probably BC to where this collision has occured
DECLARE_SOA_COLUMN(PosX, posX, float);   //! X Vertex position in cm
DECLARE_SOA_COLUMN(PosY, posY, float);   //! Y Vertex position in cm
DECLARE_SOA_COLUMN(PosZ, posZ, float);   //! Z Vertex position in cm
DECLARE_SOA_COLUMN(CovXX, covXX, float); //! Vertex covariance matrix
DECLARE_SOA_COLUMN(CovXY, covXY, float); //! Vertex covariance matrix
DECLARE_SOA_COLUMN(CovXZ, covXZ, float); //! Vertex covariance matrix
DECLARE_SOA_COLUMN(CovYY, covYY, float); //! Vertex covariance matrix
DECLARE_SOA_COLUMN(CovYZ, covYZ, float); //! Vertex covariance matrix
DECLARE_SOA_COLUMN(CovZZ, covZZ, float);
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);
DECLARE_SOA_COLUMN(IsInelGt0, isInelGt0, bool);
DECLARE_SOA_COLUMN(IsInelGt1, isInelGt1, bool);               //! Vertex covariance matrix
DECLARE_SOA_COLUMN(Flags, flags, uint16_t);                    //! Run 2: see CollisionFlagsRun2 | Run 3: see Vertex::Flags
DECLARE_SOA_COLUMN(Chi2, chi2, float);                         //! Chi2 of vertex fit
DECLARE_SOA_COLUMN(NumContrib, numContrib, uint16_t);          //! Number of tracks used for the vertex
DECLARE_SOA_COLUMN(CollisionTime, collisionTime, float);       //! Collision time in ns relative to BC stored in bc()
DECLARE_SOA_COLUMN(CollisionTimeRes, collisionTimeRes, float); //! Resolution of collision time
} // namespace drvcollision

DECLARE_SOA_TABLE(DrvCollisions, "AOD", "DRVCOLLISION", //! Time and vertex information of collision
                  o2::soa::Index<>, drvcollision::BCId,
                  drvcollision::PosX, drvcollision::PosY, drvcollision::PosZ,
                  drvcollision::CovXX, drvcollision::CovXY, collision::CovXZ, drvcollision::CovYY, drvcollision::CovYZ, drvcollision::CovZZ,
                  drvcollision::CentFT0C, drvcollision::CentFT0M, drvcollision::IsInelGt0, drvcollision::IsInelGt1, drvcollision::Flags, drvcollision::Chi2, drvcollision::NumContrib,
                  drvcollision::CollisionTime, drvcollision::CollisionTimeRes);
} // namespace o2::aod
