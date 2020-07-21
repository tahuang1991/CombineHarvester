#ifndef TAMUHHbbWW_HttSystematics_MSSMRun2_h
#define TAMUHHbbWW_HttSystematics_MSSMRun2_h
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"

namespace ch {
// Run2 MSSM analysis systematics
// Implemented in src/HttSystematics_MSSMRun2.cc
void AddMSSMRun2Systematics(CombineHarvester& cb, int control_region);
}

#endif
