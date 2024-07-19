#ifndef detector_construction_HH
#define detector_construction_HH 1 

#include <unordered_map>
#include <string>
#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4UserLimits;
class G4Material;
class G4GlobalMagFieldMessenger;

class detector_construction : public G4VUserDetectorConstruction{

  typedef detector_construction self_t;
  typedef G4VUserDetectorConstruction base_t;

public:
  detector_construction() = default;
  virtual ~detector_construction() = default;

  virtual G4VPhysicalVolume* Construct() override;
  virtual void ConstructSDandField() override;

private:
  double m_world_size;
  G4UserLimits* m_step_limits;
  std::unordered_map<std::string,G4Material*> m_materials;

  void init_materials();

  static G4ThreadLocal G4GlobalMagFieldMessenger* m_s_magfield_messenger;
};
#endif
