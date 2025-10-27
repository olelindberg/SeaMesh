#ifndef PHYSICAL_FIELDS_CREATOR_H
#define PHYSICAL_FIELDS_CREATOR_H

#include "logger.h"

class PhysicalFieldsCreator
{
public:
  PhysicalFieldsCreator()
  { // Creating all prognostic fields:
    Logger::info("Creating prognostic fields ...");

    // Create surface elevation:
    Logger::info("Surface elevation");

    // Create horizontal velocity field:
    Logger::info("Horizontal velocity");

    // Create salinity:
    Logger::info("Salinity");

    // Create temperature:
    Logger::info("Temperature");

    // Create diagnostic fields:
    Logger::info("Creating diagnostic fields ...");

    // Create vertical velocity:
    Logger::info("Vertical velocity");

    // Surface vind stress:
    Logger::info("Surface wind stress");

    // Surface pressure:
    Logger::info("Surface pressure");

    // Create coriolis field:
    Logger::info("Coriolis field");

    // Create density field:
    Logger::info("Density field");

    // Create turbulent kininetic energy:
    Logger::info("Turbulent kinetic energy");

    // Create turbulent dissipation rate:
    Logger::info("Turbulent dissipation rate");
  }
  ~PhysicalFieldsCreator() {}
};

#endif // PHYSICAL_FIELDS_CREATOR_H