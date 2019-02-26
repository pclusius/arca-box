!!!! module stores all constants here

module Constants

USE second_Precision

implicit none

REAL(dp),parameter ::   &
      Boltzmann_cons    = 5.67051e-8,                &   ! Boltzman constant (W/m2/K4)
      Stefan_boltz_con  = 1.3807d-23,                &   ! Stefan Boltzman constant (J/K)
      Karm              = 0.4,                       &   ! Von Karmans constant [KV]
      Omega             = 7.292E-5,                  &   ! Earth rotation speed (radians/second) [OMEGA1] 
      PI                = 3.141592654,               &   ! No comments
      Grav              = 9.80665,                   &   ! Gravity constant (m2/s^2)
      Rgas              = 8.31451E3,                 &   ! Universal gas constant (J/kmol/K)
      RGSCNS            = 287.05,                    &   ! Univ. gas constant of air/mean molecular weight (J/kg/K)
      Md                = 28.966,                    &   ! Mean molecular mass of air (kg/kmol)
      Avog              = 6.0221E23,                 &   ! Avogadro-number [molecules/mol]
      Clight            = 2.9979E8,                  &   ! Speed of light [m/s]                  
      Hplanck           = 6.6236E-34                     ! Placks Constant [J s]



  

end module
