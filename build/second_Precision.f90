!!!! module to define second precision

module second_precision

   integer, parameter :: dp = selected_real_kind (p=14,r=30)
   INTEGER, PARAMETER :: ik = SELECTED_INT_KIND(24)

end module second_precision
