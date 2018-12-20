!!!!!! Module contains all the inputs used in "Superbox-model"

module Inputs


!!!!!!!!! Input flag subroutine. Contains flags for all the processes we would like to include
subroutine input_flags(Aerosol_flag, emission_flag, nucleation_flag, condensation_flag,coagulation_flag, &
                       Wall_loss_flag, chemistry_flag)
 

 ! Value for flags 0 = included, 1 = not included
 !!!!!!!!
 
 implicit none
 
 integer :: Aerosol_flag, emission_flag, nucleation_flag, condensation_flag,coagulation_flag, &
            Wall_loss_flag, chemistry_flag


end subroutine input_flags

end module
