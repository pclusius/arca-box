Variables =

ifdef fdx
	Variables += -Dfdx=$(fdx)
endif
ifdef nVerticalLayers
	Variables += -DfnVerticalLayers=$(nVerticalLayers)
endif
ifdef nHorisontalColumns
	Variables += -DfnHorisontalColumns=$(nHorisontalColumns)
endif
ifdef zIndexPlume
	Variables += -DfzIndexPlume=$(zIndexPlume)
endif
ifdef yIndexPlume
	Variables += -DfyIndexPlume=$(yIndexPlume)
endif
ifdef mixingTimestep
	Variables += -DfmixingTimestep=$(mixingTimestep)
endif
ifdef totalRuntime
	Variables += -DftotalRuntime=$(totalRuntime)
endif
ifdef nSpec
	Variables += -DfnSpec=$(nSpec)
endif
ifdef Path
	Variables += -DfPath=\"$(Path)\"
endif
ifdef use_KyKz
	Variables += -Dfuse_KyKz=$(KyKz)
endif
ifdef use_randomize_k
	Variables += -Dfuse_randomize_k=$(randomize_k)
endif
ifdef use_friction_vel
	Variables += -Dfuse_friction_vel=$(friction_vel)
endif
ifdef use_pblh
	Variables += -Dfuse_pblh=$(pblh)
endif
ifdef use_namelist
	Variables += -Dfuse_namelist=1
endif


all:
	gfortran -O2 -cpp $(Variables) seq_file_names.f90 variables.f90 gtritest.f90 -o mixer2d.exe

variables:
	@echo $(Variables)
