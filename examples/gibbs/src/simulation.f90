!> Example to read in YAML thermodynamic file
module simulation
   use precision, only: WP
   use string,    only: str_medium
   use param, only: param_read
   
   ! YAML-specific
   use, intrinsic :: iso_fortran_env, only:  output_unit
   use fortran_yaml_c, only: YamlFile

   implicit none
   private
   
   public :: simulation_init,simulation_run,simulation_final
   
contains

   !> Nothing to initialize
   subroutine simulation_init
      implicit none  

      type(YamlFile) :: file
      character(:), allocatable :: err

      call file%parse("test1.yaml", err)
      if (allocated(err)) then
         print*,err
         stop 1
      endif
    
      call file%dump(unit=output_unit, indent=0)

   end subroutine simulation_init


   !> NOthing here yet
   subroutine simulation_run
      implicit none

      ! Nothing to do
      
   end subroutine simulation_run
   

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      
   end subroutine simulation_final


   
   
end module simulation
