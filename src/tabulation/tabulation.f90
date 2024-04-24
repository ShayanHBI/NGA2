!> Tabulation tools for chemistry:
module tabulation
   use string,   only: str_medium
   use param,    only: param_final,param_read
   use messager, only: die
   implicit none
   private

   ! Make public what is useful outside
   public :: tabulate_flamelet


   contains

   
   !> Creates a chemistry table and stores it in a binary file
   subroutine tabulate_flamelet(model,chfname)
      use, intrinsic :: iso_fortran_env, only: output_unit
      use flameletlib_class,             only: flameletLib
      use diffusionTable_class,          only: diffusionTable
      implicit none
      integer, intent(in) :: model
      character(len=str_medium), intent(in) :: chfname
      type(flameletLib), target :: flmlib
      type(diffusionTable) :: dfftbl
      integer :: ifile

      ! Log
      write(output_unit,'("*** Creating flamelet library ***")')

      ! Construct the flameletLib object
      flmlib=flameletLib(model=model)

      ! Construct the diffusionTable object
      dfftbl=diffusionTable(flmlib=flmlib,filename=chfname)

      ! Loop over files
      write(output_unit,'(" ")')
      write(output_unit,'("** Files in the table **")')
      do ifile=1,dfftbl%flmlib%nfiles
         write(*,'(a)') trim(dfftbl%flmlib%files(ifile))
         ! Read the file
         call dfftbl%flmlib%readfile(ifile)
         ! Convolute with PDF
         call dfftbl%convolute(ifile)
         ! Deallocate the data array
         call dfftbl%flmlib%cleanup()
      end do

      ! Change the variables names
      call dfftbl%convert_names()
      ! Compute the table
      call dfftbl%setup()
      ! Print some statistics
      call dfftbl%stats()
      ! Write the table
      call dfftbl%write()
   end subroutine tabulate_flamelet


end module tabulation
