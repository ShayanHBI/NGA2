program createDiffusionChemTable
   use param,                only: param_init,param_final,param_read
   use flameletlib_class,    only: flameletLib,sfm
   use diffusionTable_class, only: diffusionTable
   use messager,             only: die
   implicit none

   ! Declarations
   type(flameletLib), target :: flmlib
   type(diffusionTable)      :: dfftbl
   integer :: n

   ! Get the input file
   call param_init

   ! Construct the flameletLib object
   flmlib=flameletLib(model=sfm)

   ! Construct the diffusionTable object
   dfftbl=diffusionTable(flmlib=flmlib)

   ! Loop over files
   print*,''
   print*,'** Files in the table **'
   do n=1,dfftbl%flmlib%nfiles
      write(*,'(a)') trim(dfftbl%flmlib%files(n))
      ! Read the file
      call dfftbl%flmlib%readfile(n)
      ! Convolute with PDF
      call dfftbl%convolute(n)
      ! Deallocate the data array
      call dfftbl%flmlib%cleanup
   end do

   ! Change the variables names
   call dfftbl%convert_names
   ! Compute the table
   call dfftbl%setup
   ! Print some statistics
   call dfftbl%stats
   ! Write the table
   call dfftbl%write
   ! Clean up user parameters
   call param_final
end program createDiffusionChemTable
