program sphere
   use omg_ui, only: OMG,                     & ! mesh generator
                     nDomains,                & ! )
                     ClipBoxLow, ClipBoxHigh, & ! ) input variables of OMG
                     GiantSize, nRefinements, & ! )
                     PreferredAxis
   nDomains = 1
   ClipBoxLow  = [0d0, -1d0, -1d0]  ! low limit of x,y,z in the Clip Box
   ClipBoxHigh = [0.9d0,  1d0, 1d0]    ! high limit of x,y,z in the Clip Box
   GiantSize = 0.5                   ! largest cell size in the mesh
   nRefinements = 5                  ! number of mesh refinement passes
   PreferredAxis = 'x'
   call OMG(DomainIndicator, RefineCell, SurfaceName)

contains

   integer function DomainIndicator(xyz); real(8),intent(in):: xyz(3)

      if (xyz(1)**2 + xyz(2)**2 + xyz(3)**2 > 1) then
          DomainIndicator = 0
      else
          DomainIndicator = 1       ! xyz is in Inner Domain
      end if
   end function DomainIndicator

   logical function RefineCell()
      use omg_ui,only: cell_onSIC
      RefineCell = cell_onSIC()
   end function RefineCell

   character(16) function SurfaceName()
      SurfaceName='*' ! i.e. use a default name
   end function SurfaceName
end program sphere
