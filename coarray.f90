program coarray
   implicit none

!***********************************************!
!						!
!               Declaration of variables	!
!						!
!***********************************************!

   integer :: i, j, k, m, l, img 
! i,j,k are indices used in loops
! m,l are the array boundariesarray boundaries
! img is the image ID 
!-----------------------------------------------!

   integer, allocatable, target :: iarray(:,:)[:,:], jarray(:,:)[:,:] 
! 2D arrays with rank irank in the space of
! images.
!-----------------------------------------------!

   integer, allocatable :: larray(:)[:,:] 
! 1D array with rank irank in the space of 
! images.
!-----------------------------------------------!
 
   integer :: irank 
! Rank of the image array, this could be and
! input data. 
!-----------------------------------------------!

! Future...
!   integer, pointer :: e => null(), n => null()
!   integer, pointer :: w => null(), s => null()
!   integer, pointer :: ne => null(), nw => null()
!   integer, pointer :: sw => null(), se => null()
!-----------------------------------------------!

   integer :: sub(2) 
! coordinates of the sub-array
!-----------------------------------------------!

   character(len=50) :: format = '(6(i2," "),"  ",6(i2," "))'
! Format variable for output, needs to be 
! combined with the input of array dimension
!**********************************************

   m = 6 ! move this to imput data in the future
   l = 6 ! move this to imput data in the future

   img = this_image()

   irank = 2 ! move this to imput data in the future
   allocate(iarray(m,l)[irank,*])
   allocate(jarray(m,l)[irank,*])
   allocate(larray(m*l)[irank,*])

   larray = img*10
   sub = [this_image(iarray,1), this_image(iarray,2)]
!   write(6,'("Image ",i2," holds the subarray ",2i2)') img, sub

   print *, " "

  !  iarray = reshape(larray, [m, l])

  !  if (img .eq. 1) then
  !    write(6,*)"After larray:"
  !    do j=1,l
  !      write(unit=6, fmt=format) iarray(:,j)[1,1], iarray(:,j)[2,1]
  !    enddo
  !    print *, " "
  !    do j=1,l
  !      write(unit=6, fmt=format) iarray(:,j)[1,2], iarray(:,j)[2,2]
  !    enddo
  !    print *, " "
  !    print *, " "
  !  endif
  !  sync all

   iarray(2:m-1,1) = img*10+2
   iarray(2:m-1,l) = img*10+4
   iarray(m,2:l-1) = img*10+1
   iarray(1,2:l-1) = img*10+3
   iarray(1,1)     = img*10+6
   iarray(m,1)     = img*10+5
   iarray(1,l)     = img*10+7
   iarray(m,l) = (img-1)*10+8

   if (img .eq. 1) then
     write(6,*)"The initial setup:"
     do j=1,l
       write(unit=6, fmt=format) iarray(:,j)[1,1], iarray(:,j)[2,1]
     enddo
     print *, " "
     do j=1,l
       write(unit=6, fmt=format) iarray(:,j)[1,2], iarray(:,j)[2,2]
     enddo
     print *, " "
     print *, " "
   endif


   jarray = cshift(iarray, -1, 1)
   sync all

   if (img .eq. 1) then
     write(6,*)"After streaming to the east:"
     do j=1,l
       write(unit=6, fmt=format) jarray(:,j)[1,1], jarray(:,j)[2,1]
     enddo
     print *, " "
     do j=1,l
       write(unit=6, fmt=format) jarray(:,j)[1,2], jarray(:,j)[2,2]
     enddo
     print *, " "
     print *, " "
   endif
   sync all

   jarray = cshift(iarray, -1, 2)
   sync all

   if (img .eq. 1) then
     write(6,*)"After streaming to the north:"
     do j=1,l
       write(unit=6, fmt=format) jarray(:,j)[1,1], jarray(:,j)[2,1]
     enddo
     print *, " "
     do j=1,l
       write(unit=6, fmt=format) jarray(:,j)[1,2], jarray(:,j)[2,2]
     enddo
     print *, " "
     print *, " "
   endif
   sync all

   print *, " "
   print *, " "

   jarray = cshift(iarray, -1, 1)
   jarray = cshift(jarray,  1, 2)
   sync all

   if (img .eq. 1) then
     write(6,*)"After streaming to the north-east:"
     do j=1,l
       write(unit=6, fmt=format) jarray(:,j)[1,1], jarray(:,j)[2,1]
     enddo
     print *, " "
     do j=1,l
       write(unit=6, fmt=format) jarray(:,j)[1,2], jarray(:,j)[2,2]
     enddo
   endif

end program coarray
