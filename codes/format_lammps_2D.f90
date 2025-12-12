Program formatlammps

  implicit none 
  integer::i,j,k,N,alp,nalp,pi,flag,dum,Npt
  real(kind=8)::box,hbox,ibox,tmp,sigma_avg,rho
  real(kind=8),allocatable::sigma(:),x(:),y(:),z(:)
  integer,allocatable::pt(:)
  CHARACTER*500:: Buf,filetag
  IF(IARGC() < 3) THEN
     WRITE(6,*)"Incorrect syntax: should be 8 arguments :: N-Npt-rho"
     CALL EXIT(8)
  ENDIF

  CALL GETARG(1,Buf) 
  READ(Buf,*)N
  CALL GETARG(2,Buf) 
  READ(Buf,*)Npt
  CALL GETARG(3,Buf) 
  READ(Buf,*)rho

  if(rho>1.0d0)then
   write(*,*)"ARE YOU SURE, PLEASE CHECK THE DENSITY VALUE FOR 2D SYSTEM" 
   stop
  end if
  allocate(sigma(Npt),pt(N),x(N),y(N),z(N))
  sigma_avg = 0
  box = (real(N)/rho)**(1.0d0/2.0d0)
  open(unit=100,file='sample-file',action='read')
  read(100,*)
  read(100,*)
  do i =1,N
     read(100,*)pt(i),x(i),y(i),tmp
     sigma(pt(i)) = tmp 
     sigma_avg = sigma_avg + tmp
     z(i) =0.0d0
  end do
  open(unit=111,file='lammps-file')
  open(unit=222,file='para-file')
  sigma_avg = sigma_avg/real(N)
  box  = box*sigma_avg
  ibox = 1.0d0/box ; hbox = 0.50d0*box

  write(111,*)"LAMMPS data file via write_data, version 15 May 2015, timestep = 0"
  write(111,*)" "
  write(111,*)N," atoms"
  write(111,*)Npt," atom types"
  write(111,*)""
  write(111,'(F20.16,4X,F20.16,4X,10A)')-0.50d0*box, 0.50d0*box," xlo xhi"
  write(111,'(F20.16,4X,F20.16,4X,10A)')-0.50d0*box, 0.50d0*box," ylo yhi"
  write(111,'(F20.16,4X,F20.16,4X,10A)')-0.50d0, 0.50d0, " zlo zhi"
  write(111,*)"0.0000000000000000 0.0000000000000000 0.0000000000000000 xy xz yz"
  write(111,*)" "
  write(111,*)"Masses"
  write(111,*)" "
  do i = 1, Npt
     write(111,*)i,"1.0"
  end do
  write(111,*)" "
  write(111,*)"Atoms # atomic"
  write(111,*)" "

  do i =1,N
     x(i) =x(i) - box*anint(x(i)*ibox)
     y(i) =y(i) - box*anint(y(i)*ibox)
     z(i) = 0.0d0
     write(111,'(I5,4X,I5,4X,F20.16,4X,F20.16,4X,F20.16,4X,3(4A))')i,pt(i),x(i),y(i),z(i)," 0 ", " 0 ", " 0 "  
  end do
  write(111,*) "" 
  write(111,*) "Velocities" 
  write(111,*) "" 
  do i =1,N
     write(111,*)i," 0 ", " 0 ", " 0 " 
  end do
  j = 1
  do i =1,Npt
     write(222,'(A15, 4X, I5, 4X, I5, 4X, I5, 4X, F20.16, 4X, F20.16)')"pair_coeff",i,i,j,sigma(i),sigma(i)*1.250d0
  end do
  write(222,*)" "
  do i =1,Npt
     write(222,*)"mass",i,"1.0"
  end do
  write(222,*)" "
end program formatlammps

