!***************************************************************************
!***************  MODUL ***************************************************
!***************************************************************************

module aleas
!del alea3d
integer, public                ::nvalo,nvaloq  !en el fons tenim nvalo**2
real*8 , public, allocatable   ::particions_esfera(:,:) !s'hauria de fer nomes de out
integer, public, parameter     ::tamstocas=100000   !es per fer la matriu de ran2
integer, public                ::nunualea,cotran !numero de numeros aleas k realment tenim
real*8 , public                ::stocas(tamstocas)

contains

!**************************************************************************
subroutine inialea3d(un)
  integer i,j,k
  real*8 a,b,c,e,x,y,z
  integer un
  nvalo=un
  nvaloq=un*un
  allocate(particions_esfera(nvaloq,3))
  e=2d0*3.141592/nvalo
  do i=1,nvalo
    do j=1,nvalo
      a=dcos(i*e)
      c=sqrt(1d0-a**2)
      x=c*dcos(j*e)
      y=c*dsin(j*e)
      z=a
      k=(i-1)*nvalo+j
      particions_esfera(k,1)=x
      particions_esfera(k,2)=y
      particions_esfera(k,3)=z
    end do
  end do
end subroutine inialea3d

!**************************************************************************
subroutine llaleat
    integer ui,i,ii,ios
    open(10,file='alea.dat',status='old',iostat=ios)
    if (ios/=0) then
      print *,"";      print *,"";      print *,"";
      print *,"there is no alea.dat file with aleatory numbers"
      print *,"";      print *,"";      print *,"";
      print *,"I can't tolerate that: f**ck you I'm quitting"
      print *,"";      print *,"";      print *,"";
      stop 
    end if
    nunualea=0
    cotran=1
    do ui=1,tamstocas
      read (10,*,ERR=100,END=90) stocas(ui)
      nunualea=nunualea+1
    end do
    return
100 print *,"";      print *,"";      print *,"";
    print *,"error in the random file alea.dat"
    print *,"";      print *,"";      print *,"";
    print *," I can't tolerate that: f**ck you I'm quitting"
    print *,"";      print *,"";      print *,"";
    stop 
90  close(10)
    open(10,file='alea.dat')
	print*,"end alea.dat"
	rewind(10)
	ii=int(a*1000)+1
print *,a,ii,"rew"
	do i=1,ii-1
          read(10,*)
	end do
end subroutine

end module
