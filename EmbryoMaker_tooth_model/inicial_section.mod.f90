!    EmbryoMaker software (General Node Model)
!    Computational model to simulate morphogenetic processes in living organs and tissues.
!    Copyright (C) 2014 Miquel Marin-Riera, Miguel Brun-Usan, Roland Zimm, Tommi Välikangas & Isaac Salazar-Ciudad

!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.




!***************************************************************************
!***************  MODUL INICIAL ********************************************
!***************************************************************************
module inicial
use general
use genetic    !>>>>>>>>>>>> Is 29-4-13
use aleas
use neighboring
use ic
use io
use nexus !>>Miquel24-9-14

character*8, public :: cdate

contains

!**************************************************************************
subroutine initials
character*200 kk
  itvi=1
  itviactual=1
  call getarg(1,carg)
  if (len_trim(carg)==0.or.carg=="0") then
    !call default_param_values
    call default_ic
    erep=0.0d0
    erepcel=0.0d0
    eyou=0.0d0
    eadh=0.0d0
    etor=0.0d0
    espring=0.0d0
    !<<< Is 14-3-15
    if (allocated(nodeo)) deallocate(nodeo)
    allocate(nodeo(nda))
    !print *,nda,"nda"
    nodeo=node  
    !<<< Is 14-3-15
  else
    call iniread
    call readsnap(carg)
    if(nd>1) call neighbor_build  !>>Miquel24-9-14
    !call nexe  !>>Miquel24-9-14
    if (errorlec==1.and.eva==0) stop  !>>> Is 22-1-14
  end if
  call iniio              ! this is just to allocate and inicializate the matrices for the variable names and stuff
  call inialea3d(nparti)  ! this is to inicialize the partition of random numbers in a sphere
  !call llaleat            ! this is to read a bunch of random numbers
  if(ffu(13)==0)then
    call iniboxes           ! this to inicialize the boxes
  end if
  if(nd>1) call neighbor_build
  call put_param_to_matrix(param)
  paramo=param


  call sectionate


  if(aut>0)  call writesnapini    !here we write the output file with the initial conditions              !>>>>Miquel18-11-13

  !call OPC

end subroutine

!**************************************************************************

subroutine default_ic

   call default_values          !IS 2-1-14 this is just some default values that get overrided by the ic
                                !but I put them here in case you forget
!  call epi_sphere              !OK

!  call mes_shell               !doesn't crash, but doesn't work
!  call mes_ecm                 !WTF it's not even mesenchyme, it's epithelium, DELETE?


!  call grn_test                !doesn't crash, but doesn't seem to work

!print*,"aut",aut
  !if(aut/=1 .and. aut/=5)then
  !  !IC for the mechanisms
  !  print*,"ATENTION,*************"
  !  print*,"you didn't load any input file,"
  !  print*,"1 Cell death mechanism"
  !  print*,"2 Differential adhesion mechanism"
  !  print*,"3 Cell contraction mechanism"
  !  print*,"4 Polarized cell growth and division mechanism"
  !  print*,"5 Extracellular Matrix secretion"
  !  print*,"6 Cell migration"
  !  read(*,*) ic_load
  !end if
!print*,"aut",aut,"icload",ic_load

!ic_load=3

  select case(ic_load)
    !case(1);  call epi_apoptosis            !OK ; CHECKED 25-8-14
    !case(1);  call scale            !OK ; CHECKED 25-8-14
    case(2);  call epidermal_organ_reduc         !OK ; CHECKED 25-8-14
    !case(2);  call torque_test         !OK ; CHECKED 25-8-14
    case(3);  call epidermal_organ_in_tooth_bud            !OK ; CHECKED 25-8-14
    case(4);  call epidermal_organ_in_tooth_bud_s2_small            !OK ; CHECKED 25-8-14
    case(5);  call epidermal_organ_in_tooth_bud_overall            !OK ; CHECKED 25-8-14
    case(6);  call hair_placode
    case(7);  call epi_folding_growth_test
    case(8);  call dros_wing_vein_2d
    !case(3);  call torque_test           !OK ; CHECKED 25-8-14
    !case(4);  call epi_polar_growth         !OK ; CHECKED 25-8-14
    !case(5);  call epi_mes_ecm         !OK ; CHECKED 25-8-14
    !case(6);  call migration              !OK ; CHECKED 25-8-14
    !case(7);  call mes_polar_growth                !OK ; CHECKED 25-8-14
    !case(8);  call blastula_ensemble 
    !case(9);  call father
  end select


!COOL SUBROUTINES

!  call polarized                !OK ; clean ; no ; CHECKED 2-7-14
!  call epi_growth_and_division(2,2)  !OK ; clean ; on ; CHECKED 7-7-14
! PROBLEM in mitosis call epi_mes_growth_and_division  !OK ; clean ; no ; CHECKED 2-7-14
!  call epi_mes                  !OK ; clean ; no ; CHECKED 2-7-14
!  call mes_ecm                  !OK ; clean ; si ; CHECKED 7-7-14
!  call neg_adh
!  call differential_growth_epi  !OK ; clean ; si ; CHECKED 7-7-14
!  call diffusion_test           !OK ; clean ; no
!  call mesenchyme               !OK ; clean ; no
!  call mesenchyme_apop          !OK ; clean ; si ; CHECKED 7-7-14
!  call teloblast                 !                ; CHECKED 1-7-14 NO !growth does weird things, adds nodes in a single row and outside the cell
!  call differential_growth_mes  !OK ; clean ; si ; CHECKED 4-7-14
!  call mes_polar_growth         !OK ; clean ; si  ; CHECKED 4-7-14
!  call directed_mit_mes         !OK ; clean ; si ; CHECKED 7-7-14 !kinda works,supposed to do assym. div. , but the gradient it's not strong enough I guess... 
!  call mes_shell                ! miguel4-11-13
!  call ic_emt                   !OK ; clean ;si  ; CHECKED 4-7-14
!  call twoeps
!  call diffusion_test
!  call epi_active_transport
!  call mes_ecm_degradation
!  call epi_mes_primordium
!  call mes_primordium
!  call epi_mes_bud
!  call hair_placode
!  call feather_placode
!  call epi_mes_bud_ingrowth
!  call tooth_bud
!  call delta_notch


!call blastuloid

!call founding_father

!call epi_growth_and_division

!call invaginacio_diff

!call blastula

!call blastula_ensemble

!call father

end subroutine default_ic


subroutine sectionate

integer::i,j,k

  m_ywall=2.0d0
  mi_ywall=-2.0d0

  i=1
  do while(i<=nd)
    if(gex(i,3)>0.1)then
        k=node(i)%icel
        call apoptosis(i)
        call eliminate_cell(k)
    else if(node(i)%tipus<3)then
      j=node(i)%altre
      if((node(i)%y<mi_ywall .or. node(i)%y>m_ywall).or.(node(j)%y<mi_ywall .or. node(j)%y>m_ywall))then
        k=node(i)%icel
        call apoptosis(i) !; call apoptosis(j)
        call eliminate_cell(k)
      else
        i=i+1
      end if
    else
      if((node(i)%y<mi_ywall .or. node(i)%y>m_ywall))then
        k=node(i)%icel
        call apoptosis(i)
        call eliminate_cell(k)
      else
        i=i+1
      end if
    end if
  end do

  do i=1,nd
    !print*,"i",i,"tipus",node(i)%tipus,"node altre",node(i)%altre
    nodeo(i)%altre=node(i)%altre
    nodeo(i)%tor=node(i)%tor
    nodeo(i)%stor=node(i)%stor
    nodeo(i)%kplast=node(i)%kplast
    nodeo(i)%kvol=node(i)%kvol
  end do

  !LETS REMOVE THE MESENCHYMAL HOLDS IN ONE SIDE

!  do i=1,nd
!    if(node(i)%tipus==3)then
!      if(gex(i,3)>0.1)then
!        if(node(i)%hold>0 .and. node(i)%x>0)then
!          node(i)%hold=0
!        end if
!      end if
!    end if
!  end do


  !LETS SUPRESS EPITHELIAL MESENCHYMAL ADHESION

  kadh(1,3)=0d0 ; kadh(3,1)=kadh(1,3)

  !we turn off growth

  do i=1,ng  
    gen(i)%wa(nparam_per_node+2)=0d0
  end do
  call update_npag


end subroutine sectionate

!si si si
end module inicial

