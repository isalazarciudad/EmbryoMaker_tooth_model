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




! Is 20-5-13
module ic	!the mesenchimal and epithelial i.c. subroutines are now unified. 
use general
use genetic
use shell ! miguel4-11-13
use death
!use nexus

integer  :: nodecela
integer  :: l, ciclo, cont
real*8   :: x, xx, y, yy, z, zz, de, di, hip              
real*8   :: dx1, dx2, dx3, dy1, dy2, dy3
integer, public  :: radi,radicel,layer,mradi,mradicel,xradi,xlayer
real*8, public   ::  zepi,zmes,radius,zx
integer, public  :: packed    !>>>Miquel10-1-14
integer          :: iccag
integer,public,allocatable :: mesradicel(:)

contains

!*************************************************************************************

subroutine default_values  !IS 2-1-14

    ndmax=9000
    ecmmax=0.25d0
    realtime=0
    ttalone=10
    mnn=500
    dif_req=1.0d0
    min_comp=-1.0d-1
    screen_radius=1.0d0
    reqmin=0.05d0
    df_reqmax=0.50d0
    deltamin=1.0d-3
    prec=0.1  ! Is 26-8-14
    angletor=0.00 !Miquel15-9-14
    ldi=epsilod

end subroutine

!*************************************************************************************



!***********************************************************************************************************************


!***********************************************************************************************************************

subroutine epiteli(radi,radicel,zepi,dreq,dreqs)		!zepi is the z-position of the bottom layer of nodes
integer            ::radi,radicel,valcel
real*8,dimension(:)::p1(3),p2(3),vector(3)
real*8             ::beta,angle,alt,de,di,zepi,dreq,dreqs

    di=dreqs*2
    de=dreq*2  !miguel 14-10-13
    vector=0.0d0
    
    node(:)%marge=1
 

	beta=2d0*pi/6d0


	do i=1,ncelsepi
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do


	ii=0
	do i=1,2
		ii=ii+1
!		alt=real(i-1)+zepi
        if(i==1)then
          alt=zepi
        else
          alt=zepi+di
        end if
!		print*,"alt",alt
		node(ii)%x=0.0d0;node(ii)%y=0d0;node(ii)%z=alt
        if(i==1) node(ii)%marge=0
		if(i==1)then
			node(ii)%tipus=2
		else if(i==2)then
			node(ii)%tipus=1
		end if
		cels(1)%nunodes=cels(1)%nunodes+1
		cels(1)%node(cels(1)%nunodes)=ii
!                node(ii)%marge=1
		node(ii)%icel=1
		do j=2,radi
			angle=beta
!			print*,"angle",angle
			d=de*(j-1d0)
			ii=ii+1
			node(ii)%x=d;node(ii)%y=0;node(ii)%z=alt
			if(i==1)then
				node(ii)%tipus=2
			else
				node(ii)%tipus=1
			end if
			cels(1)%nunodes=cels(1)%nunodes+1
			cels(1)%node(cels(1)%nunodes)=ii
			node(ii)%icel=1
			p1(1)=d;p1(2)=0;p1(3)=alt
			p2(1)=d*dcos(angle);p2(2)=d*dsin(angle);p2(3)=alt
			jj=j-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
			if(jj>0)then
!				print*,"jj",jj
				vector=p2-p1
				modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
				iii=jj+1
				vector=vector/iii
				do kk=1,jj
					ii=ii+1
					node(ii)%x=p1(1)+vector(1)*kk;node(ii)%y=p1(2)+vector(2)*kk;node(ii)%z=alt
					if(i==1)then
						node(ii)%tipus=2
					else
						node(ii)%tipus=1
					end if
					cels(1)%nunodes=cels(1)%nunodes+1
					cels(1)%node(cels(1)%nunodes)=ii
					node(ii)%icel=1
				end do
			end if
			do k=2,5
				angle=angle+beta
				p1=p2
				ii=ii+1
				node(ii)%x=p1(1);node(ii)%y=p1(2);node(ii)%z=alt
				if(i==1)then
					node(ii)%tipus=2
				else
					node(ii)%tipus=1
				end if
				cels(1)%nunodes=cels(1)%nunodes+1
				cels(1)%node(cels(1)%nunodes)=ii
				node(ii)%icel=1
				p2(1)=d*dcos(angle);p2(2)=d*dsin(angle);p2(3)=alt
				jj=j-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
				if(jj>0)then
					vector=p2-p1
					modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
					iii=jj+1
					vector=vector/iii
					do kk=1,jj
						ii=ii+1
						node(ii)%x=p1(1)+vector(1)*kk;node(ii)%y=p1(2)+vector(2)*kk;node(ii)%z=alt
						if(i==1)then
							node(ii)%tipus=2
						else
							node(ii)%tipus=1
						end if
						cels(1)%nunodes=cels(1)%nunodes+1
						cels(1)%node(cels(1)%nunodes)=ii
						node(ii)%icel=1
					end do
				end if
			end do
			angle=angle+beta
!			print*,"angle",angle
			p1=p2
			ii=ii+1
			node(ii)%x=p1(1);node(ii)%y=p1(2);node(ii)%z=alt
			if(i==1)then
				node(ii)%tipus=2
			else
				node(ii)%tipus=1
			end if
			cels(1)%nunodes=cels(1)%nunodes+1
			cels(1)%node(cels(1)%nunodes)=ii
			node(ii)%icel=1
			p2(1)=d;p2(2)=0;p2(3)=alt
			jj=j-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
			if(jj>0)then
				vector=p2-p1
				modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
				iii=jj+1
				vector=vector/iii
				do kk=1,jj
					ii=ii+1
					node(ii)%x=p1(1)+vector(1)*kk;node(ii)%y=p1(2)+vector(2)*kk;node(ii)%z=alt
					if(i==1)then
						node(ii)%tipus=2
					else
						node(ii)%tipus=1
					end if
					cels(1)%nunodes=cels(1)%nunodes+1
					cels(1)%node(cels(1)%nunodes)=ii
					node(ii)%icel=1	
				end do
			end if
		end do
	end do

!hem fet la 1a cel al centre 

!ara farem la resta copiant i pegant



if(radicel>1)then
	a=2*sqrt(((radi-1)*de)**2-(0.5*de*(radi-1))**2)+de	!distancia de centre de cel·lula a centre de cel·lula
	valcel=1
	do i=2,radicel
!		angle=0.5*pi		!"ok" hexagonal configuration
		angle=0.5*pi+pi/12	!perfect hexagonal configuration
!		print*,"i",i,"angle",angle
!		vector(1)=0;vector(2)=a*(i-1);vector(3)=0			!"ok" hexagonal configuration
		vector(1)=-a*(i-1)*sin(pi/12);vector(2)=a*(i-1)*cos(pi/12);vector(3)=0	!perfect hexagonal configuration
		valcel=valcel+1
		do j=1,cels(1)%nunodes
			ii=ii+1
            if (j==1) node(ii)%marge=0
			node(ii)%x=node(j)%x+vector(1);node(ii)%y=node(j)%y+vector(2);node(ii)%z=node(j)%z
			node(ii)%tipus=node(j)%tipus
			cels(valcel)%nunodes=cels(valcel)%nunodes+1
			cels(valcel)%node(cels(valcel)%nunodes)=ii
			node(ii)%icel=valcel
		end do

		p1(1)=a*(i-1)*dcos(angle-beta);p1(2)=a*(i-1)*dsin(angle-beta);p2(3)=0
		p2(1)=a*(i-1)*dcos(angle);p2(2)=a*(i-1)*dsin(angle);p2(3)=0
		jj=i-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
		if(jj>0)then
			vector=p2-p1
			modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
			iii=jj+1
			vector=vector/iii
			do kk=1,jj
				valcel=valcel+1
				do j=1,cels(1)%nunodes
					kkk=cels(valcel-1)%node(j)
					ii=ii+1
                    if (j==1) node(ii)%marge=0
					node(ii)%x=node(kkk)%x-vector(1);node(ii)%y=node(kkk)%y-vector(2);node(ii)%z=node(kkk)%z
					node(ii)%tipus=node(j)%tipus
					cels(valcel)%nunodes=cels(valcel)%nunodes+1
					cels(valcel)%node(cels(valcel)%nunodes)=ii
					node(ii)%icel=valcel
				end do
			end do
		end if
		p1=p2
!if(1==2)then
		do k=2,5
			angle=angle+beta
!			print*,"angle",angle
			p2(1)=a*(i-1)*dcos(angle);p2(2)=a*(i-1)*dsin(angle);p2(3)=0
			valcel=valcel+1
			do j=1,cels(1)%nunodes
				ii=ii+1
                if (j==1) node(ii)%marge=0
				node(ii)%x=node(j)%x+p2(1);node(ii)%y=node(j)%y+p2(2);node(ii)%z=node(j)%z
				node(ii)%tipus=node(j)%tipus
				cels(valcel)%nunodes=cels(valcel)%nunodes+1
				cels(valcel)%node(cels(valcel)%nunodes)=ii
				node(ii)%icel=valcel
			end do
			jj=i-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
			if(jj>0)then
				vector=p2-p1
				modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
!				iii=jj+1
				vector=vector/(jj+1)
				do kk=1,jj
					valcel=valcel+1
					do j=1,cels(1)%nunodes
						kkk=cels(valcel-1)%node(j)
						ii=ii+1
                                                if (j==1) node(ii)%marge=0
						node(ii)%x=node(kkk)%x-vector(1);node(ii)%y=node(kkk)%y-vector(2);node(ii)%z=node(kkk)%z
						node(ii)%tipus=node(j)%tipus
						cels(valcel)%nunodes=cels(valcel)%nunodes+1
						cels(valcel)%node(cels(valcel)%nunodes)=ii
						node(ii)%icel=valcel
					end do
				end do
			end if
			p1=p2
		end do
		angle=angle+beta
!		p2(1)=0;p2(2)=a*(i-1);p2(3)=0
		p2(1)=a*(i-1)*dcos(angle);p2(2)=a*(i-1)*dsin(angle);p2(3)=0
		valcel=valcel+1
		do j=1,cels(1)%nunodes
			ii=ii+1
            if (j==1) node(ii)%marge=0
			node(ii)%x=node(j)%x+p2(1);node(ii)%y=node(j)%y+p2(2);node(ii)%z=node(j)%z
			node(ii)%tipus=node(j)%tipus
			cels(valcel)%nunodes=cels(valcel)%nunodes+1
			cels(valcel)%node(cels(valcel)%nunodes)=ii
			node(ii)%icel=valcel
		end do
		jj=i-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
		if(jj>0)then
			vector=p2-p1
			modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
!			iii=jj+1
			vector=vector/(jj+1)
			do kk=1,jj
				valcel=valcel+1
				do j=1,cels(1)%nunodes
					kkk=cels(valcel-1)%node(j)
					ii=ii+1
                                        if (j==1) node(ii)%marge=0
					node(ii)%x=node(kkk)%x-vector(1);node(ii)%y=node(kkk)%y-vector(2);node(ii)%z=node(kkk)%z
					node(ii)%tipus=node(j)%tipus
					cels(valcel)%nunodes=cels(valcel)%nunodes+1
					cels(valcel)%node(cels(valcel)%nunodes)=ii
					node(ii)%icel=valcel
				end do
			end do
		end if
!end if
	end do
end if




	!define the cell's centroid

	do i=1,ncelsepi
		cels(i)%ctipus=1
		cels(i)%cex=0;cels(i)%cey=0;cels(i)%cez=0
		do j=1,cels(i)%nunodes
			k=cels(i)%node(j)
            if(node(k)%tipus==1)then
	 			cels(i)%cex=cels(i)%cex+node(k)%x
 				cels(i)%cey=cels(i)%cey+node(k)%y
				cels(i)%cez=cels(i)%cez+node(k)%z
			end if
		end do
		cels(i)%cex=2*cels(i)%cex/real(cels(i)%nunodes)
		cels(i)%cey=2*cels(i)%cey/real(cels(i)%nunodes)
		cels(i)%cez=2*cels(i)%cez/real(cels(i)%nunodes)
	end do


	do i=1,ndepi		!veins i parametres mecanics
		ii=node(i)%icel
		if(node(i)%tipus==1)then
   		  node(i)%altre=i-cels(1)%nunodes/2
		else
		  node(i)%altre=i+cels(1)%nunodes/2
		end if
	end do

!        node(1)%marge=0.0d0

        !now we want the nucleus to be the first node in cels%node: this allows diffusion to be faster >>> Is 14-9-13 
	do i=1,ncelsepi
   	  do j=1,cels(i)%nunodes							
            k=cels(i)%node(j)								
            if (node(k)%marge==0) then
              ii=k
              exit
            end if
	  end do										
          jjj=cels(i)%node(1)
          cels(i)%node(1)=ii
   	  do j=1,cels(i)%nunodes							
            k=cels(i)%node(j)								
            if (k==ii) then
              jj=j
              exit
            end if
	  end do									
          cels(i)%node(jj)=jjj
	end do											


    !some funny rotation to make a less biased diffusion lost
    do i=1,nd
      c=0.75
      a=node(i)%x*cos(c)-node(i)%y*sin(c)
      b=node(i)%x*sin(c)+node(i)%y*cos(c)
      node(i)%x=a
      node(i)%y=b
    end do

    node(:)%talone=0.0d0

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111

subroutine mesenq(radi,radicel,layer,zmes,dreq)									!!! miguel 20.11
integer  ::i,j,k,ii,jj,radi,radicel,layer,signo,ent,ont       ! number of layers and concentric hexagons 
real*8   ::rad,der,zmes,de,di,dreq
real*8   :: xx,yy,zz                                                        ! miguel 4-6-13
rad=pi/180d0
de=dreq*2                    ! call radius
!di=2.0d0*de                 ! distance between cells and layers (it has to be >2 to avoid cell contacts)
di=2*de+2*de*cos(60d0*rad)

	do i=ncelsepi+1,ncels               
!		cels(i)%nunodes=radi+1
		cels(i)%nunodes=radi				!>>>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
!		allocate(cels(i)%node(radi+1))
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%node=0
    end do

node(ndepi+1:)%marge=1

!	print*,"node",size(node),"cels",size(cels),"celsnode",size(cels(1)%node)


	kkk=ndmes/layer
! radi=radi+1
	cont=ncelsepi+1

	ii=ndepi		!node counter
	jj=ncelsepi		!cel counter


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! central nodes (=nucleus) !!!!

   do l=1,layer  
!	jj=(kkk*(l-1))+radi+ndepi+1
	ii=ii+1
	jj=jj+1
!	node(ii)%x=(sqrt(di/2d0))*mod(l+1,2) ; node(ii)%y=(sqrt(di/2d0))*mod(l+1,2) ! Euler packing algorhytm
	node(ii)%x=0d0 ; node(ii)%y=0d0 ! Euler packing algorhytm
	node(ii)%z=zmes-di*(l-1)	!zmes marks the z-position of the uppermost layer
!	print*,"ii",ii,node(ii)%x,node(ii)%y,node(ii)%z
    node(ii)%tipus=3 ; node(ii)%icel=jj                        ! origin	
	cels(jj)%node(1)=ii
!	print*,"cel central, node central. ii:",ii
	ent=ii
	ont=ii
	!fill with the external nodes
	if(radi>1)then
        signo=1                                                              ! miguel 4-6-13
		do k=2,radi
			ii=ii+1
            !print*,"cel perif. ent:",ent,ii,k-1, 
             call random_number(a)
			der=de!(sqrt(de)*sqrt(de*a))                               ! miguel 4-6-13
121         call random_number(a)
            xx=der*(1d0-2*a) ; 
            call random_number(a)
            yy=der*(1d0-2*a)            ! miguel 4-6-13 
            zz=(der**2)-(xx**2)-(yy**2) 									 ! miguel 4-6-13
            if(zz.lt.0)then;goto 121;endif                                   ! miguel 4-6-13
            zz=signo*sqrt(zz) ; signo=-1*signo                               ! miguel 4-6-13  
            if(mod(k-1,3).eq.0)then      ; b=xx ; c=yy ; a=zz                ! miguel 4-6-13
            else if(mod(k-1,3).eq.1)then ; a=xx ; c=yy ; b=zz                ! miguel 4-6-13
  	        else if(mod(k-1,3).eq.2)then ; a=xx ; b=yy ; c=zz ; end if       ! miguel 4-6-13       
			!a=(2*ran2(idum)-1d0)*de;b=(2*ran2(idum)-1d0)*de;c=(2*ran2(idum)-1d0)*de   ! miguel 4-6-13     
			node(ii)%x=node(ent)%x+a
			node(ii)%y=node(ent)%y+b
			node(ii)%z=node(ent)%z+c
		    node(ii)%tipus=3 ; node(ii)%icel=jj
			cels(jj)%node(k)=ii
		end do
	end if

	do j=2,radicel                                                        ! "cicle"
		dx1=0.0 ; dx2=0.0 ; dx3=0.0 ; dy1=0.0 ; dy2=0.0 ; dy3=0.0 
		do i=1,6                                                        ! "sector" of the hexagon
			if(i.eq.1)then   

				dx2=0!*cos(real(i)*60d0*rad) 
				dy2=di*(real(j)-1d0)!*sin(real(i)*60d0*rad);print*,"dx2",dx2,"dy2",dy2
				dx1=di*(real(j)-1d0)*sin(-60d0*rad) ; dy1=di*(real(j)-1d0)*cos(-60d0*rad)!;print*,"dx1",dx1

			else
				hip=di*(real(j)-1d0)                                    ! hipotenusa
				dx1=dx2 ; dy1=dy2
				dx2=hip*sin(real(i-1)*60d0*rad) ; dy2=hip*cos(real(i-1)*60d0*rad)          
			end if            
			ii=ii+1
			jj=jj+1
!			print*,"cel perif, node central. ent:",ent,ii
			node(ii)%x=node(ont)%x+dx2 ; node(ii)%y=node(ont)%y+dy2 ; node(ii)%z=node(ont)%z
!			print*,"vertex i",i,"ii",ii,node(ii)%x,node(ii)%y,node(ii)%z
!			print*,i,"ii",node(ii)%x,"ent",node(ent)%x
			node(ii)%icel=jj ; node(ii)%tipus=3
			cels(jj)%node(1)=ii
			ent=ii
			if(radi>1)then
                signo=1                                                                  ! miguel 4-6-13
				do k=2,radi
					ii=ii+1
!                            print*,"cel perif, node extern. ent:",ent,ii
                            call random_number(a)
                            der=de!(sqrt(de)*sqrt(de*a))                               ! miguel 4-6-13
122                         call random_number(a)
                            xx=der*(1d0-2*a) ; 
                            call random_number(a)
                            yy=der*(1d0-2*a)            ! miguel 4-6-13 
                            zz=(der**2)-(xx**2)-(yy**2) 									 ! miguel 4-6-13
                            if(zz.lt.0)then;goto 122;endif                                   ! miguel 4-6-13
                            zz=signo*sqrt(zz) ; signo=-1*signo                               ! miguel 4-6-13  
                            if(mod(k-1,3).eq.0)then      ; b=xx ; c=yy ; a=zz               ! miguel 4-6-13
                            else if(mod(k-1,3).eq.1)then ; a=xx ; c=yy ; b=zz               ! miguel 4-6-13
  	                        else if(mod(k-1,3).eq.2)then ; a=xx ; b=yy ; c=zz ; end if      ! miguel 4-6-13

!					a=(2*ran2(idum)-1d0)*de;b=(2*ran2(idum)-1d0)*de;c=(2*ran2(idum)-1d0)*de    ! miguel 4-6-13     
					node(ii)%x=node(ent)%x+a
					node(ii)%y=node(ent)%y+b
					node(ii)%z=node(ent)%z+c
				    node(ii)%tipus=3 ; node(ii)%icel=jj
					cels(jj)%node(k)=ii
				end do
			end if

			if(j.gt.2)then                                              ! intermediate points
				dx3=dx2-dx1       ; dy3=dy2-dy1                         ! vectors which link "extreme" points
				dx3=dx3/real(j-1) ; dy3=dy3/real(j-1)                   ! sub-vectors                    
				do k=1,j-2                               
					ii=ii+1
					jj=jj+1
					node(ii)%x=node(ont)%x+(dx1+(real(k)*dx3)) 
					node(ii)%y=node(ont)%y+(dy1+(real(k)*dy3))    
					node(ii)%z=node(ont)%z                                                
!					print*,"intra i",i,"ii",ii,node(ii)%x,node(ii)%y,node(ii)%z
					node(ii)%icel=jj ; node(ii)%tipus=3 ; cels(jj)%node(1)=ii
					ent=ii
					if(radi>1)then
                        signo=1											            ! miguel 4-6-13
						do kk=2,radi
							ii=ii+1
!							print*,"cel perif, node extern. ent:",ent,ii
                            call random_number(a)
							der=de!(sqrt(de)*sqrt(de*a))                               ! miguel 4-6-13
123                         call random_number(a) 
                            xx=der*(1d0-2*a) ; 
                            call random_number(a)
                            yy=der*(1d0-2*a)            ! miguel 4-6-13 
                            zz=(der**2)-(xx**2)-(yy**2) 									 ! miguel 4-6-13
                            if(zz.lt.0)then;goto 123;endif                                   ! miguel 4-6-13
                            zz=signo*sqrt(zz) ; signo=-1*signo                               ! miguel 4-6-13  
                            if(mod(ii,3).eq.0)then      ; b=xx ; c=yy ; a=zz               ! miguel 4-6-13
                            else if(mod(ii,3).eq.1)then ; a=xx ; c=yy ; b=zz               ! miguel 4-6-13
  	                        else if(mod(ii,3).eq.2)then ; a=xx ; b=yy ; c=zz ; end if      ! miguel 4-6-13
                        
							!a=(2*ran2(idum)-1d0)*de;b=(2*ran2(idum)-1d0)*de;c=(2*ran2(idum)-1d0)*de    ! miguel 4-6-13     
							node(ii)%x=node(ent)%x+a
							node(ii)%y=node(ent)%y+b
							node(ii)%z=node(ent)%z+c
						    node(ii)%tipus=3 ; node(ii)%icel=jj
							cels(jj)%node(kk)=ii
						end do
					end if
				end do
			end if

		end do
	end do
  end do




	!define the cell's centroid and nucleus (margin)			!>>>>>>>>>>>>>> Miquel 14-4-13

	do i=ncelsepi+1,ncels											!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%ctipus=3									!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%cex=0;cels(i)%cey=0;cels(i)%cez=0			!>>>>>>>>>>>>>> Miquel 14-4-13
!		print*,"i",i,"node",cels(i)%node(:),"ncels",ncels
		do j=1,cels(i)%nunodes								!>>>>>>>>>>>>>> Miquel 14-4-13
			k=cels(i)%node(j)								!>>>>>>>>>>>>>> Miquel 14-4-13
!			print*,"k",k,"nunodes",cels(i)%nunodes
			cels(i)%cex=cels(i)%cex+node(k)%x				!>>>>>>>>>>>>>> Miquel 14-4-13
			cels(i)%cey=cels(i)%cey+node(k)%y				!>>>>>>>>>>>>>> Miquel 14-4-13
			cels(i)%cez=cels(i)%cez+node(k)%z				!>>>>>>>>>>>>>> Miquel 14-4-13
		end do												!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%cex=cels(i)%cex/real(cels(i)%nunodes)		!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%cey=cels(i)%cey/real(cels(i)%nunodes)		!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%cez=cels(i)%cez/real(cels(i)%nunodes)		!>>>>>>>>>>>>>> Miquel 14-4-13
	end do													!>>>>>>>>>>>>>> Miquel 14-4-13

	do i=ncelsepi+1,ncels
          b=1.0d8											!>>>>>>>>>>>>>> Is 14-9-13
   	  do j=1,cels(i)%nunodes								!>>>>>>>>>>>>>> Is 14-9-13
            k=cels(i)%node(j)								!>>>>>>>>>>>>>> Is 14-9-13
   	    a=sqrt((cels(i)%cex-node(k)%x)**2+(cels(i)%cey-node(k)%y)**2+(cels(i)%cez-node(k)%z)**2)
            if (b>a) then ; b=a ; ii=k ; jj=j ; end if
	  end do												!>>>>>>>>>>>>>> Is 14-9-13
          node(ii)%marge=0

          !now we want the nucleus to be the first node in cels%node: this allows diffusion to be faster  
          jjj=cels(i)%node(1)
          cels(i)%node(1)=ii
          cels(i)%node(jj)=jjj
	end do													!>>>>>>>>>>>>>> Is 14-9-13

    node(:)%talone=0.0d0

end subroutine mesenq

!**********************************************************************************************************

subroutine mesenq_packed(radi,radicel,layer,zepi,dreq)
integer            ::radi,radicel,valcel,layer,valcelo
real*8,dimension(:)::p1(3),p2(3),vector(3)
real*8             ::beta,angle,alt,de,di,zepi,dreq

    di=2d0*dreq
    de=2d0*dreq  !miguel 14-10-13
    vector=0.0d0
    
    node(ndepi+1:nd)%marge=1
 

	beta=2d0*pi/6d0


	do i=ncelsepi+1,ncels
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do


	ii=ndepi
	do i=1,radi
		ii=ii+1
		alt=real(i-1)*di+zepi
!        if(i==1)then
!          alt=zepi
!        else
!          alt=zepi+di
!        end if
!		print*,"alt",alt
		node(ii)%x=0.0d0;node(ii)%y=0d0;node(ii)%z=alt
        if(i==1) node(ii)%marge=0
        node(ii)%tipus=3
		cels(ncelsepi+1)%nunodes=cels(ncelsepi+1)%nunodes+1
		cels(ncelsepi+1)%node(cels(ncelsepi+1)%nunodes)=ii
		node(ii)%icel=ncelsepi+1
		do j=2,radi
			angle=beta
!			print*,"angle",angle
			d=de*(j-1d0)
			ii=ii+1
			node(ii)%x=d;node(ii)%y=0;node(ii)%z=alt
			node(ii)%tipus=3
			cels(ncelsepi+1)%nunodes=cels(ncelsepi+1)%nunodes+1
			cels(ncelsepi+1)%node(cels(ncelsepi+1)%nunodes)=ii
			node(ii)%icel=ncelsepi+1
			p1(1)=d;p1(2)=0;p1(3)=alt
			p2(1)=d*dcos(angle);p2(2)=d*dsin(angle);p2(3)=alt
			jj=j-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
			if(jj>0)then
!				print*,"jj",jj
				vector=p2-p1
				modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
				iii=jj+1
				vector=vector/iii
				do kk=1,jj
					ii=ii+1
					node(ii)%x=p1(1)+vector(1)*kk;node(ii)%y=p1(2)+vector(2)*kk;node(ii)%z=alt
					node(ii)%tipus=3
					cels(ncelsepi+1)%nunodes=cels(ncelsepi+1)%nunodes+1
					cels(ncelsepi+1)%node(cels(ncelsepi+1)%nunodes)=ii
					node(ii)%icel=ncelsepi+1
				end do
			end if
			do k=2,5
				angle=angle+beta
				p1=p2
				ii=ii+1
				node(ii)%x=p1(1);node(ii)%y=p1(2);node(ii)%z=alt
				node(ii)%tipus=3
				cels(ncelsepi+1)%nunodes=cels(ncelsepi+1)%nunodes+1
				cels(ncelsepi+1)%node(cels(ncelsepi+1)%nunodes)=ii
				node(ii)%icel=ncelsepi+1
				p2(1)=d*dcos(angle);p2(2)=d*dsin(angle);p2(3)=alt
				jj=j-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
				if(jj>0)then
					vector=p2-p1
					modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
					iii=jj+1
					vector=vector/iii
					do kk=1,jj
						ii=ii+1
						node(ii)%x=p1(1)+vector(1)*kk;node(ii)%y=p1(2)+vector(2)*kk;node(ii)%z=alt
						node(ii)%tipus=3
						cels(ncelsepi+1)%nunodes=cels(ncelsepi+1)%nunodes+1
						cels(ncelsepi+1)%node(cels(ncelsepi+1)%nunodes)=ii
						node(ii)%icel=1
					end do
				end if
			end do
			angle=angle+beta
!			print*,"angle",angle
			p1=p2
			ii=ii+1
			node(ii)%x=p1(1);node(ii)%y=p1(2);node(ii)%z=alt
			node(ii)%tipus=3
			cels(ncelsepi+1)%nunodes=cels(ncelsepi+1)%nunodes+1
			cels(ncelsepi+1)%node(cels(ncelsepi+1)%nunodes)=ii
			node(ii)%icel=ncelsepi+1
			p2(1)=d;p2(2)=0;p2(3)=alt
			jj=j-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
			if(jj>0)then
				vector=p2-p1
				modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
				iii=jj+1
				vector=vector/iii
				do kk=1,jj
					ii=ii+1
					node(ii)%x=p1(1)+vector(1)*kk;node(ii)%y=p1(2)+vector(2)*kk;node(ii)%z=alt
					node(ii)%tipus=3
					cels(ncelsepi+1)%nunodes=cels(ncelsepi+1)%nunodes+1
					cels(ncelsepi+1)%node(cels(ncelsepi+1)%nunodes)=ii
					node(ii)%icel=ncelsepi+1	
				end do
			end if
		end do
	end do

!hem fet la 1a cel al centre 

!ara farem la resta copiant i pegant


	valcel=ncelsepi+1  !>>>Miquel28-1-14
if(radicel>1)then
	a=2*sqrt(((radi-1)*de)**2-(0.5*de*(radi-1))**2)+de	!distancia de centre de cel·lula a centre de cel·lula
	do i=2,radicel
!		angle=0.5*pi		!"ok" hexagonal configuration
		angle=0.5*pi+pi/12	!perfect hexagonal configuration
!		print*,"i",i,"angle",angle
!		vector(1)=0;vector(2)=a*(i-1);vector(3)=0			!"ok" hexagonal configuration
		vector(1)=-a*(i-1)*sin(pi/12);vector(2)=a*(i-1)*cos(pi/12);vector(3)=0	!perfect hexagonal configuration
		valcel=valcel+1
		do j=1,cels(ncelsepi+1)%nunodes
			ii=ii+1
            if (node(ndepi+j)%marge==0) node(ii)%marge=0
			node(ii)%x=node(ndepi+j)%x+vector(1);node(ii)%y=node(ndepi+j)%y+vector(2);node(ii)%z=node(ndepi+j)%z
			node(ii)%tipus=node(ndepi+j)%tipus
			cels(valcel)%nunodes=cels(valcel)%nunodes+1
			cels(valcel)%node(cels(valcel)%nunodes)=ii
			node(ii)%icel=valcel
		end do

		p1(1)=a*(i-1)*dcos(angle-beta);p1(2)=a*(i-1)*dsin(angle-beta);p2(3)=0
		p2(1)=a*(i-1)*dcos(angle);p2(2)=a*(i-1)*dsin(angle);p2(3)=0
		jj=i-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
		if(jj>0)then
			vector=p2-p1
			modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
			iii=jj+1
			vector=vector/iii
			do kk=1,jj
				valcel=valcel+1
				do j=1,cels(ncelsepi+1)%nunodes
					kkk=cels(valcel-1)%node(j)
					ii=ii+1
                    if (node(ndepi+j)%marge==0) node(ii)%marge=0
					node(ii)%x=node(kkk)%x-vector(1);node(ii)%y=node(kkk)%y-vector(2);node(ii)%z=node(kkk)%z
					node(ii)%tipus=node(ndepi+j)%tipus
					cels(valcel)%nunodes=cels(valcel)%nunodes+1
					cels(valcel)%node(cels(valcel)%nunodes)=ii
					node(ii)%icel=valcel
				end do
			end do
		end if
		p1=p2
		do k=2,5
			angle=angle+beta
!			print*,"angle",angle
			p2(1)=a*(i-1)*dcos(angle);p2(2)=a*(i-1)*dsin(angle);p2(3)=0
			valcel=valcel+1
			do j=1,cels(ncelsepi+1)%nunodes
				ii=ii+1
                if (node(ndepi+j)%marge==0) node(ii)%marge=0
				node(ii)%x=node(ndepi+j)%x+p2(1);node(ii)%y=node(ndepi+j)%y+p2(2);node(ii)%z=node(ndepi+j)%z
				node(ii)%tipus=node(ndepi+j)%tipus
				cels(valcel)%nunodes=cels(valcel)%nunodes+1
				cels(valcel)%node(cels(valcel)%nunodes)=ii
				node(ii)%icel=valcel
			end do
			jj=i-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
			if(jj>0)then
				vector=p2-p1
				modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
!				iii=jj+1
				vector=vector/(jj+1)
				do kk=1,jj
					valcel=valcel+1
					do j=1,cels(ncelsepi+1)%nunodes
						kkk=cels(valcel-1)%node(j)
						ii=ii+1
                        if (node(ndepi+j)%marge==0) node(ii)%marge=0
						node(ii)%x=node(kkk)%x-vector(1);node(ii)%y=node(kkk)%y-vector(2);node(ii)%z=node(kkk)%z
						node(ii)%tipus=node(ndepi+j)%tipus
						cels(valcel)%nunodes=cels(valcel)%nunodes+1
						cels(valcel)%node(cels(valcel)%nunodes)=ii
						node(ii)%icel=valcel
					end do
				end do
			end if
			p1=p2
		end do
		angle=angle+beta
!		p2(1)=0;p2(2)=a*(i-1);p2(3)=0
		p2(1)=a*(i-1)*dcos(angle);p2(2)=a*(i-1)*dsin(angle);p2(3)=0
		valcel=valcel+1
		do j=1,cels(ncelsepi+1)%nunodes
			ii=ii+1
            if (node(ndepi+j)%marge==0) node(ii)%marge=0
			node(ii)%x=node(ndepi+j)%x+p2(1);node(ii)%y=node(ndepi+j)%y+p2(2);node(ii)%z=node(ndepi+j)%z
			node(ii)%tipus=node(ndepi+j)%tipus
			cels(valcel)%nunodes=cels(valcel)%nunodes+1
			cels(valcel)%node(cels(valcel)%nunodes)=ii
			node(ii)%icel=valcel
		end do
		jj=i-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
		if(jj>0)then
			vector=p2-p1
			modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
!			iii=jj+1
			vector=vector/(jj+1)
			do kk=1,jj
				valcel=valcel+1
				do j=1,cels(ncelsepi+1)%nunodes
					kkk=cels(valcel-1)%node(j)
					ii=ii+1
                    if (node(ndepi+j)%marge==0) node(ii)%marge=0
					node(ii)%x=node(kkk)%x-vector(1);node(ii)%y=node(kkk)%y-vector(2);node(ii)%z=node(kkk)%z
					node(ii)%tipus=node(ndepi+j)%tipus
					cels(valcel)%nunodes=cels(valcel)%nunodes+1
					cels(valcel)%node(cels(valcel)%nunodes)=ii
					node(ii)%icel=valcel
				end do
			end do
		end if
	end do
end if

valcelo=valcel
if(layer>1)then
  do i=2,layer
    kk=0
    do jj=1,mesradicel(i)-1
      kk=kk+jj
    end do
    kk=(6*kk+1)  !number of mesenchymal cells
    do j=ncelsepi+1,ncelsepi+kk
      valcel=valcel+1
      do k=1,cels(j)%nunodes
        kk=cels(j)%node(k) !;print*,"cel",valcel,"node",kk
        ii=ii+1
        if (node(kk)%marge==0) node(ii)%marge=0
        node(ii)%x=node(kk)%x;node(ii)%y=node(kk)%y;node(ii)%z=node(kk)%z-radi*(i-1)*di
        node(ii)%tipus=node(kk)%tipus
        node(ii)%icel=valcel
        cels(valcel)%nunodes=cels(valcel)%nunodes+1
        cels(valcel)%node(cels(valcel)%nunodes)=ii
      end do
    end do
  end do
end if


	!define the cell's centroid

	do i=ncelsepi+1,ncels
		cels(i)%ctipus=3
		cels(i)%cex=0;cels(i)%cey=0;cels(i)%cez=0
		do j=1,cels(i)%nunodes
			k=cels(i)%node(j)
 			cels(i)%cex=cels(i)%cex+node(k)%x
			cels(i)%cey=cels(i)%cey+node(k)%y
			cels(i)%cez=cels(i)%cez+node(k)%z
		end do
		cels(i)%cex=cels(i)%cex/real(cels(i)%nunodes)
		cels(i)%cey=cels(i)%cey/real(cels(i)%nunodes)
		cels(i)%cez=cels(i)%cez/real(cels(i)%nunodes)
	end do


!	do i=1,ndepi		!veins i parametres mecanics
!		ii=node(i)%icel
!		if(node(i)%tipus==1)then
!   		  node(i)%altre=i-cels(1)%nunodes/2
!		else
!		  node(i)%altre=i+cels(1)%nunodes/2
!		end if
!	end do

!        node(1)%marge=0.0d0

        !now we want the nucleus to be the first node in cels%node: this allows diffusion to be faster >>> Is 14-9-13 
	do i=ncelsepi+1,ncels
   	  do j=1,cels(i)%nunodes							
            k=cels(i)%node(j)								
            if (node(k)%marge==0) then
              ii=k
              exit
            end if
	  end do										
          jjj=cels(i)%node(1)
          cels(i)%node(1)=ii
   	  do j=1,cels(i)%nunodes							
            k=cels(i)%node(j)								
            if (k==ii) then
              jj=j
              exit
            end if
	  end do									
          cels(i)%node(jj)=jjj
	end do											


    !some funny rotation to make a less biased diffusion lost
    do i=ndepi+1,nd
      c=0.75
      a=node(i)%x*cos(c)-node(i)%y*sin(c)
      b=node(i)%x*sin(c)+node(i)%y*cos(c)
      node(i)%x=a
      node(i)%y=b
    end do

    node(:)%talone=0.0d0


!do i=ndepi+1,nd
! print*,"tipus",node(i)%tipus
!end do

end subroutine




subroutine matrix(radi,layer,zepi)		!to place a block of ECM
integer            ::radi,radicel,valcel,layer
real*8,dimension(:)::p1(3),p2(3),vector(3)
real*8             ::beta,angle,alt,de,di,zepi

    di=1.0d0
    de=2d0*node(1)%req  !miguel 14-10-13
    vector=0.0d0
    beta=2d0*pi/6d0



    ii=ndepi+ndmes
    node(ii+1:nd)%marge=1  !so there are no nuclei among ECM

    do i=1,layer
      ii=ii+1
      alt=zepi-real(i-1)*2*node(1)%req
      node(ii)%x=0.0d0;node(ii)%y=0d0;node(ii)%z=alt
      node(ii)%tipus=4
      node(ii)%icel=-ii
      do j=2,radi
        angle=beta
        d=de*(j-1d0)
        ii=ii+1
        node(ii)%x=d;node(ii)%y=0;node(ii)%z=alt
        node(ii)%tipus=4
        node(ii)%icel=-ii
        p1(1)=d;p1(2)=0;p1(3)=alt
        p2(1)=d*dcos(angle);p2(2)=d*dsin(angle);p2(3)=alt
        jj=j-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
        if(jj>0)then
          vector=p2-p1
          modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
          iii=jj+1
          vector=vector/iii
          do kk=1,jj
            ii=ii+1
            node(ii)%x=p1(1)+vector(1)*kk;node(ii)%y=p1(2)+vector(2)*kk;node(ii)%z=alt
            node(ii)%tipus=4
            node(ii)%icel=-ii
          end do
        end if
        do k=2,5
          angle=angle+beta
          p1=p2
          ii=ii+1
          node(ii)%x=p1(1);node(ii)%y=p1(2);node(ii)%z=alt
          node(ii)%tipus=4
          node(ii)%icel=-ii
          p2(1)=d*dcos(angle);p2(2)=d*dsin(angle);p2(3)=alt
          jj=j-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
          if(jj>0)then
            vector=p2-p1
            modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
            iii=jj+1
            vector=vector/iii
            do kk=1,jj
              ii=ii+1
              node(ii)%x=p1(1)+vector(1)*kk;node(ii)%y=p1(2)+vector(2)*kk;node(ii)%z=alt
              node(ii)%tipus=4
              node(ii)%icel=-ii
            end do
          end if
        end do
        angle=angle+beta
        p1=p2
        ii=ii+1
        node(ii)%x=p1(1);node(ii)%y=p1(2);node(ii)%z=alt
        node(ii)%tipus=4
        node(ii)%icel=-ii
        p2(1)=d;p2(2)=0;p2(3)=alt
        jj=j-2;if(jj<0)jj=0	!quants nodes intermitjos en la capa
        if(jj>0)then
          vector=p2-p1
          modul=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
          iii=jj+1
          vector=vector/iii
          do kk=1,jj
            ii=ii+1
            node(ii)%x=p1(1)+vector(1)*kk;node(ii)%y=p1(2)+vector(2)*kk;node(ii)%z=alt
            node(ii)%tipus=4
            node(ii)%icel=-ii
          end do
        end if
      end do
    end do

    node(:)%talone=0.0d0

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111

subroutine mesenq_cell_sorting(radi,radicel)									!!! miguel 8-7-2013 

integer  ::i,j,k,ii,jj,radi,radicel,layer,signo,signo1,ent,ont       ! number of layers and concentric hexagons 
real*8   ::rad,der,dar,zmes,de,di
real*8   :: ka,kb,kc,xx,yy,zz                                                        ! miguel 4-6-13
rad=pi/180d0
de=0.4d0                     ! cell radius
di=5.5d0*de                    ! radius of the bulk of cells

	do i=ncelsepi+1,ncels               !	
		cels(i)%nunodes=radi				!>>>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela     !    
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%node=0
    !cels(i)%tipcel=1       
  end do
         signo1=1 ; ii=1
         do i=1,radicel						
                            call random_number(a)
						dar=(sqrt(di)*sqrt(di*a))                            ! miguel 4-6-13
122         call random_number(a)         
            xx=dar*(1d0-2*a) ; 
            call random_number(a)
            yy=dar*(1d0-2*a)         ! miguel 4-6-13 
            zz=(dar**2)-(xx**2)-(yy**2) 			               						  ! miguel 4-6-13
            if(zz.lt.0)then;goto 122;endif                                ! miguel 4-6-13
            zz=signo1*sqrt(zz) ; signo1=-1*signo1                            ! miguel 4-6-13  
            if(mod(i,3).eq.0)then      ; kb=xx ; kc=yy ; ka=zz              ! miguel 4-6-13
            else if(mod(i,3).eq.1)then ; ka=xx ; kc=yy ; kb=zz              ! miguel 4-6-13
  	        else if(mod(i,3).eq.2)then ; ka=xx ; kb=yy ; kc=zz ; end if     ! miguel 4-6-13
            signo=1
            call random_number(a) 
						do kk=1,radi			
            call random_number(a)			
							der=(sqrt(de)*sqrt(de*a))                               ! miguel 4-6-13
123                         call random_number(a)           
                            xx=der*(1d0-2*a) ; 
                            call random_number(a)
                            yy=der*(1d0-2*a)            ! miguel 4-6-13 
                            zz=(der**2)-(xx**2)-(yy**2) 									 ! miguel 4-6-13
                            if(zz.lt.0)then;goto 123;endif                                   ! miguel 4-6-13
                            zz=signo*sqrt(zz) ; signo=-1*signo                               ! miguel 4-6-13  
                            if(mod(ii,3).eq.0)then      ; b=xx ; c=yy ; a=zz               ! miguel 4-6-13
                            else if(mod(ii,3).eq.1)then ; a=xx ; c=yy ; b=zz               ! miguel 4-6-13
  	                        else if(mod(ii,3).eq.2)then ; a=xx ; b=yy ; c=zz ; end if      ! miguel 4-6-13
							node(ii)%x=ka+a
							node(ii)%y=kb+b
							node(ii)%z=kc+c
					    node(ii)%tipus=3 ; node(ii)%icel=i !;write(*,*)'cel',i,'node',ii,'en orden',kk,'r,rcel',radi,radicel
							cels(i)%node(kk)=ii ; 	ii=ii+1
						end do
         end do

	!define the cell's centroid								!>>>>>>>>>>>>>> Miquel 14-4-13

	do i=ncelsepi+1,ncels 							!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%ctipus=3									!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%cex=0;cels(i)%cey=0;cels(i)%cez=0			!>>>>>>>>>>>>>> Miquel 14-4-13
		do j=1,cels(i)%nunodes								!>>>>>>>>>>>>>> Miquel 14-4-13
			k=cels(i)%node(j)								!>>>>>>>>>>>>>> Miquel 14-4-13
			cels(i)%cex=cels(i)%cex+node(k)%x				!>>>>>>>>>>>>>> Miquel 14-4-13
			cels(i)%cey=cels(i)%cey+node(k)%y				!>>>>>>>>>>>>>> Miquel 14-4-13
			cels(i)%cez=cels(i)%cez+node(k)%z				!>>>>>>>>>>>>>> Miquel 14-4-13
		end do												!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%cex=cels(i)%cex/real(cels(i)%nunodes)		!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%cey=cels(i)%cey/real(cels(i)%nunodes)		!>>>>>>>>>>>>>> Miquel 14-4-13
		cels(i)%cez=cels(i)%cez/real(cels(i)%nunodes)		!>>>>>>>>>>>>>> Miquel 14-4-13
	end do													!>>>>>>>>>>>>>> Miquel 14-4-13

    prop_noise=0.1d0    
    realtime=0
    reqmin=0.05d0 

    mmae=node(1)%req

    cels(:)%fase=0.0d0
    cels(:)%minsize_for_div=cels(:)%nunodes*2
    cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
    dmax=1
    screen_radius=1.0d0
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%da)*3
    node(:)%diffe=0.0d0
end subroutine mesenq_cell_sorting
!*************************************************************************************************************

!***********************************************************************

subroutine escribe
 open(666,file="nodoscels.dat",action='write')
    do i=1,size(node)
    write(666,*)'node',i,'tipus',node(i)%tipus,'cel',node(i)%icel,'xyz',node(i)%x,node(i)%y,node(i)%z
    end do
    do i=1,size(cels)
    write(666,*)'cel',i,'nodes',cels(i)%node(:)
    end do
  close(666)
end subroutine escribe

!***********************************************************************

subroutine migration

!print *,""
!print *,"WARNING THIS IS GOING TO TAKE A WHILE"
!print *,""


!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=3    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=14     !number of nodes per cell
	mradicel=2   !number of radial cell layers
	layer=1      !number of planar cell layers
	zmes=1.7d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=1
  !End initializing dimension parameters

  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.0 !1d-3 !low value is low temperature	
    desmax=0.0d0
    resmax=1d-3
    prop_noise=0.1d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=0.7d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !to quite if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !eggshell miguel4-1-13
    ffu(12)=1
    ffu(22)=0 !0 = noise biased by energies , 1 = unbiased noise

    ffu(9)=1   !euler or runge-kutta
    ffu(19)=0  !adaptive delta


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1d1;node(i)%adh=0d0	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%repcel=1d1 !>>Miquel 8-10-12
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req !;node(i)%reqcel=0.28d0	!>>Miquel 26-10-12!de
        node(i)%da=node(i)%req*1.40 
        node(i)%reqs=0.25d0
        node(i)%ke=1d1
        node(i)%tor=1d1
        node(i)%stor=1d1
        node(i)%mo=0d0
        node(i)%dmo=0d0
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=6d1;node(i)%adh=0d0 	!>>Miquel 26-10-12
        node(i)%rep=2d1;node(i)%repcel=3.0d1	!>>Miquel 8-10-12
        node(i)%req=0.25d0 ;node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*2.0d0 
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=1d0  !only for epithelium
        !node(i)%stor=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=1d0
        node(i)%dmo=5d-2
      end do
    end if




    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
      !Special modifications
      do i=ndepi+1,nd
        node(i)%y=node(i)%y+1.50d0
      end do


    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%hold=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=3
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=
       gen(1)%kindof=1 ; gen(1)%diffu=0d0 ; gen(1)%mu=0d0
       gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=0d0   !mesenchyme
       gen(3)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=0d0   !epithelium


    !Gene-behavior interactions
      gen(1)%wa(1)=1  !adhesion
      gen(2)%wa(1)=2  !adhesion
      gen(3)%wa(1)=3  !adhesion

      !gen(2)%wa(6)=0.1d0           !this is da (keep expression levels at one so this will be the actual value
      !gen(2)%wa(16)=0.1d0   !this is dmo (keep expression levels at one so this will be the actual value

    !Gene-gene interactions

       !****a rellenar****

    !Adhesion molecules

	ntipusadh=3
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

       kadh(1,2)=3d1 !adhesion table
       kadh(2,1)=3d1 !adhesion table THIS HAS TO BE SYMMETRIC
       kadh(3,3)=5.0d1
       kadh(2,2)=1d0  ! HI HA ALGO QUE NO VA
       !kadh(3,2)=1d1  ! HI HA ALGO QUE NO VA
       !kadh(2,3)=1d1  ! HI HA ALGO QUE NO VA


    end if

    !Gene expression on nodes
       do i=1,ndepi
         gex(i,3)=1.0d0
       end do
       do i=ndepi+1,nd
         gex(i,2)=1d0
         !gex(i,3)=0d0
       end do

       !adhesive concentration gradient in epithelium
       a=minval(node(:)%x)
       do i=1,ndepi
         if (node(i)%tipus==1) then
           gex(i,1)=1d0*((node(i)%y+a))**2
         end if
       end do

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%da)*3
    node(:)%diffe=0.0d0
 
    do i=1,ndepi
      node(i)%hold=2
    end do

!    node(ndepi+1:nd)%z=node(ndepi+1:nd)%z-0.3

end subroutine

!************************************************************************************************************************

subroutine epi_apoptosis   !>>>>>>>>>>>> Miquel 18-6-13


!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=4    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=0     !number of nodes per cell
	mradicel=0   !number of radial cell layers
	layer=0      !number of planar cell layers
	zmes=0.0d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=1
  !End initializing dimension parameters




    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d1 !low value is low temperature	
    desmax=0.001
    resmax=1d-3
    prop_noise=0.1d0
    deltamax=1d-2 ! miguel 14-10-13
    dmax=1
    screen_radius=1.0d0


    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !to quite if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !eggshell miguel4-1-13
    ffu(12)=1
    ffu(13)=0 !neighboring by triangulation
    ffu(22)=0 !0 = noise biased by energies , 1 = unbiased noise

    ffu(9)=1 !integration method 0=euler , 1=runge-kutta forces , 2=runge-kutta forces+genes
    ffu(19)=0 !adaptive time step 0=no , 1=yes for forces , 2=yes for forces+genes

    
   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=2.0d1;node(i)%adh=8d0	!>>Miquel 26-10-12
        node(i)%rep=1.0d1;node(i)%repcel=1d1	!>>Miquel 8-10-12
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req !;node(i)%reqcel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%reqs=0.5d0
		node(i)%da=node(i)%req*1.50 
        node(i)%ke=1d1
        node(i)%tor=1d0
        node(i)%stor=1d0
        node(i)%mo=1d0
        node(i)%dmo=5d-3
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=0d0;node(i)%adh=0d0 	!>>Miquel 26-10-12
        node(i)%rep=0d0;node(i)%repcel=0d0	!>>Miquel 8-10-12
        node(i)%req=0d0 !;node(i)%reqcel=0d0	!>>Miquel 26-10-12!de
        node(i)%reqs=0d0
        node(i)%da=0d0
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=1d0  !only for epithelium
        node(i)%stor=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%hold=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    !Number of genes
    ng=2
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=
       gen(1)%kindof=2 ; gen(1)%diffu=0d0 ; gen(1)%mu=0d0
       gen(1)%npost=1 ; allocate(gen(1)%post(gen(1)%npost))
       gen(1)%post(1)=2       
       gen(2)%kindof=3 ; gen(2)%diffu=0d0 ; gen(2)%mu=0d0
       !gen(2)%npost=1 ; allocate(gen(2)%post(gen(2)%npost))
       !gen(2)%post(1)=2       
       gen(2)%npre=1 ; allocate(gen(2)%pre(gen(2)%npre))
       gen(2)%pre(1)=1       


    !Gene-behavior interactions
      gen(1)%wa(nparam_per_node+3)=5d-1   !apoptosis

    !Gene-gene interactions

      !gen(1)%w(2)=1d0
      !gen(2)%w(1)=1d0
    
    
       gen(1)%nww=1
       gen(1)%ww(1,1)=1
       gen(1)%ww(1,2)=2
       gen(1)%ww(1,3)=1d0

    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes
      do i=1,cels(1)%nunodes ; gex(cels(1)%node(i),1)=1.0d0 ; end do
      do i=1,cels(8)%nunodes ; gex(cels(8)%node(i),1)=1.0d0  ; end do
      do i=1,cels(10)%nunodes ; gex(cels(10)%node(i),1)=1.0d0  ; end do
      do i=1,cels(12)%nunodes ; gex(cels(12)%node(i),1)=1.0d0  ; end do
      do i=1,cels(14)%nunodes ; gex(cels(14)%node(i),1)=1.0d0  ; end do
      do i=1,cels(16)%nunodes ; gex(cels(16)%node(i),1)=1.0d0  ; end do
      do i=1,cels(18)%nunodes ; gex(cels(18)%node(i),1)=1.0d0  ; end do

      gex(:14,1)=1.0d0
      do i=1,cels(1)%nunodes
        j=cels(1)%node(i)
        gex(j,1)=1.0d0
      end do



    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%da)*3
    node(:)%diffe=0.0d0
end subroutine epi_apoptosis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine epiteli_sphere(radi,radicel,radius)		!radi is the number of cell rows from one pole to the equator
integer            ::radi,radicel,valcel

integer::trobat,cont,rep,n,i,j,k,ii,jj,kk,lev,numdos,ncels2,val,val2,ll,val3,val4,val5,val6,iiii,jjjj,kkkk,suma,cocels
real*8::modul,pescab,pescac,pescbc,modula,modulb,modulc,aax,aay,aaz,bbx,bby,bbz,ccx,ccy,ccz,costat,rhex,lad,ax,ay,az
real*8::alf,bet,gam,a,b,c,l,aa,bb,cc,aaa,bbb,ccc,angle,rpent,modmin,pesc,radius,radius2,bx,by,bz,cost,ucost,sint,thet
real*8,dimension(:)::minu(3),vec1(3),vec2(3),vec3(3),vec4(3),vec5(3)
real*8,allocatable::cmalla(:,:),cveci(:,:),primers(:)




    di=1.0d0
    de=0.5d0

    radicel=(radicel-1)*2

    !SPHERE CODE, FIRST HEMISPHERE
	cont=0 ; cocels=0
	ii=radicel/2
!	radius=5d0  !radius of the sphere
    radius2=radius+di

	alf=pi/radicel

!!!pole cell
	cont=cont+1 ; cocels=cocels+1
    node(1)%x=0d0 ; node(1)%y=0d0 ; node(1)%z=radius
    node(1)%altre=2 ; node(1)%tipus=2 ; node(1)%icel=cocels
    cont=cont+1
    node(2)%x=0d0 ; node(2)%y=0d0 ; node(2)%z=radius2
    node(2)%altre=1 ; node(2)%tipus=1 ; node(2)%icel=cocels

    !el radi 2 de la cèlula

    cont=cont+1
    node(cont)%x=0d0 ; node(cont)%y=de ; node(cont)%z=radius
    node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
    cont=cont+1
    node(cont)%x=0d0 ; node(cont)%y=de ; node(cont)%z=radius2
    node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels

    gam=2*pi/6d0
    do i=1,5
      cont=cont+1
      node(cont)%x=de*cos(i*gam+pi/2) ; node(cont)%y=de*sin(i*gam+pi/2) ; node(cont)%z=radius
      node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
      cont=cont+1
      node(cont)%x=de*cos(i*gam+pi/2) ; node(cont)%y=de*sin(i*gam+pi/2) ; node(cont)%z=radius2
      node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels
    end do

!!!!! end pole cell


	do i=1,ii

        !!!!! spine cell
!  print*,"*********espina******",i
		cont=cont+1 ; cocels=cocels+1
		suma=0
		jj=cont
		angle=alf*i
!		print*,"cont",cont
		vec1(1)=radius*sin(angle)
		vec1(2)=0
		vec1(3)=radius*cos(angle)
		node(cont)%x=vec1(1) ; node(cont)%y=vec1(2) ; node(cont)%z=vec1(3)
        node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
!        print*,"altre",node(cont)%altre

        cont=cont+1
		vec2(1)=radius2*sin(angle)
		vec2(2)=0
		vec2(3)=radius2*cos(angle)
		node(cont)%x=vec2(1) ; node(cont)%y=vec2(2) ; node(cont)%z=vec2(3)
        node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels

        ux=vec2(1)-vec1(1) ; uy=vec2(2)-vec1(2) ; uz=vec2(3)-vec1(3)
        a=1d0/sqrt(ux**2+uy**2+uz**2) ; ux=ux*a ; uy=uy*a ; uz=uz*a

        ax=-de*sin(0d0) ; ay=de*cos(0d0) ; az=0d0
!        print*,"vora"," a",ax,ay,az


        do j=1,6
          cont=cont+1
          thet=j*gam
!   print*,"thet",thet
          cost=cos(thet); ucost=1-cost ; sint=sin(thet)
!  print*,"cos",cost,"ucost",ucost,"sint",sint
          bx=(cost+ux**2*ucost)*ax+(ux*uy*ucost-uz*sint)*ay+(ux*uz*ucost+uy*sint)*az
          by=(uy*ux*ucost+uz*sint)*ax+(cost+uy**2*ucost)*ay+(uy*uz*ucost-ux*sint)*az
          bz=(ux*uz*ucost-uy*sint)*ax+(uy*uz*ucost+ux*sint)*ay+(cost+uz**2*ucost)*az
          node(cont)%x=vec1(1)+bx ; node(cont)%y=vec1(2)+by ; node(cont)%z=vec1(3)+bz
          node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
!   print*,"b",bx,by,bz
!  print*,"vec1",vec1
!  print*,"nodecont",node(cont)%x,node(cont)%y,node(cont)%z
          cont=cont+1

          node(cont)%x=vec2(1)+bx ; node(cont)%y=vec2(2)+by ; node(cont)%z=vec2(3)+bz
          node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels
        end do

        !!!!!!!!! end spine cell



!        print*,"altre",node(cont)%altre

		if(i+1<=ii+1)then
			kk=6+6*(i-1)
		else
			kk=6+6*(radicel-i-1)
		end if
		bet=2*pi/kk

!		print*,"i",i,"nombre de cels del paralel",kk

		do j=1,kk-1

            !!!!!!!!!!! rib cell

			angle=alf

			cont=cont+1 ; cocels=cocels+1
			a=node(1)%z*sin(i*angle)!+node(1)%x*cos(i*angle)!fem els 2 girs respecte al node 1 (pol nord)
			b=node(1)%y
			c=node(1)%z*cos(i*angle)!-node(1)%x*sin(i*angle)
			node(cont)%x=a
			node(cont)%y=b
			node(cont)%z=c
			a=node(cont)%x*cos(j*bet)-node(cont)%y*sin(j*bet)
			b=node(cont)%x*sin(j*bet)+node(cont)%y*cos(j*bet)
			c=node(cont)%z
			node(cont)%x=a
			node(cont)%y=b
			node(cont)%z=c
            node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
            vec1(1)=node(cont)%x ; vec1(2)=node(cont)%y ; vec1(3)=node(cont)%z

			cont=cont+1
			a=node(2)%z*sin(i*angle)!+node(1)%x*cos(i*angle)!fem els 2 girs respecte al node 1 (pol nord)
			b=node(2)%y
			c=node(2)%z*cos(i*angle)!-node(1)%x*sin(i*angle)
			node(cont)%x=a
			node(cont)%y=b
			node(cont)%z=c
			a=node(cont)%x*cos(j*bet)-node(cont)%y*sin(j*bet)
			b=node(cont)%x*sin(j*bet)+node(cont)%y*cos(j*bet)
			c=node(cont)%z
			node(cont)%x=a
			node(cont)%y=b
			node(cont)%z=c
            node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels

            ux=node(cont)%x-node(cont-1)%x
            uy=node(cont)%y-node(cont-1)%y
            uz=node(cont)%z-node(cont-1)%z
            a=1d0/sqrt(ux**2+uy**2+uz**2) ; ux=ux*a ; uy=uy*a ; uz=uz*a !;print*,"spring",a
            vec2(1)=node(cont)%x ; vec2(2)=node(cont)%y ; vec2(3)=node(cont)%z

            ax=-de*sin(j*bet) ; ay=de*cos(j*bet) ; az=0d0

!        ax=de*sin(0d0) ; ay=de*cos(0d0) ; az=0d0
!        print*,"bet",j*bet,"a",ax,ay,"modul",sqrt(ax**2+ay**2)
            do k=1,6
              cont=cont+1
              thet=k*gam
              cost=cos(thet); ucost=1-cost ; sint=sin(thet)
              bx=(cost+ux**2*ucost)*ax+(ux*uy*ucost-uz*sint)*ay+(ux*uz*ucost+uy*sint)*az
              by=(uy*ux*ucost+uz*sint)*ax+(cost+uy**2*ucost)*ay+(uy*uz*ucost-ux*sint)*az
              bz=(ux*uz*ucost-uy*sint)*ax+(uy*uz*ucost+ux*sint)*ay+(cost+uz**2*ucost)*az
              node(cont)%x=vec1(1)+bx ; node(cont)%y=vec1(2)+by ; node(cont)%z=vec1(3)+bz
              node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels

!              print*,"b rib",bx,by,bz,"modul",sqrt(bx**2+by**2+bz**2)

              cont=cont+1

              node(cont)%x=vec2(1)+bx ; node(cont)%y=vec2(2)+by ; node(cont)%z=vec2(3)+bz
              node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels
           end do


            !!!!!!!!!!! end rib cell



		end do
	end do

!print*,"cont",cont,"cocels",cocels


!the 2ond hemisphere

    j=0
    do i=1,(radicel/2+1)-2
      j=j+i
    end do
    ii=(6*j+1)*2*7   !number of nodes on the 2ond hemisphere

    print*,"ii",ii
    do i=1,ii
      cont=cont+1
      node(cont)%x=node(i)%x
      node(cont)%y=node(i)%y
      node(cont)%z=-node(i)%z
      node(cont)%tipus=node(i)%tipus
      if(node(i)%tipus==2)then
        node(cont)%altre=cont+1
      else
        node(cont)%altre=cont-1
      end if
      node(cont)%icel=node(i)%icel+cocels
    end do

!print*,"cont complet",cont,"cocels",cocels


!	cont=cont+1
!	malla(cont,1)=0;malla(cont,2)=0;malla(cont,3)=-radius
!	primers(radi+1)=cont
!	print*,i,"malla",malla(cont,:)
!	print*,"nombre de cels creades",cont
!	print*,"nombre de paral·lels",lev





	!define the cell's centroid

    do i=1,ncelsepi
      kk=0
      cels(i)%ctipus=1
      cels(i)%nunodes=nodecel
      cels(i)%nodela=nodecela
      cels(i)%cex=0d0 ; cels(i)%cey=0d0 ; cels(i)%cez=0d0
      cels(i)%polx=0d0 ; cels(i)%poly=0d0 ; cels(i)%polz=0d0
      allocate(cels(i)%node(nodecela))
      do j=1,ndepi
        k=node(j)%icel
        if(k==i)then
          kk=kk+1
          cels(i)%node(kk)=j
          if(node(j)%tipus==1)then
            cels(i)%cex=cels(i)%cex+node(j)%x
            cels(i)%cey=cels(i)%cey+node(j)%y
            cels(i)%cez=cels(i)%cez+node(j)%z
          end if
        end if
      end do
      cels(i)%cex=2*cels(i)%cex/cels(i)%nunodes
      cels(i)%cey=2*cels(i)%cey/cels(i)%nunodes
      cels(i)%cez=2*cels(i)%cez/cels(i)%nunodes
    end do


!	do i=1,ncelsepi
!		cels(i)%ctipus=1
!		cels(i)%cex=0;cels(i)%cey=0;cels(i)%cez=0
!		do j=1,cels(i)%nunodes
!			k=cels(i)%node(j)
!            if(node(k)%tipus==1)then
!	 			cels(i)%cex=cels(i)%cex+node(k)%x
! 				cels(i)%cey=cels(i)%cey+node(k)%y
!				cels(i)%cez=cels(i)%cez+node(k)%z
!			end if
!		end do
!		cels(i)%cex=2*cels(i)%cex/real(cels(i)%nunodes)
!		cels(i)%cey=2*cels(i)%cey/real(cels(i)%nunodes)
!		cels(i)%cez=2*cels(i)%cez/real(cels(i)%nunodes)
!	end do


	do i=1,ndepi		!veins i parametres mecanics
		node(i)%you=1d1;node(i)%adh=1d1
		node(i)%rep=1d2;node(i)%repcel=1d2
		node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req !;node(i)%reqcel=0.25d0 !de
            node(i)%reqs=0.5d0
		node(i)%da=node(i)%req*1.25!+da 
        node(i)%ke=1d1
!        node(i)%tipus=3
!        node(i)%icel=i
	end do


    realtime=0
    dmax=1
    screen_radius=1.0d0
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%da)*3
    node(:)%diffe=0.0d0
end subroutine epiteli_sphere

!**********************************************************************************************

subroutine epi_sphere


	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=5    !number of radial layers of cells
    radius=3d0   !radius of the sphere
!	zepi=0d0     !z-position of basal layer


	j=0
	do i=1,radi-1
	  j=j+i
	end do
	nodecel=(6*j+1)*2	!number of nodes per cell

	nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

!    nodecel=1
!    nodecela=1


	j=0					!cell count
	do i=1,radicel-1
  	  j=j+i
	end do

    k=0
	do i=1,radicel-2
  	  k=k+i
	end do



	ncelsepi=(6*j+1)+(6*k+1) ; print*,"ncelsepi",ncelsepi,"j",j
    ndepi=nodecel*ncelsepi

	ncels=ncelsepi+ncelsmes

	ncelsmes=0

    nx=0  !number of ECM nodes

	ndmes=0

	nd=ndepi+ndmes+ndx ; print*,"nd",nd

    nda=nd+10
	ncals=ncels+10



	print*,"nd",nd,"ncels",ncels,"nodecela",nodecela

    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 

    call iniarrays


    call epiteli_sphere(radi,radicel,radius)


  ng=0       !>>>>>> Is 29-4-13

  !physical
  getot=0
  temp=0.1d-5 !low value is low temperature	
  !mathematical
  itacc=0
  nparti=1000
  desmax=0.01
        resmax=1d-3
  idum=-11111
  idumoriginal=idum
  rv=0.5d0*1.5d0	!el maxim rang d'abast que tindra un node   !it is redefined below as function of da
	print*,"rv",rv
  interfmax=1d0
    dmax=1
    screen_radius=1.0d0

  !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !to quite if there are too many cells
    ffu(3)=0 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source !miguel4-11-13
    ffu(6)=0 !eggshell miguel4-1-13
    ffu(22)=0 !0 = noise biased by energies , 1 = unbiased noise


    !number of params
    !nparam=32
    nvarglobal_out=5

	geu=0;rtime=0
    tacre=5	!growth rate and differential	>>Miquel 8-10-12

                    !no es el parametre en si si no una forma de veure'l

    ntipusadh=1


    if (allocated(kadh)) deallocate(kadh)
    allocate(kadh(ntipusadh,ntipusadh))
    kadh=0d0
    gen(1)%wa(1)=1

    do i=1,nd		!veins i parametres mecanics
      node(i)%you=1.3d1;node(i)%adh=8d0	!>>Miquel 26-10-12
      node(i)%rep=1.3d1;node(i)%repcel=1.8d1	!>>Miquel 8-10-12
      node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req !;node(i)%reqcel=0.25d0	!>>Miquel 26-10-12!de
      node(i)%reqs=0.5d0
      node(i)%da=node(i)%req*1.35  
      node(i)%ke=1d2
      node(i)%tor=1d0
      !node(i)%stor=1d1
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
      ii=node(i)%icel

      if(node(i)%tipus==1)then
        node(i)%req=0.35d0
!        node(i)%reqcel=0.35d0
        node(i)%da=node(i)%req*1.35d0
      end if
      node(i)%dmo=desmax
      node(i)%mo=temp

    end do

    rv=2*node(1)%da;  print*,"rv",rv

    ng=1

    call initiate_gene
    gen(1)%wa(15)=1d-4
    do i=1,nd
      gex(i,1)=0.0d0	!no growth
    end do

    call update_npag

    prop_noise=0.1d0    


    realtime=0
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%da)*3
    node(:)%diffe=0.0d0
end subroutine epi_sphere

!***************************************************************************************************************

subroutine epi_mes_ecm


!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=4    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=15     !number of nodes per cell
	mradicel=4   !number of radial cell layers
	layer=1      !number of planar cell layers
	zmes=-0.65d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************

    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=mradi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=1
  !End initializing dimension parameters


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-4 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    ecmmax=0.25d0
    dmax=1
    screen_radius=1.0d0
    deltamin=1d-3

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !to quite if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(8)=0
    ffu(11)=1
    ffu(12)=1
    ffu(13)=0
    ffu(17)=1
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise

    ffu(9)=0 !integration method 0=euler , 1=runge-kutta forces , 2=runge-kutta forces+genes
    ffu(19)=0 !adaptive time step 0=no , 1=yes for forces , 2=yes for forces+genes


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1d1;node(i)%adh=0d0	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%repcel=2d1	!>>Miquel 8-10-12
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req !;node(i)%reqcel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%reqs=0.25d0
        node(i)%da=node(i)%req*1.50
        node(i)%ke=1d1
        node(i)%tor=5d0
        node(i)%stor=1d1
        node(i)%mo=temp
        node(i)%dmo=desmax
        node(i)%kplast=1d0
        node(i)%kvol=1d0
        node(i)%acecm=0d0
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d1;node(i)%adh=0d0	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%repcel=2d1	!>>Miquel 8-10-12
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req !;node(i)%reqcel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%reqs=0d0
		node(i)%da=node(i)%req*1.90; 
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=desmax
        node(i)%acecm=0d0
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)

    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%hold=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=4
    call initiate_gene
    !Gene parameters

      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=
       gen(1)%kindof=1 ; gen(1)%diffu=1d1 ; gen(1)%mu=0d0 ; gen(1)%npost=0 ; gen(1)%npre=0
!       gen(1)%npre=1   ; allocate(gen(1)%pre(1)) ; gen(1)%pre(1)=5
       gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=0d0 ; gen(2)%npost=0 ; gen(2)%npre=0
       gen(3)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=0d0 ; gen(3)%npost=0 ; gen(3)%npre=0
       gen(4)%kindof=1 ; gen(4)%diffu=0d0 ; gen(4)%mu=0d0 ; gen(4)%npost=0 ; gen(4)%npre=0
!       gen(5)%kindof=1 ; gen(5)%diffu=5d0 ; gen(5)%mu=0d0 ; 
!       gen(5)%npost=1   ; allocate(gen(5)%post(1)) ; gen(5)%post(1)=1

    !Gene-behavior interactions
      gen(1)%wa(1)=1
      gen(2)%wa(1)=2
      !gen(3)%wa(1)=3
      gen(4)%wa(nparam_per_node+4)=2d-1 
      gen(1)%wa(nparam_per_node+5)=1d0
      gen(4)%wa(nparam_per_node+6)=2d1
      gen(4)%wa(nparam_per_node+7)=0.10d0

    !Gene-gene interactions
      gen(1)%w(4)=1.0d4
      !gene 4 induces the synthesis of the secreted gene

    !Adhesion molecules

	ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecule interactions

        kadh(1,1)=1d2
        kadh(1,2)=1d0 ; kadh(2,1)=kadh(1,2)
        !kadh(1,3)=1d0;kadh(3,1)=kadh(1,3)
        kadh(2,2)=1d1
        !kadh(2,3)=1d0 ; kadh(3,2)=kadh(2,3)
        !kadh(3,3)=1d1


    end if

    !Gene expression on nodes
      do i=1,nd
        if(node(i)%tipus==3)then
          gex(i,2)=1d0
          !gex(i,1)=1.0d0
          !gex(i,4)=1.0d0;!print*,"i",i
        elseif(node(i)%icel<8.and.node(i)%tipus==2) then
          !gex(i,2)=1d0
          !gex(i,1)=1.0d0
          gex(i,4)=1.0d0;!print*,"i",i
        end if
        if(node(i)%tipus<3)then
          gex(i,2)=1d0
        end if
      end do

    call update_npag

    node(:)%talone=0.0d0
    ramax=maxval(node(:)%da)*3
    node(:)%diffe=0.0d0
end subroutine epi_mes_ecm

!***************************************************************************************************


subroutine epi_polar_growth			!>>Miquel 14-10-12

!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=2    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=0     !number of nodes per cell
	mradicel=0   !number of radial cell layers
	layer=0      !number of planar cell layers
	zmes=0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=1
  !End initializing dimension parameters

!  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-5 !low value is low temperature	
    desmax=0.00
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    min_comp=-1d-1
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !to quite if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=1 !external signal source
    ffu(6)=0 !eggshell miguel4-1-13
    ffu(12)=1
    ffu(9)=0
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise


!ffu(13)=1
   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1.3d1;node(i)%adh=6d0	!>>Miquel 26-10-12
        node(i)%rep=1.3d1;node(i)%repcel=1.5d1	!>>Miquel 8-10-12
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req !;node(i)%reqcel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%reqs=0.25d0
        node(i)%da=node(i)%req*1.60  
        node(i)%ke=1d1
        node(i)%tor=1d0
        node(i)%stor=1d1
        node(i)%mo=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1.3d1;node(i)%adh=6d0 	!>>Miquel 26-10-12
        node(i)%rep=1.3d1;node(i)%repcel=1.5d1	!>>Miquel 8-10-12
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req !;node(i)%reqcel=0.25d0	!>>Miquel 26-10-12!de
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.40
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%hold=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=2
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu
      gen(1)%diffu=0d0 ; gen(1)%kindof=0 ; gen(1)%mu=0.0d0 ! miguel 14-10-13
      gen(2)%diffu=1d1 ; gen(2)%kindof=0 ; gen(2)%mu=0.0d0 ! miguel 14-10-13

    !Gene-behavior interactions
      gen(2)%wa(nparam_per_node+1)=5.0d-1
      gen(2)%wa(nparam_per_node+2)=1.0d0  
      gen(1)%wa(nparam_per_node+8)=1.0d0  
      gen(2)%wa(nparam_per_node+9)=0.5d0   
      gen(2)%wa(nparam_per_node+11)=1.0d0   

    !Gene-gene interactions
      gen(2)%w(2)=0.0d0

    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes
      gex(:,1)=abs(node(:)%x)
      gex(:,2)=1d0

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%da)*3
    node(:)%diffe=0.0d0
end subroutine

!***********************************************************************



subroutine epi_reaction_diffusion



!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=4    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
	mradi=0     !number of nodes per cell
	mradicel=0   !number of radial cell layers
	layer=0      !number of planar cell layers
	zmes=-1.0d0   !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=mradi*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=1
  !End initializing dimension parameters

!  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-5 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.01d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !to quite if there are too many cells
    ffu(3)=0 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(22)=0 !0 = noise biased by energies , 1 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        node(i)%you=1.2d1;node(i)%adh=8d1    	!>>Miquel 26-10-12
        node(i)%rep=1d1;node(i)%repcel=1d5	!>>Miquel 8-10-12
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req !;node(i)%reqcel=0.25d0;	!>>Miquel 26-10-12!de
        node(i)%reqs=0.5d0
        node(i)%da=node(i)%req*1.25; 
        node(i)%ke=1d2
        node(i)%tor=1d1
        !node(i)%stor=1d1
        node(i)%mo=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=0d0;node(i)%adh=0d0 	!>>Miquel 26-10-12
        node(i)%rep=0d0;node(i)%repcel=0d0	!>>Miquel 8-10-12
        node(i)%req=0d0 !;node(i)%reqcel=0d0	!>>Miquel 26-10-12!de
        node(i)%reqs=0d0
        node(i)%da=0d0
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=1d0  !only for epithelium
        !node(i)%stor=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    if(mradi>0.and.mradicel>0.and.layer>0) call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%hold=0


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=9
    call initiate_gene


    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=
       gen(1)%kindof=1 ; gen(1)%diffu=1d0 ; gen(1)%mu=1d-2 ; gen(1)%post=2  !the activator intracell form
       gen(2)%kindof=4 ; gen(2)%diffu=1d0 ; gen(2)%mu=1d-2 ; gen(2)%pre=1   !the activator extracell form
       gen(3)%kindof=1 ; gen(3)%diffu=1d0 ; gen(3)%mu=1d-2 ; gen(3)%post=4  !the activator receptor inactive form
       gen(4)%kindof=2 ; gen(4)%diffu=1d0 ; gen(4)%mu=1d-2 ; gen(4)%pre=3   !the activator receptor active form
       gen(5)%kindof=1 ; gen(5)%diffu=1d0 ; gen(5)%mu=1d-2 ; gen(5)%post=6  !the inhibitor intracell form
       gen(6)%kindof=4 ; gen(6)%diffu=5d0 ; gen(6)%mu=1d-2 ; gen(6)%pre=5   !the inhibitor extracell form
       gen(7)%kindof=1 ; gen(7)%diffu=1d0 ; gen(7)%mu=1d-2 ; gen(7)%post=8  !the inhibitor receptor inactive form
       gen(8)%kindof=2 ; gen(8)%diffu=1d0 ; gen(8)%mu=1d-2 ; gen(8)%pre=7   !the inhibitor receptor active form
       gen(9)%kindof=0 ; gen(9)%diffu=1d0 ; gen(9)%mu=0d0                   !house-keeping gene, induces expression of receptors

    !Gene-behavior interactions

       !****a rellenar****

    !Gene-gene interactions
      gen(1)%w(2)=1.0d0  !activator induces own secretion
      gen(2)%w(4)=1.0d0  !activator extracell activates activator receptor
      gen(4)%w(1)=5.0d0  !activator active receptor activates transcription of activator
      gen(4)%w(5)=1.0d0  !activator active receptor activates transcription of inhibitor
      gen(5)%w(6)=1.0d0  !inhibitor induces own secretion
      gen(6)%w(8)=1.0d0  !inhibitor extracell activates inhibitor receptor
      gen(8)%w(1)=-5.0d0 !inhibitor active receptor inhibits activator transcription
      gen(9)%w(3)=1.0d0  !house-keeping induces expression of the activator receptor
      gen(9)%w(7)=1.0d0  !house-keeping induces expression of the inhibitor receptor


    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes

      !one cell expressing activator in the centre
      gex(1,1)=1.0d0
      gex(:,3)=1.0d0
      gex(:,7)=1.0d0
      gex(:,9)=1.0d0

    
    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%da)*3
    node(:)%diffe=0.0d0
end subroutine

!***********************************************************************************

subroutine epi_mes_primordium



!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=2       !number of radial layers of nodes per cell
	radicel=4    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
	mradi=2      !number of nodes per cell    !if packed this is number of radial nodes per cell
	mradicel=4   !number of radial cell layers
	layer=2      !number of planar cell layers
	zmes=-1.0d0  !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi	!number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=1
  !End initializing dimension parameters

  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-2 !low value is low temperature	
    desmax=0.01
    resmax=1d-2
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    khold=1d6

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !to quite if there are too many cells
    ffu(3)=0 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !eggshell miguel4-1-13
    ffu(22)=0 !0 = noise biased by energies , 1 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(node(i)%tipus==2)then  !basal
          node(i)%you=1d1 ; node(i)%adh=1d1
          node(i)%rep=1d1 ; node(i)%repcel=1d1
        else                      !apical
          node(i)%you=1d1 ; node(i)%adh=1d1
          node(i)%rep=1d1 ; node(i)%repcel=1d1
        end if
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.25d0
        node(i)%da=node(i)%req*1.30; 
        node(i)%ke=1d1
        node(i)%tor=5d0
        !node(i)%stor=5d1
        node(i)%mo=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d1 ; node(i)%adh=1d1
        node(i)%rep=1d1 ; node(i)%repcel=1d1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.5d0
        node(i)%da=node(i)%req*1.15 
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=1d0  !only for epithelium
        !node(i)%stor=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=desmax
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if


    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%hold=0
    !let's do it for the most external layer of cells
      j=0
      do i=1,radicel-2
        j=j+i
      end do
      j=(6*j+1)*nodecel !this is the number of epithelial nodes wich are not external
      do i=j+1,ndepi
        node(i)%hold=1
        node(i)%rep=1d5;node(i)%repcel=1d5
        node(i)%ke=1d5
        node(i)%tor=1d5
        !node(i)%stor=1d5
      end do

!      j=0
!      do i=1,mradicel-2
!        j=j+i
!      end do
!      j=(6*j+1)*mradi !this is the number of mesenchymal nodes wich are not external
      do i=ndepi+j+1,ndepi+ndmes
        node(i)%hold=1
        node(i)%rep=1d5;node(i)%repcel=1d5
!        node(i)%req=0.30 ; node(i)%da=0.31 !I make them bigger so they make a wall and don't let anyone pass
      end do
    !end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=10
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=2 ; gen(1)%diffu=1d2 ; gen(1)%mu=5d-1 !activator transcript
      gen(1)%npost=1 ; allocate(gen(1)%post(gen(1)%npost))
      gen(1)%post(1)=2

      gen(2)%kindof=4 ; gen(2)%diffu=1d1 ; gen(2)%mu=5d-1 !activator signal
      gen(2)%npre=1 ; allocate(gen(2)%pre(gen(2)%npre))
      gen(2)%pre(1)=1

      gen(3)%kindof=2 ; gen(3)%diffu=1d2 ; gen(3)%mu=0d0 !free receptor
      gen(3)%npost=1 ; allocate(gen(3)%post(gen(3)%npost))
      gen(3)%post(1)=4

      gen(4)%kindof=3 ; gen(4)%diffu=1d2 ; gen(4)%mu=5d-1 !activated receptor
      gen(4)%npre=2 ; allocate(gen(4)%pre(gen(4)%npre))
      gen(4)%pre(1)=3 ; gen(4)%pre(2)=2

      gen(5)%kindof=1 ; gen(5)%diffu=1d2 ; gen(5)%mu=0d0 !TF for activator (inactive)
      gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost))
      gen(5)%post(1)=6

      gen(6)%kindof=3 ; gen(6)%diffu=1d2 ; gen(6)%mu=5d-1 !activated TF
      gen(6)%npre=1 ; allocate(gen(6)%pre(gen(6)%npre))
      gen(6)%pre(1)=5

!      gen(7)%kindof=1 ; gen(7)%diffu=1d1 ; gen(7)%mu=5d-1 !cell proliferation factor

      gen(7)%kindof=1 ; gen(7)%diffu=1d2 ; gen(7)%mu=0d0 !inactive cell proliferation factor (mesenchyme)
      gen(7)%npost=1 ; allocate(gen(7)%post(gen(7)%npost))
      gen(7)%post(1)=8

      gen(8)%kindof=3 ; gen(8)%diffu=1d2 ; gen(8)%mu=5d-1 !active cell proliferation factor (mesenchyme)
      gen(8)%npre=1 ; allocate(gen(8)%pre(gen(8)%npre))
      gen(8)%pre(1)=7

      gen(9)%kindof=1 ; gen(9)%diffu=1d2 ; gen(9)%mu=0d0 !inactive cell proliferation factor (epithelium)
      gen(9)%npost=1 ; allocate(gen(9)%post(gen(9)%npost))
      gen(9)%post(1)=10

      gen(10)%kindof=3 ; gen(10)%diffu=1d2 ; gen(10)%mu=5d-1 !active cell proliferation factor (epithelium)
      gen(10)%npre=1 ; allocate(gen(10)%pre(gen(10)%npre))
      gen(10)%pre(1)=9


    !Gene-behavior interactions

       gen(8)%wa(nparam_per_node+1)=1d-3  !growth
       gen(8)%wa(nparam_per_node+2)=1d-3  !cell division

       gen(10)%wa(nparam_per_node+1)=1d-5  !growth
       gen(10)%wa(nparam_per_node+2)=1d-5  !cell division


    !Gene-gene interactions
       gen(1)%nww=1
       gen(1)%ww(1,2)=1d1 !1 autocatalyzes transcription
       gen(2)%nww=1
       gen(2)%ww(3,4)=1d1 !signal activates receptor
       gen(4)%nww=3
       gen(4)%ww(5,6)=1d1 !active receptor activates TF
       gen(4)%ww(7,8)=1d1 !active receptor activates cell prolif. factor
       gen(4)%ww(9,10)=1d1 !active receptor activates cell prolif. factor
       gen(6)%w(1)=1d1 !active TF promotes transcirption


    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

         !****a rellenar****

    end if

    !Gene expression on nodes

    do i=1,nd
      gex(i,3)=1d3
      if(node(i)%tipus==3)then
        gex(i,5)=1d3
        gex(i,7)=1d3
        if(node(i)%marge==0)then
          gex(i,1)=1d0
        end if
      else
        gex(i,9)=1d3
      end if
    end do

    call update_npag

    node(:)%talone=0.0d0
    ramax=maxval(node(:)%da)*3
    node(:)%diffe=0.0d0
end subroutine


!*************************************************************************************

subroutine mes_primordium
integer:: val


!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
	radi=0       !number of radial layers of nodes per cell
	radicel=0    !number of radial layers of cells
	zepi=0d0     !z-position of basal layer

	!mesenchyme's dimension parameters
        packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
	mradi=2      !number of nodes per cell    !if packed this is number of radial nodes per cell
	mradicel=4   !number of radial cell layers
	layer=3      !number of planar cell layers
	zmes=0.0d0  !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi ; print*,"nodecel",nodecel	!number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells      
      val=(6*j+1)*nodecel ;print*,"val",val    !number of nodes per layer
      ndmes=nodecel*ncelsmes  !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        if(packed==0)then
          nodecel=mradi
        else
    !      nodecel=(6*j+1)*mradi
        end if
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
    nd=ndepi+ndmes
    ncels=ncelsepi+ncelsmes
    nda=nd+10
    ncals=ncels+10
    !End initializing dimension parameters

  print*,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes,"nd",nd,"ncels",ncels


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-2 !low value is low temperature	
    desmax=0.01
    resmax=1d-2
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    khold=1d2

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !to quite if there are too many cells
    ffu(3)=0 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !eggshell miguel4-1-13
    ffu(22)=0 !0 = noise biased by energies , 1 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(node(i)%tipus==2)then  !basal
          node(i)%you=1d1 ; node(i)%adh=1d1
          node(i)%rep=5d0 ; node(i)%repcel=5d0
        else                      !apical
          node(i)%you=5d0 ; node(i)%adh=5d0
          node(i)%rep=1d1 ; node(i)%repcel=1d1
        end if
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.25d0
        node(i)%da=node(i)%req*1.30; 
        node(i)%ke=1d2
        node(i)%tor=1d0
        !node(i)%stor=1d1
        node(i)%mo=temp
        node(i)%dmo=desmax
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d1;node(i)%adh=0d0 	
        node(i)%rep=1d2;node(i)%repcel=1d2	
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.5d0
        node(i)%da=node(i)%req*1.15 
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=1d0  !only for epithelium
        !node(i)%stor=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=desmax
      end do
      do i=1,val !this is the "epithelial like part"
        node(i)%da=node(i)%req*1.35
      end do
    end if
    
    
    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if


    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%hold=0
    !let's do it for the most external layer of cells
      !j=0
      !do i=1,radicel-2
      !  j=j+i
      !end do
      !j=(6*j+1)*nodecel !this is the number of epithelial nodes wich are not external
      !do i=j+1,ndepi
      !  node(i)%hold=1
      !end do

      j=0
      do i=1,mradicel-2
        j=j+i
      end do
      j=(6*j+1)*nodecel !this is the number of mesenchymal nodes wich are not external

      do i=ndepi+j+1,ndepi+val
        node(i)%hold=1
        node(i)%rep=1d6;node(i)%repcel=1d6
!        node(i)%req=0.30 ; node(i)%da=0.31 !I make them bigger so they make a wall and don't let anyone pass
      end do
      do i=ndepi+val+j+1,ndepi+ndmes
        node(i)%hold=1
        node(i)%rep=1d6;node(i)%repcel=1d6
!        node(i)%req=0.30 ; node(i)%da=0.31 !I make them bigger so they make a wall and don't let anyone pass
      end do


    !end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=12
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=2 ; gen(1)%diffu=1d1 ; gen(1)%mu=5d-1 !activator transcript
      gen(1)%npost=1 ; allocate(gen(1)%post(gen(1)%npost))
      gen(1)%post(1)=2

      gen(2)%kindof=4 ; gen(2)%diffu=1d1 ; gen(2)%mu=5d-1 !activator signal
      gen(2)%npre=1 ; allocate(gen(2)%pre(gen(2)%npre))
      gen(2)%pre(1)=1

      gen(3)%kindof=2 ; gen(3)%diffu=1d1 ; gen(3)%mu=0d0 !free receptor
      gen(3)%npost=1 ; allocate(gen(3)%post(gen(3)%npost))
      gen(3)%post(1)=4

      gen(4)%kindof=3 ; gen(4)%diffu=1d1 ; gen(4)%mu=5d-1 !activated receptor
      gen(4)%npre=2 ; allocate(gen(4)%pre(gen(4)%npre))
      gen(4)%pre(1)=3 ; gen(4)%pre(2)=2

      gen(5)%kindof=1 ; gen(5)%diffu=1d1 ; gen(5)%mu=0d0 !TF for activator (inactive)
      gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost))
      gen(5)%post(1)=6

      gen(6)%kindof=3 ; gen(6)%diffu=1d1 ; gen(6)%mu=5d-1 !activated TF
      gen(6)%npre=1 ; allocate(gen(6)%pre(gen(6)%npre))
      gen(6)%pre(1)=5

      gen(7)%kindof=1 ; gen(7)%diffu=1d1 ; gen(7)%mu=0d0 !inactive cell proliferation factor (mesenchyme)
      gen(7)%npost=1 ; allocate(gen(7)%post(gen(7)%npost))
      gen(7)%post(1)=8

      gen(8)%kindof=3 ; gen(8)%diffu=1d1 ; gen(8)%mu=5d-1 !active cell proliferation factor (mesenchyme)
      gen(8)%npre=1 ; allocate(gen(8)%pre(gen(8)%npre))
      gen(8)%pre(1)=7

      gen(9)%kindof=1 ; gen(9)%diffu=1d1 ; gen(9)%mu=0d0 !inactive cell proliferation factor (epithelium)
      gen(9)%npost=1 ; allocate(gen(9)%post(gen(9)%npost))
      gen(9)%post(1)=10

      gen(10)%kindof=3 ; gen(10)%diffu=1d1 ; gen(10)%mu=5d-1 !active cell proliferation factor (epithelium)
      gen(10)%npre=1 ; allocate(gen(10)%pre(gen(10)%npre))
      gen(10)%pre(1)=9

      gen(11)%kindof=1 ; gen(11)%diffu=1d1 ; gen(11)%mu=0d0 !epithelial adhesion molecule
      
      gen(12)%kindof=1 ; gen(12)%diffu=1d1 ; gen(12)%mu=0d0 !mesenchymal adhesion molecule

    !Gene-behavior interactions

       gen(8)%wa(nparam_per_node+1)=1d-2  !growth
       gen(8)%wa(nparam_per_node+2)=5d-2  !cell division

       gen(10)%wa(nparam_per_node+1)=5d-4  !growth
       gen(10)%wa(nparam_per_node+2)=1d-4  !cell division
       
       gen(11)%wa(1)=1
       gen(12)%wa(1)=2


    !Gene-gene interactions

       gen(1)%nww=1
       gen(1)%ww(1,2)=1d1 !1 autocatalyzes transcription
       gen(2)%nww=1
       gen(2)%ww(3,4)=1d1 !signal activates receptor
       gen(4)%nww=3
       gen(4)%ww(5,6)=1d1 !active receptor activates TF
       gen(4)%ww(7,8)=1d1 !active receptor activates cell prolif. factor
       gen(4)%ww(9,10)=1d1 !active receptor activates cell prolif. factor
       gen(6)%w(1)=1d1 !active TF promotes transcirption



    !Adhesion molecules

    ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=1d1
      kadh(1,2)=-1d1 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=1d1


    end if

    !Gene expression on nodes

    do i=1,nd
      gex(i,3)=1d3
    !  if(node(i)%tipus==3)then
    !    gex(i,5)=1d3
    !    gex(i,7)=1d3
    !    if(node(i)%marge==0)then
    !      gex(i,1)=1d0
    !    end if
    !  else
    !    gex(i,9)=1d3
    !  end if
    end do
    
    do i=1,val
      gex(i,9)=1d3
      gex(i,11)=1d0
    end do
    do i=val+1,nd
      gex(i,5)=1d3
      gex(i,7)=1d3
      gex(i,12)=1d0
      if(node(i)%marge==0) gex(i,1)=1d0
    end do

    call update_npag

    node(:)%talone=0.0d0

    !print*,"w de 1",gen(1)%w(:)
    !print*,"w de 2",gen(2)%w(:)
    !print*,"w de 3",gen(3)%w(:)
    !print*,"w de 4",gen(4)%w(:)
    !print*,"w de 5",gen(5)%w(:)
    !print*,"w de 6",gen(6)%w(:)
    !print*,"w de 7",gen(7)%w(:)
    !print*,"w de 8",gen(8)%w(:)
    !print*,"w de 9",gen(9)%w(:)
    !print*,"w de 10",gen(10)%w(:)
    !print*,"w de 11",gen(11)%w(:)
    !print*,"w de 12",gen(12)%w(:)
    ramax=maxval(node(:)%da)*3    
    node(:)%diffe=0.0d0
end subroutine

!!!!!!!!!***********************************SUBROUTINE********************************

subroutine epi_mes_bud



!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
    radi=2       !number of radial layers of nodes per cell
    radicel=4    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=2      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=2   !number of radial cell layers
    layer=2      !number of planar cell layers
    zmes=-1.0d0  !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi	!number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    if(radi==1.or.mradi==1) ffu(1)=1
  !End initializing dimension parameters

  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=0.1d-2 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    khold=1d0


    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !to quite if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !eggshell miguel4-1-13

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=0 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(22)=0 !0 = noise biased by energies , 1 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(i<=ndepi/2)then  !basal
          node(i)%you=5d0 ; node(i)%adh=0d0
          node(i)%rep=5d0 ; node(i)%repcel=5d0
        else                      !apical
          node(i)%you=5d0 ; node(i)%adh=0d0
          node(i)%rep=5d0 ; node(i)%repcel=5d0
        end if
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.25
        node(i)%da=node(i)%req*1.40; 
        node(i)%ke=1d2
        node(i)%tor=1d0
        !node(i)%stor=1d0
        node(i)%mo=temp
        node(i)%dmo=desmax
        node(i)%kplast=1d-2 ; node(i)%kvol=3d-5
        node(i)%khold=khold
        node(i)%diffe=0d0
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=1d0 ; node(i)%adh=0d0
        node(i)%rep=1d0 ; node(i)%repcel=5d0
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.30
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=1d0  !only for epithelium
        !node(i)%stor=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=desmax
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if


!    do i=2,7
!      a=node(1)%x-node(i)%x ; b=node(1)%y-node(i)%y ; c=node(1)%z-node(i)%z
!      d=sqrt(a**2+b**2+c**2) ; a=a*0.15/d ; b=b*0.15/d ; c=c*0.15/d
!      node(i)%x=node(i)%x+a ; node(i)%y=node(i)%y+b ; node(i)%z=node(i)%z+c
!    end do

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%hold=0
!    node(2:7)%hold=1 !; node(2:7)%rep=1d3 ; node(2:7)%you=1d3
    !let's do it for the most external layer of cells
      j=0
      do i=1,radicel-2
        j=j+i
      end do
      j=(6*j+1)*nodecel !this is the number of epithelial nodes wich are not external
      do i=j+1,ndepi
        node(i)%hold=1 ;node(i)%khold=khold
        node(i)%oriz=node(i)%oriz-20
        if(node(i)%tipus==2)then;node(i)%orix=0d0 ; node(i)%oriy=0d0;end if
        !node(i)%rep=1d1;node(i)%repcel=1d1
        !node(i)%ke=1d1
        !node(i)%tor=1d1
        !!node(i)%stor=1d1
      end do
      j=0
      do i=1,mradicel-2
        j=j+i
      end do
      !j=(6*j+1)*mradi !this is the number of mesenchymal nodes wich are not external
      !do i=ndepi+j+1,ndepi+ndmes
      !  node(i)%hold=1
      !  !node(i)%rep=1d1;node(i)%repcel=1d1
      !end do
    
    
!    do i=ndepi+1,nd ; node(i)%hold=1;enddo
    
    
    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=5
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=1d1 ; gen(1)%mu=0d0 !activator transcript
      
      gen(2)%kindof=1 ; gen(2)%diffu=1d1 ; gen(2)%mu=0d0 !epithelial basal adhesion molecule
      gen(3)%kindof=1 ; gen(3)%diffu=1d1 ; gen(3)%mu=0d0 !epithelial apical adhesion molecule
      gen(4)%kindof=1 ; gen(4)%diffu=1d1 ; gen(4)%mu=0d0 !mesenchymal adhesion molecule

      gen(5)%kindof=1 ; gen(5)%diffu=1d1 ; gen(5)%mu=0d0 !activator transcript



    !Gene-behavior interactions
 
       gen(2)%wa(1)=1  !epithelial adhesion molecule
       gen(3)%wa(1)=2  !epithelial-mesenchymal adhesion molecule
       gen(4)%wa(1)=3 !mesenchymal adhesion molecule

    !Gene-gene interactions



    !Adhesion molecules

	ntipusadh=3
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=5d0
      kadh(1,2)=5d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=5d0
      kadh(1,3)=5d0 ; kadh(3,1)=kadh(1,3)
      kadh(2,2)=5d0
      kadh(2,3)=0d0 ; kadh(3,2)=kadh(2,3)
      kadh(3,3)=1d1


    end if

    !Gene expression on nodes

    do i=1,nd
      !gex(i,3)=1d3 !;gex(i,10)=1d-5
      if(node(i)%tipus==3)then
        !gex(i,5)=1d3
        !gex(i,7)=1d3
        gex(i,4)=1d0 ;gex(i,5)=1d0
        !if(node(i)%marge==0)then
        !  gex(i,1)=1d0
        !end if
         gex(i,1)=1d0
      else
        !gex(i,9)=1d3
        gex(i,1)=1d0
        if(node(i)%tipus==2)then
          gex(i,2)=1d0
        else
          gex(i,3)=1d0
        end if
      end if
    end do

    call update_npag

    node(:)%talone=0.0d0
    ramax=maxval(node(:)%da)*3

    
    
end subroutine


!!!!!!!!!***********************************SUBROUTINE********************************

subroutine hair_placode



!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=7    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=7   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=0.75d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    mesradicel(1)=mradicel
    !mesradicel(2)=mradicel
    !mesradicel(3)=2

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      ncelsmes=0
      do k=1,layer
        j=0 ;print*,"mesradicel",mesradicel(k)
        do i=1,mesradicel(k)-1
          j=j+i
        end do
        ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
      end do
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
    allocate(ffu(nfu))
    if(radi==1.or.mradi==1) ffu(1)=1
  !End initializing dimension parameters

  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.5d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=0.7d0
    khold=1d0

    k_press=0d0
    m_xwall=zmes !tooth
    mi_xwall=0d0 !tooth


    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=1
    ffu(2)=0 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=3 !external signal source
    ffu(6)=0 !eggshell miguel4-1-13

    ffu(11)=0 !epithelial node plastic deformation
    ffu(12)=0 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(14)=1 !physical boundaries (walls)
    ffu(22)=1 !0 = unbiased random noise / 1 = noise biased by energies

   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(i<=ndepi/2)then  !basal
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d0
        else                      !apical
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d0
        end if
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.125d0
        node(i)%da=node(i)%req*2.0; 
        node(i)%ke=1d1
        node(i)%tor=5d0
        node(i)%stor=5d0
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=1d1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.70
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=1d0  !only for epithelium
        !node(i)%stor=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if



    do i=1,nd
      if(node(i)%tipus==1) node(i)%marge=0
    end do

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
    do i=1,nd
      call random_number(a)
      node(i)%x=node(i)%x+2*a*desmax-desmax
      call random_number(a)
      node(i)%y=node(i)%y+2*a*desmax-desmax
      call random_number(a)
      node(i)%z=node(i)%z+2*a*desmax-desmax
    end do



    
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%hold=0
    !!let's do it for the most external layer of cells
      j=0
      do i=1,radicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
      do i=j+1,ndepi
        node(i)%hold=0
    !    !node(i)%da=node(i)%req*2.0
    !    !node(i)%orix=0 ; node(i)%oriy=0
    !    !node(i)%oriz=node(i)%oriz-100
    !    !node(i)%rep=1d1;node(i)%repcel=1d1
    !    !node(i)%ke=1d1
    !    !node(i)%tor=1d1
    !    !!node(i)%stor=1d1
      end do
      j=0
      do i=1,mradicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      k=0
      do i=1,mradicel-1
        k=k+i
      end do
      k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
        node(i)%hold=0
        !node(i)%da=node(i)%req*2.0
        !node(i)%orix=0 ; node(i)%oriy=0
        !node(i)%rep=1d1;node(i)%repcel=1d1
      end do
      !do i=ndepi+k+j+1,ndepi+2*k !middle mesenchymal layer
      !  node(i)%hold=1
        !node(i)%da=node(i)%req*2.0
        !node(i)%orix=0 ; node(i)%oriy=0
        !node(i)%rep=1d1;node(i)%repcel=1d1
      !end do
      !do i=ndepi+2*k+1,nd  !lower mesenchymal layer
      !  node(i)%hold=1
      !  !node(i)%da=node(i)%req*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%repcel=1d1
      !end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels
      call random_number(a)
      cels(i)%fase=a
    end do



   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=6
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=0d0 ; gen(1)%mu=5d-1 ;gen(1)%name="interplacodal marker"!E-cadherin
      
      gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=5d-1 ;gen(2)%name="placodal marker"!P-cadherin

      gen(3)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=0d0 !basal lamina

      gen(4)%kindof=1 ; gen(4)%diffu=0d0 ; gen(4)%mu=5d-1 ;gen(4)%name="migration gene"

      gen(5)%kindof=2 ; gen(5)%diffu=0d0 ; gen(5)%mu=0d0 !morphogen transcript
      !gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost))
      !gen(5)%post(1)=6

      gen(6)%kindof=4 ; gen(6)%diffu=0d0 ; gen(6)%mu=0d0 ;gen(6)%name="morphogen gradient"
      !gen(6)%npre=1 ; allocate(gen(6)%pre(gen(6)%npre))
      !gen(6)%pre(1)=5

      
    !Gene-behavior interactions
    

       gen(1)%wa(1)=1  !interplacodal marker
       gen(2)%wa(1)=2  !placodal marker
       gen(3)%wa(1)=3  !basal lamina
       !gen(4)%wa(5)=-0.05 !this is req (emulating a migratory behavior)
       !gen(4)%wa(6)=0.10 !this is da (emulating a migratory behavior)
       gen(4)%wa(16)=0.01 !this is dmo (migratory cells)
       
       gen(1)%wa(nparam_per_node+2)=1d-1
       gen(2)%wa(nparam_per_node+2)=2d-1
       gen(3)%wa(nparam_per_node+2)=1d-1

       gen(6)%wa(nparam_per_node+8)=1d0
       gen(6)%wa(nparam_per_node+16)=5d-1 !this makes random noise biased towards the gradient of the gene



    !Gene-gene interactions

      gen(1)%w(1)=1d0
      gen(2)%w(2)=1d0
      gen(4)%w(4)=1d0
      !gen(5)%w(5)=1.0d3
      !gen(5)%nww=1
      !gen(5)%ww(1,1)=5
      !gen(5)%ww(1,2)=6
      !gen(5)%ww(1,3)=1d2

    !Adhesion molecules

	ntipusadh=3
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=1d0
      kadh(1,2)=1d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=1d0
      kadh(1,3)=1d0 ; kadh(3,1)=kadh(1,3)
      kadh(2,3)=1d0 ; kadh(3,2)=kadh(2,3)
      kadh(3,3)=1d0
!
    end if

    !Gene expression on nodes

    j=0
    do ii=1,radicel-3
      j=j+ii
    end do
    j=(6*j+1) ;print*,"jota",j

    do i=1,nd
      !a=sqrt(node(i)%x**2+node(i)%y**2)
      if(i<=ndepi)then
        if(node(i)%icel<=j)then
          !gex(i,2)=1d0!(3.5-a)/3.5 !; print*,"gex(i,2)",gex(i,2)
          gex(i,3)=1d0
          !gex(i,4)=1d0
          !gex(i,4)=1/(1+sqrt(node(i)%x**2+node(i)%y**2))
        else
          gex(i,1)=1d0
        end if
        !if(node(i)%tipus==2)then
        !  gex(i,3)=1d0
        !  gex(i,1)=0d0 ; gex(i,2)=0d0
        !end if
      else
        j=0
        do ii=1,mradicel-3
          j=j+ii
        end do
        j=(6*j+1) !territory 1 (inner ring)
        k=0
        do ii=1,mradicel-1
          k=k+ii
        end do
        k=(6*k+1) !whole layer
        if(node(i)%icel<=ncelsepi+j)then !upper layer - inner ring
          gex(i,2)=1d0!(3.5-a)/3.5
          gex(i,4)=1d0
          !gex(i,4)=1/(1+sqrt(node(i)%x**2+node(i)%y**2))
        elseif(node(i)%icel>ncelsepi+j.and.node(i)%icel<=ncelsepi+k)then !upper layer - outer ring
          gex(i,1)=1d0
        elseif(node(i)%icel>ncelsepi+k .and. node(i)%icel<=ncelsepi+k+j)then !midle layer - inner ring
          gex(i,2)=1d0
          gex(i,4)=1d0
          !gex(i,4)=1/(1+sqrt(node(i)%x**2+node(i)%y**2))
        elseif(node(i)%icel>ncelsepi+k+j.and.node(i)%icel<=ncelsepi+2*k)then !midle layer - outer ring
          gex(i,1)=1d0
        !else           !bottom layer
        !  gex(i,5)=1d0
        !  gex(i,3)=1d0
        end if
      end if
    end do

    call update_npag

    node(:)%talone=0.0d0

    ramax=maxval(node(:)%da)*3
   
end subroutine

!************************************************************************************


subroutine feather_placode



!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=5    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=5   !number of radial cell layers
    layer=2      !number of planar cell layers
    zmes=-0.5d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    mesradicel(1)=mradicel
    mesradicel(2)=mradicel
    !mesradicel(3)=2

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      ncelsmes=0
      do k=1,layer
        j=0 !;print*,"mesradicel",mesradicel(k)
        do i=1,mesradicel(k)-1
          j=j+i
        end do
        ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
      end do
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    

  !End initializing dimension parameters

  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=1
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    deltamin=1d-3
    khold=1d0
    angletor=0d0

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=1
    ffu(2)=0 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !eggshell miguel4-1-13

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(14)=0 !physical boundaries (walls)
    ffu(17)=1 !volume conservation (epithelial)
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(i<=ndepi/2)then  !basal
          node(i)%you=5d0 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%repcel=1d1
        else                      !apical
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%repcel=1d1
        end if
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.25
        node(i)%da=node(i)%req*1.60; 
        node(i)%ke=5d1
        node(i)%tor=5d0
        node(i)%stor=5d0
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=5d-2
        node(i)%kvol=5d-2
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=5d0 ; node(i)%adh=0d0
        node(i)%rep=1d1 ; node(i)%repcel=3d1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.30
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=1d0  !only for epithelium
        node(i)%stor=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
    do i=1,nd
      call random_number(a)
      node(i)%x=node(i)%x+2*a*desmax-desmax
      call random_number(a)
      node(i)%y=node(i)%y+2*a*desmax-desmax
      call random_number(a)
      node(i)%z=node(i)%z+2*a*desmax-desmax
    end do
    
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%hold=0
    !!let's do it for the most external layer of cells
      j=0
      do i=1,radicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
      do i=j+1,ndepi
        node(i)%hold=1 ;node(i)%repcel=1d2
      end do
      j=0
      do i=1,mradicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      k=0
      do i=1,mradicel-1
        k=k+i
      end do
      k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
        node(i)%hold=1;node(i)%repcel=1d2
        node(i)%border=1
        !node(i)%da=node(i)%req*2.0
        !node(i)%orix=0 ; node(i)%oriy=0
        !node(i)%rep=1d1;node(i)%repcel=1d1
      end do
      do i=ndepi+k+1,nd
        node(i)%hold=1;node(i)%repcel=1d2
        node(i)%border=1
      end do
      !do i=ndepi+k+j+1,ndepi+2*k !middle mesenchymal layer
      !  node(i)%hold=1
      !  !node(i)%da=node(i)%req*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%repcel=1d1
      !end do
      !do i=ndepi+2*k+1,nd  !lower mesenchymal layer
      !  node(i)%hold=1
      !  !node(i)%da=node(i)%req*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%repcel=1d1
      !end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=14
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=1d1 ; gen(1)%mu=1d0 !Epithelial-cadherin
      
      gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=1d0 !Mesenchymal-cadherin

      gen(3)%kindof=1 ; gen(3)%diffu=1d0 ; gen(3)%mu=1d0 !housekeeping gene epithelial

      gen(4)%kindof=2 ; gen(4)%diffu=0d0 ; gen(4)%mu=0d0 !epi. morphogen transcript (like SHH)
      gen(4)%npost=1 ; allocate(gen(4)%post(gen(4)%npost))
      gen(4)%post(1)=5

      gen(5)%kindof=4 ; gen(5)%diffu=1d0 ; gen(5)%mu=1d0 !epi. morphogen diffusible form (like SHH)
      gen(5)%npre=2 ; allocate(gen(5)%pre(gen(5)%npre))
      gen(5)%pre(1)=4 ; gen(5)%pre(2)=9
      gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost))
      gen(5)%post(1)=9

      gen(6)%kindof=2 ; gen(6)%diffu=0d0 ; gen(6)%mu=0d0 !mes. morphogen transcript (like BMP)
      gen(6)%npost=1 ; allocate(gen(6)%post(gen(6)%npost))
      gen(6)%post(1)=7

      gen(7)%kindof=4 ; gen(7)%diffu=1d0 ; gen(7)%mu=5d-1 !mes. morphogen diffusible form (like BMP)
      gen(7)%npre=2 ; allocate(gen(7)%pre(gen(7)%npre))
      gen(7)%pre(1)=6 ; gen(7)%pre(2)=11
      gen(7)%npost=1 ; allocate(gen(7)%post(gen(7)%npost))
      gen(7)%post(1)=11
      
      gen(8)%kindof=2 ; gen(8)%diffu=0d0 ; gen(8)%mu=5d0 !SHH receptor
      gen(8)%npost=1 ; allocate(gen(8)%post(gen(8)%npost))
      gen(8)%post(1)=9
      
      gen(9)%kindof=8 ; gen(9)%diffu=0d0 ; gen(9)%mu=0d0 !SHH receptor activated
      gen(9)%npre=2 ; allocate(gen(9)%pre(gen(9)%npre))
      gen(9)%pre(1)=8 ; gen(9)%pre(2)=5
      gen(9)%npost=2 ; allocate(gen(9)%post(gen(9)%npost))
      gen(9)%post(1)=8 ; gen(9)%post(2)=5
      
      gen(10)%kindof=2 ; gen(10)%diffu=0d0 ; gen(10)%mu=5d0 !BMP receptor
      gen(10)%npost=1 ; allocate(gen(10)%post(gen(10)%npost))
      gen(10)%post(1)=11

      gen(11)%kindof=8 ; gen(11)%diffu=0d0 ; gen(11)%mu=0d0 !BMP receptor activated
      gen(11)%npre=2 ; allocate(gen(11)%pre(gen(11)%npre))
      gen(11)%pre(1)=10 ; gen(11)%pre(2)=5
      gen(11)%npost=2 ; allocate(gen(11)%post(gen(11)%npost))
      gen(11)%post(1)=10 ; gen(11)%post(2)=5

      gen(12)%kindof=1 ; gen(12)%diffu=0d0 ; gen(12)%mu=1d0 !housekeeping gene mesenchymal
      gen(13)%kindof=1 ; gen(13)%diffu=1d0 ; gen(13)%mu=1d0 !housekeeping gene signalling center
      gen(14)%kindof=1 ; gen(14)%diffu=0d0 ; gen(14)%mu=0d0 !adhesion molecule hold
      
      
    !Gene-behavior interactions
    

       gen(1)%wa(1)=1  !epithelial adhesion molecule
       gen(2)%wa(1)=2  !epithelial-mesenchymal adhesion molecule
       gen(14)%wa(1)=3  !basal lamina
       !gen(4)%wa(5)=-0.05 !this is req (emulating a migratory behavior)
       !gen(4)%wa(6)=0.10 !this is da (emulating a migratory behavior)
       !gen(4)%wa(16)=0.03 !this is dmo (migratory cells)
       !gen(6)%wa(nparam_per_node+8)=1d0
       !gen(6)%wa(nparam_per_node+16)=1d-2 !this makes random noise biased towards the gradient of the gene
       gen(9)%wa(nparam_per_node+2)=1d-2 !SHH effect on epithelial growth
       gen(11)%wa(nparam_per_node+2)=5d-3 !BMP effect on mesenchymal growth



    !Gene-gene interactions

      !gen(11)%nww=2      !BMP activated receptor induces production of BMP
      !gen(11)%ww(1,1)=6  
      !gen(11)%ww(1,2)=7
      !gen(11)%ww(1,3)=1d1
      !gen(11)%ww(2,1)=4   !BMP activated receptor induces production of SHH
      !gen(11)%ww(2,2)=5
      !gen(11)%ww(2,3)=1d1
      
      !gen(7)%nww=1      !BMP morphogen activates BMP receptor
      !gen(7)%ww(1,1)=10  
      !gen(7)%ww(1,2)=11
      !gen(7)%ww(1,3)=1d0

      gen(3)%w(3)=1d0 !housekeeping activates itself
      gen(1)%w(3)=1d0  !activates epi cadherin
      gen(8)%w(3)=1d0  !activates ssh receptor

      gen(12)%w(12)=1d0 !housekeeping activates itself
      gen(2)%w(12)=1d0   !activates mes cadherin
      gen(10)%w(12)=1d0  !activates BMP receptor

      gen(13)%w(13)=1d0 !housekeeping activates itself
      gen(1)%w(13)=1d0  !activates epi cadherin
      gen(4)%w(13)=1d1  !activates ssh transcript



      gen(13)%nww=1    !housekeeping epi mediates secretion of shh
      gen(13)%ww(1,1)=4
      gen(13)%ww(1,2)=5
      gen(13)%ww(1,3)=1d1

      
      gen(9)%nww=4      !SHH morphogen activates SHH receptor
      gen(9)%ww(1,1)=5  
      gen(9)%ww(1,2)=9
      gen(9)%ww(1,3)=1d0
      gen(9)%ww(2,1)=9  
      gen(9)%ww(2,2)=5
      gen(9)%ww(2,3)=1d0

      gen(9)%ww(3,1)=8  
      gen(9)%ww(3,2)=9
      gen(9)%ww(3,3)=1d0
      gen(9)%ww(4,1)=9  
      gen(9)%ww(4,2)=8
      gen(9)%ww(4,3)=1d0


      gen(11)%nww=4      !SHH morphogen activates SHH receptor
      gen(11)%ww(1,1)=5  
      gen(11)%ww(1,2)=11
      gen(11)%ww(1,3)=1d0
      gen(11)%ww(2,1)=11  
      gen(11)%ww(2,2)=5
      gen(11)%ww(2,3)=1d0

      gen(11)%ww(3,1)=10
      gen(11)%ww(3,2)=11
      gen(11)%ww(3,3)=1d0
      gen(11)%ww(4,1)=11  
      gen(11)%ww(4,2)=10
      gen(11)%ww(4,3)=1d0



                        !SHH morphogen activates BMP receptor
      !gen(5)%ww(3,1)=5  
      !gen(5)%ww(3,2)=11
      !gen(5)%ww(3,3)=1d0
      !gen(5)%ww(4,1)=11  
      !gen(5)%ww(4,2)=5
      !gen(5)%ww(4,3)=1d0

      !gen(8)%nww=2      !ssh receptor activates SHH receptor
      !gen(8)%ww(1,1)=8  
      !gen(8)%ww(1,2)=9
      !gen(8)%ww(1,3)=1d0
      !gen(8)%ww(2,1)=9  
      !gen(8)%ww(2,2)=8
      !gen(8)%ww(2,3)=1d0

      !gen(10)%nww=2               !bmp receptor activates BMP receptor
      !gen(10)%ww(1,1)=10  
      !gen(10)%ww(1,2)=11
      !gen(10)%ww(1,3)=1d0
      !gen(10)%ww(2,1)=11  
      !gen(10)%ww(2,2)=10
      !gen(10)%ww(2,3)=1d0



    !Adhesion molecules

	ntipusadh=3
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=2d1
      kadh(1,2)=1d1 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=2d1
      kadh(1,3)=1d0 ; kadh(3,1)=kadh(1,3)
      kadh(2,3)=1d0 ; kadh(3,2)=kadh(2,3)
      kadh(3,3)=1d2

    end if

    !Gene expression on nodes

    j=0
    do ii=1,radicel-2
      j=j+ii
    end do
    j=(6*j+1) !;print*,"jota",j

    do i=1,nd
      if(node(i)%hold==1) gex(i,14)=1d0
      !a=sqrt(node(i)%x**2+node(i)%y**2)
      if(i<=2)then; gex(i,4)=1d0; gex(i,13)=1;end if
      if(i<=ndepi)then
        gex(i,1)=1d0
        if(node(i)%icel>1)then
          gex(i,3)=1d0
          gex(i,8)=1d0
        !else
        !  gex(i,1)=1d0
        end if
      else
        gex(i,12)=1d0
        gex(i,10)=1d0
        j=0
        do ii=1,mradicel-2
          j=j+ii
        end do
        j=(6*j+1) !territory 1 (inner ring)
        k=0
        do ii=1,mradicel-1
          k=k+ii
        end do
        k=(6*k+1) !whole layer
        if(node(i)%icel<=ncelsepi+j)then !upper layer - inner ring
          gex(i,2)=1d0!(3.5-a)/3.5
          !gex(i,6)=1d0
          !gex(i,7)=1d0
          !gex(i,10)=1d0
          !gex(i,4)=1/(1+sqrt(node(i)%x**2+node(i)%y**2))
        elseif(node(i)%icel>ncelsepi+j.and.node(i)%icel<=ncelsepi+k)then !upper layer - outer ring
          !gex(i,1)=1d0
        elseif(node(i)%icel>ncelsepi+k .and. node(i)%icel<=ncelsepi+k+j)then !midle layer - inner ring
          gex(i,2)=1d0
          !gex(i,6)=1d0
          !gex(i,7)=1d0
          !gex(i,10)=1d0
          !gex(i,4)=1/(1+sqrt(node(i)%x**2+node(i)%y**2))
        elseif(node(i)%icel>ncelsepi+k+j.and.node(i)%icel<=ncelsepi+2*k)then !midle layer - outer ring
          !gex(i,1)=1d0
        !else           !bottom layer
        !  gex(i,5)=1d0
        !  gex(i,3)=1d0
        end if
      end if
    end do

    


    call update_npag

    node(:)%talone=0.0d0

    ramax=maxval(node(:)%da)*3

    
end subroutine

!************************************************************************************



subroutine epi_mes_bud_ingrowth



!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=5    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=5   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=1.0d0  !z-position of uppermost layer

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+30
    
    if(radi==1.or.mradi==1) ffu(1)=1
  !End initializing dimension parameters

  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-4
    prop_noise=0.3d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    khold=5d0
    mnn=700

    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(2)=0 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !eggshell miguel4-1-13

    ffu(11)=0 !epithelial node plastic deformation
    ffu(12)=0 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(22)=0 !0 = noise biased by energies , 1 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(node(i)%tipus==2)then  !basal
          node(i)%you=2d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%repcel=1d1
          node(i)%tor=3d1
          !node(i)%stor=1d1
        else                      !apical
          node(i)%you=2d1 ; node(i)%adh=0d0
          node(i)%rep=2d1 ; node(i)%repcel=2d1
          node(i)%tor=1d1
          !node(i)%stor=1d1
        end if
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.25
        node(i)%da=node(i)%req*1.40; 
        node(i)%ke=1d1
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=5d-2
        node(i)%kvol=1d-2
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=2d1 ; node(i)%adh=0d0
        node(i)%rep=1d1 ; node(i)%repcel=1d1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.60
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=1d0  !only for epithelium
        !node(i)%stor=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if



    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%hold=0
    !let's do it for the most external layer of cells
      j=0
      do i=1,radicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
      do i=j+1,ndepi
        node(i)%hold=1
        !node(i)%da=node(i)%req*2.0
        !node(i)%orix=0 ; node(i)%oriy=0
        !node(i)%oriz=node(i)%oriz-100
        !node(i)%rep=1d2;node(i)%repcel=1d2
        !node(i)%ke=1d1
        !node(i)%tor=1d1
        !!node(i)%stor=1d1
      end do
      j=0
      do i=1,mradicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      do i=ndepi+j+1,ndepi+ndmes
        node(i)%hold=1
        !node(i)%da=node(i)%req*2.0
        !node(i)%orix=0 ; node(i)%oriy=0
        !node(i)%rep=1d2;node(i)%repcel=1d2
      end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=4
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=0d0 ; gen(1)%mu=0d0 !E-cadherin
      
      gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=0d0 !P-cadherin

      gen(3)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=0d0 !basal lamina

      gen(4)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=0d0 !motility gene

      !gen(5)%kindof=1 ; gen(5)%diffu=0d0 ; gen(5)%mu=0d0 !growth gene epithelium

    !Gene-behavior interactions
    

       gen(1)%wa(1)=1  !interplacodal adhesion molecule
       gen(2)%wa(1)=2  !placodal adhesion molecule
       gen(3)%wa(1)=3  !basal lamina
       !gen(4)%wa(nparam_per_node+1)=4d-4 !growth mesenchyme
       !gen(4)%wa(nparam_per_node+2)=1d-3 !cell cycle
       !gen(5)%wa(nparam_per_node+1)=1d-5 !growth epithelium
       !gen(5)%wa(nparam_per_node+2)=1d-3 !cell cycle

       gen(4)%wa(5)=0.00 !this is req (emulating a migratory behavior)
       gen(4)%wa(6)=0.15 !this is da (emulating a migratory behavior)
       !gen(4)%wa(9)=1d1 ; gen(4)%wa(10)=1d1 !this is rep (so the nodes don't collapse)
       gen(4)%wa(16)=0.01 !this is dmo (migratory cells)       
       

    !Gene-gene interactions



    !Adhesion molecules

	ntipusadh=3
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=8d1
      kadh(1,2)=5d1 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=8d1
      kadh(1,3)=0d0 ; kadh(3,1)=kadh(1,3)
      kadh(2,3)=0d0 ; kadh(3,2)=kadh(2,3)
      kadh(3,3)=2d2


    end if

    !Gene expression on nodes

    j=0
    do ii=1,radicel-2
      j=j+ii
    end do
    j=(6*j+1) ;print*,"jota",j

    do i=1,nd
      !gex(i,1)=1d0
      if(node(i)%tipus==2) gex(i,3)=1d0
      !if(i<=ndepi)then
        if(node(i)%icel<=j.or.(i>ndepi.and.node(i)%icel<=ncelsepi+j))then
          gex(i,4)=1d0
          gex(i,2)=1d0 !P-cadherin
        !if(node(i)%tipus==2)then
          !gex(i,3)=1d0
        else
          gex(i,1)=1d0 !E-cadherin
        end if
      !else
        !gex(i,1)=1d0
        !if(node(i)%icel<=ncelsepi+j) gex(i,4)=1d0
      !end if

      !if(i<=ndepi)then
      !  if(node(i)%icel<=61)then
      !    gex(i,2)=1d0!(3.5-a)/3.5 !; print*,"gex(i,2)",gex(i,2)
      !    gex(i,4)=1d0
      !  else
      !    gex(i,1)=1d0
      !  end if
      !  !if(node(i)%tipus==2)then
      !  !  gex(i,3)=1d0
      !  !  gex(i,1)=0d0 ; gex(i,2)=0d0
      !  !end if
      !else
      !  j=0
      !  do ii=1,mradicel-2
      !    j=j+ii
      !  end do
      !  j=(6*j+1) !;print*,"jota",j
      !  if(node(i)%icel<=ncelsepi+j)then
      !    gex(i,2)=1d0!(3.5-a)/3.5
      !    gex(i,4)=1d0
      !  elseif(node(i)%icel>ncelsepi+ncelsmes/2 .and.node(i)%icel<=ncelsepi+ncelsmes/2+j)then
      !    gex(i,2)=1d0!(3.5-a)/3.5
      !    gex(i,4)=1d0
      !  else
      !    gex(i,1)=1d0
      !  end if
      !end if
    end do

    call update_npag

    node(:)%talone=0.0d0

    ramax=maxval(node(:)%da)*3

    do i=1,nd
      call random_number(a)
      node(i)%x=node(i)%x+2*a*desmax-desmax
      call random_number(a)
      node(i)%y=node(i)%y+2*a*desmax-desmax
      call random_number(a)
      node(i)%z=node(i)%z+2*a*desmax-desmax
    end do
    
end subroutine


!****************************************************************************************************



!**************************************************************************************************

subroutine tooth_bud



!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=5    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=5   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=1d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    mesradicel(1)=mradicel
    !mesradicel(2)=mradicel
    !mesradicel(3)=2

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      ncelsmes=0
      do k=1,layer
        j=0 !;print*,"mesradicel",mesradicel(k)
        do i=1,mesradicel(k)-1
          j=j+i
        end do
        ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
      end do
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
	nd=ndepi+ndmes
	ncels=ncelsepi+ncelsmes
    nda=nd+10
	ncals=ncels+10
    
  !End initializing dimension parameters

  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d-1 !low value is low temperature	
    !desmax=0.001
    resmax=1d-3
    prop_noise=0.1d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    khold=5d0

    !biological
    !functions used
    
    
    ffu=0
    if(radi==1.or.mradi==1) ffu(1)=1

     !spring of the ellipse
    ffu(1)=1
    ffu(2)=0 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !eggshell miguel4-1-13

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=0 !0 dynamic delta, 1 fixed delta
    ffu(13)=1 !neighboring by triangulation
    ffu(14)=1 !physical boundaries (walls)
    ffu(15)=1 !fixed gene 1 on borders
    ffu(22)=0 !0 = noise biased by energies , 1 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(i<=ndepi/2)then  !basal
          node(i)%you=5d0 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%repcel=1d1
        else                      !apical
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%repcel=1d1
        end if
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.25
        node(i)%da=node(i)%req*1.60; 
        node(i)%ke=5d1
        node(i)%tor=5d0
        !node(i)%stor=5d0
        node(i)%mo=temp
        node(i)%dmo=desmax
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=1d-2
        node(i)%kvol=1d-1
      end do
    end if

	!Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=5d0 ; node(i)%adh=0d0
        node(i)%rep=1d1 ; node(i)%repcel=3d1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.30
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=1d0  !only for epithelium
        !node(i)%stor=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=desmax
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
    !do i=1,nd
    !  call random_number(a)
    !  node(i)%x=node(i)%x+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%y=node(i)%y+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%z=node(i)%z+2*a*desmax-desmax
    !end do
    
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%hold=0
    !!let's do it for the most external layer of cells
      j=0
      do i=1,radicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
      do i=j+1,ndepi
        node(i)%hold=1 ;node(i)%repcel=5d1
      end do
      j=0
      do i=1,mradicel-2
        j=j+i
      end do
      j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      k=0
      do i=1,mradicel-1
        k=k+i
      end do
      k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
        node(i)%hold=1;node(i)%repcel=5d1
        node(i)%da=node(i)%req*1.50
        !node(i)%orix=0 ; node(i)%oriy=0
        !node(i)%rep=1d1;node(i)%repcel=1d1
      end do
      !do i=ndepi+k+1,nd !lower mesenchymal layer all
      !  node(i)%hold=1;node(i)%repcel=5d1
      !  node(i)%da=node(i)%req*1.50
      !end do
      !do i=ndepi+k+j+1,ndepi+2*k !middle mesenchymal layer
      !  node(i)%hold=1
      !  !node(i)%da=node(i)%req*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%repcel=1d1
      !end do
      !do i=ndepi+2*k+1,nd  !lower mesenchymal layer
      !  node(i)%hold=1
      !  !node(i)%da=node(i)%req*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%repcel=1d1
      !end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=11
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=4 ; gen(1)%diffu=1d-1 ; gen(1)%mu=1d0 !Wnt7b ligand, comes from the epithelial borders
      gen(1)%npre=1 ; allocate(gen(1)%pre(gen(1)%npre)) ; gen(1)%pre(1)=2

      gen(2)%kindof=2 ; gen(2)%diffu=1d0 ; gen(2)%mu=5d-1 !wnt7b transcript
      gen(2)%npost=1 ; allocate(gen(2)%post(gen(2)%npost)) ; gen(2)%post(1)=2
      
      gen(3)%kindof=1 ; gen(3)%diffu=1d0 ; gen(3)%mu=5d-1 !Housekeeping gene (it keeps itself at stable levels, may activate other genes)

      gen(4)%kindof=1 ; gen(4)%diffu=1d0 ; gen(4)%mu=5d-1 !Adhesion molecule epithelial (activated by housekeeping)

      gen(5)%kindof=4 ; gen(5)%diffu=1d-1 ; gen(5)%mu=1d0 !Wnt10b ligand
      gen(5)%npre=1 ; allocate(gen(5)%pre(gen(5)%npre)) ; gen(5)%pre(1)=6
      
      gen(6)%kindof=2 ; gen(6)%diffu=1d0 ; gen(6)%mu=5d-1 !Wnt10b transcript
      gen(6)%npost=1 ; allocate(gen(6)%post(gen(6)%npost)) ; gen(6)%post(1)=5
      
      gen(7)%kindof=1 ; gen(7)%diffu=1d0 ; gen(7)%mu=0d0 !Adhesion molecule external

      gen(8)%kindof=2 ; gen(8)%diffu=1d0 ; gen(8)%mu=5d-1 !Wnt receptor inactive
      gen(8)%npost=2 ; allocate(gen(8)%post(gen(8)%npost)) ; gen(8)%post(1)=9 ; gen(8)%post(2)=10
      
      gen(9)%kindof=3 ; gen(9)%diffu=1d0 ; gen(9)%mu=5d-1 !Wnt receptor activated by Wnt7b
      gen(9)%npre=1 ; allocate(gen(9)%pre(gen(9)%npre)) ; gen(9)%pre(1)=8

      gen(10)%kindof=3 ; gen(10)%diffu=1d0 ; gen(10)%mu=5d-1 !Wnt receptor activated by Wnt10b
      gen(10)%npre=1 ; allocate(gen(10)%pre(gen(10)%npre)) ; gen(10)%pre(1)=8
      
      gen(11)%kindof=1 ; gen(11)%diffu=1d0 ; gen(11)%mu=1d-1 !intermediate mediator for Wnt7b (slow but resilient)
     

      
    !Gene-behavior interactions
    

       gen(4)%wa(1)=1  !Adhesion molecule epithelial
       gen(7)%wa(1)=2  !Adhesion molecule external

       gen(10)%wa(nparam_per_node+2)=4d-4 !Wnt10b promotes division



    !Gene-gene interactions

      gen(3)%w(3)=1d0  !housekeeping activates himself
    
      gen(4)%w(3)=1d0  !housekeeping produces adhesion molecule
      gen(8)%w(3)=1d0  !housekeeping produces Wnt receptor inactive

      gen(2)%w(10)=-4.0d0  !receptor activated by Wnt10b inhibits Wnt7b production
      gen(6)%w(10)=4.0d0  !receptor activated by Wnt10b promotes Wnt10b production

      gen(11)%w(9)=1d-1  !wnt7b activated receptor produces mediator at a slow rate

      gen(2)%w(11)=1d0  !mediator produced by Wnt7b promotes Wnt7b production
      gen(6)%w(11)=-1d0  !mediatior produced by Wnt7b inhibits Wnt10b production

      gen(3)%nww=2
      gen(3)%ww(1,1)=2  !housekeeping processes Wnt7b secretion
      gen(3)%ww(1,2)=1
      gen(3)%ww(1,3)=1d0
      gen(3)%ww(2,1)=6  !housekeeping processes Wnt10b secretion
      gen(3)%ww(2,2)=5
      gen(3)%ww(2,3)=1d0

      gen(1)%nww=1
      gen(1)%ww(1,1)=8  !Wnt7b activates Wnt receptor into specific activated form
      gen(1)%ww(1,2)=9
      gen(1)%ww(1,3)=1d1

      gen(5)%nww=1
      gen(5)%ww(1,1)=8  !Wnt10b activates Wnt receptor into specific activated form
      gen(5)%ww(1,2)=10
      gen(5)%ww(1,3)=1d1


    !Adhesion molecules

	ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=2d1
      kadh(1,2)=1d1 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=5d1  !external adhesion between hold nodes

    end if

    !Gene expression on nodes

    j=0
    do ii=1,radicel-2
      j=j+ii
    end do
    j=(6*j+1) !;print*,"jota",j

    do i=1,nd
      if(node(i)%hold==1)then
        gex(i,1)=1d0
        gex(i,7)=1d0
      elseif(i<=ndepi)then
        gex(i,3)=1d0 !housekeeping
        gex(i,4)=1d0 !adhesion molecule
        gex(i,5)=5d-2 !Wnt10b ligand
        gex(i,6)=5d-2 !Wnt10b transcript
        gex(i,8)=1d0 !receptor
        
      else
        gex(i,3)=1d0
        gex(i,4)=1d0
        gex(i,5)=5d-2 !Wnt10b ligand
        gex(i,6)=5d-2 !Wnt10b transcript
        gex(i,8)=1d0
        !j=0
        !do ii=1,mradicel-2
        !  j=j+ii
        !end do
        !j=(6*j+1) !territory 1 (inner ring)
        !k=0
        !do ii=1,mradicel-1
        !  k=k+ii
        !end do
        !k=(6*k+1) !whole layer
        !if(node(i)%icel<=ncelsepi+j)then !upper layer - inner ring
        !elseif(node(i)%icel>ncelsepi+j.and.node(i)%icel<=ncelsepi+k)then !upper layer - outer ring
        !elseif(node(i)%icel>ncelsepi+k .and. node(i)%icel<=ncelsepi+k+j)then !midle layer - inner ring
        !elseif(node(i)%icel>ncelsepi+k+j.and.node(i)%icel<=ncelsepi+2*k)then !midle layer - outer ring
        !end if
      end if
    end do

    call update_npag

    node(:)%talone=0.0d0

    ramax=maxval(node(:)%da)*3

    
end subroutine

subroutine delta_notch



!******* #1 DEFINING SPATIAL DIMENSIONS *******

    !epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=2    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=0   !number of radial cell layers
    layer=3      !number of planar cell layers
    zmes=-0.5d0  !z-position of uppermost layer
   
    allocate(mesradicel(layer))
    mesradicel(1)=mradicel
    mesradicel(2)=mradicel
    mesradicel(3)=mradicel

    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2    ;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1    !cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      ncelsmes=0
      do k=1,layer
        j=0 !;print*,"mesradicel",mesradicel(k)
        do i=1,mesradicel(k)-1
          j=j+i
        end do
        ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
      end do
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
    nd=ndepi+ndmes
    ncels=ncelsepi+ncelsmes
    nda=nd+10
    ncals=ncels+10    
    
  !End initializing dimension parameters

  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals))
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature   
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13  
    dmax=1
    screen_radius=1.0d0
    khold=1d0

    !biological
    !functions used
   
   
    ffu=0 ; if (radi==1.or.mradi==1) ffu(1)=1
    ffu(2)=0 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=1 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !eggshell miguel4-1-13

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=0 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(14)=0 !physical boundaries (walls)
    ffu(17)=1 !volume conservation (epithelial)
    ffu(22)=0 !0 = noise biased by energies , 1 = unbiased noise


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(i<=ndepi/2)then  !basal
          node(i)%you=5d0 ; node(i)%adh=1d1
          node(i)%rep=1d1 ; node(i)%repcel=1d1
        else                      !apical
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%repcel=1d1
        end if
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.25
        node(i)%da=node(i)%req*1.60;
        node(i)%ke=5d1
        node(i)%tor=5d0
        node(i)%stor=5d0
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=5d-2
        node(i)%kvol=5d-2
      end do
    end if

    !Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=5d0 ; node(i)%adh=0d0
        node(i)%rep=1d1 ; node(i)%repcel=3d1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.30
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=1d0  !only for epithelium
        node(i)%stor=1d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)


    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
   
    !do i=1,nd  !ACTIVATE THIS IF YOU WANT THAT THE INITIAL POSITIONS OF NODES ARE A BIT NOISY
    !  call random_number(a)
    !  node(i)%x=node(i)%x+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%y=node(i)%y+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%z=node(i)%z+2*a*desmax-desmax
    !end do
   
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%hold=0
    !!let's do it for the most external layer of cells
!      j=0
!      do i=1,radicel-2
!        j=j+i
!      end do
!      j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
!      do i=j+1,ndepi
!        node(i)%hold=1 ;node(i)%repcel=1d2
!      end do
!      j=0
!      do i=1,mradicel-2
!        j=j+i
!      end do
!      j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
!      k=0
!      do i=1,mradicel-1
!        k=k+i
!      end do
!      k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
!      do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
!        node(i)%hold=1;node(i)%repcel=1d2
        !node(i)%da=node(i)%req*2.0
        !node(i)%orix=0 ; node(i)%oriy=0
        !node(i)%rep=1d1;node(i)%repcel=1d1
!      end do
!      do i=ndepi+k+1,nd
!        node(i)%hold=1;node(i)%repcel=1d2
!      end do
      !do i=ndepi+k+j+1,ndepi+2*k !middle mesenchymal layer
      !  node(i)%hold=1
      !  !node(i)%da=node(i)%req*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%repcel=1d1
      !end do
      !do i=ndepi+2*k+1,nd  !lower mesenchymal layer
      !  node(i)%hold=1
      !  !node(i)%da=node(i)%req*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%repcel=1d1
      !end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
   
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=7
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=7 ; gen(1)%diffu=0d0 ; gen(1)%mu=0d0!5d-1 !activated notch
      gen(1)%npre=1 ; allocate(gen(1)%pre(gen(1)%npre)); gen(1)%pre(1)=4
      gen(1)%npost=1 ; allocate(gen(1)%post(gen(1)%npost)); gen(1)%post(1)=4

     
      gen(2)%kindof=7 ; gen(2)%diffu=0d0 ; gen(2)%mu=0d0!5d-1 !activated delta
      gen(2)%npre=1 ; allocate(gen(2)%pre(gen(2)%npre)); gen(2)%pre(1)=5
      gen(2)%npost=1 ; allocate(gen(2)%post(gen(2)%npost)); gen(2)%post(1)=5


      gen(3)%kindof=1 ; gen(3)%diffu=5d0 ; gen(3)%mu=3d-1 !Housekeeping gene epithelial

      gen(4)%kindof=2 ; gen(4)%diffu=0d0 ; gen(4)%mu=0d0!5d-1 !inactive notch
      gen(4)%npost=1 ; allocate(gen(4)%post(gen(4)%npost)); gen(4)%post(1)=1
      gen(4)%npre=1 ; allocate(gen(4)%pre(gen(4)%npre)); gen(4)%pre(1)=1

     
      gen(5)%kindof=2 ; gen(5)%diffu=0d0 ; gen(5)%mu=0d0!5d-1 !inactive delta
      gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost)); gen(5)%post(1)=2
      gen(5)%npre=1 ; allocate(gen(5)%pre(gen(5)%npre)); gen(5)%pre(1)=2

     
      gen(6)%kindof=2 ; gen(6)%diffu=5d0 ; gen(6)%mu=5d-1 !effector gene inactive
      gen(6)%npost=1 ; allocate(gen(6)%post(gen(6)%npost)); gen(6)%post(1)=7
     
      gen(7)%kindof=3 ; gen(7)%diffu=5d0 ; gen(7)%mu=5d-1 !effector gene active
      gen(7)%npre=1 ; allocate(gen(7)%pre(gen(7)%npre)); gen(7)%pre(1)=6
     
     
    !Gene-behavior interactions
   

       !gen(1)%wa(1)=1  !epithelial adhesion molecule
       !gen(2)%wa(1)=2  !mesenchymal adhesion molecule
       !gen(3)%wa(1)=3  !basal lamina


    !Gene-gene interactions

      !gen(3)%w(3)=1d0 !housekeeping epi. maintains its levels of expressions
      !gen(4)%w(3)=1d0
      !gen(5)%w(3)=1d0
      !gen(6)%w(3)=1d0

      gen(4)%nww=2    !inactive notch mediates binding of delta
      gen(4)%ww(1,1)=5
      gen(4)%ww(1,2)=2
      gen(4)%ww(1,3)=1d0
      gen(4)%ww(2,1)=2
      gen(4)%ww(2,2)=5
      gen(4)%ww(2,3)=1d0

      gen(5)%nww=2    !inactive delta mediates binding of notch
      gen(5)%ww(1,1)=4
      gen(5)%ww(1,2)=1
      gen(5)%ww(1,3)=1d0
      gen(5)%ww(2,1)=1
      gen(5)%ww(2,2)=4
      gen(5)%ww(2,3)=1d0

     
      !gen(1)%nww=2
      !gen(1)%ww(1,1)=6     !activated notch activates effector
      !gen(1)%ww(1,2)=7
      !gen(1)%ww(1,3)=1d-1
      !gen(1)%ww(2,1)=2     !activated notch mediates unbinding of delta
      !gen(1)%ww(2,2)=5
      !gen(1)%ww(2,3)=1d0

      !gen(2)%nww=1     !active delta mediates unbinding of notch
      !gen(2)%ww(1,1)=1
      !gen(2)%ww(1,2)=4
      !gen(2)%ww(1,3)=1d0    
     
     
      !gen(5)%w(7)=-1d1


    !Adhesion molecules

    ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      !kadh(1,1)=2d1
      !kadh(1,2)=1d1 ; kadh(2,1)=kadh(1,2)
      !kadh(2,2)=2d1
      !kadh(1,3)=1d0 ; kadh(3,1)=kadh(1,3)
      !kadh(2,3)=1d0 ; kadh(3,2)=kadh(2,3)
      !kadh(3,3)=1d2
!
    end if

    !Gene expression on nodes

    j=0
    do ii=1,radicel-2
      j=j+ii
    end do
    j=(6*j+1) !;print*,"jota",j

    do i=3,ndepi
      !if(node(i)%tipus==1)  gex(i,1)=1d0
      !if(node(i)%tipus==2)  gex(i,3)=1d0

      gex(i,4)=1d0
      !gex(i,5)=1d0
      !gex(i,3)=1d0
      !gex(i,6)=1d0
    end do

    gex(1:2,5)=1d0 ; !gex(1,5)=1d0
    !gex(3,5)=1d0 ; !gex(3,4)=1d0

    
    !do i=ndepi+1,nd
    !  gex(i,2)=1d0
    !  gex(i,5)=1d0
    !end do



    call update_npag

    node(:)%talone=0.0d0

    ramax=maxval(node(:)%da)*3
mmae=node(1)%req
   
end subroutine


!*****************************************************************
!**********************************************************************************************************************************************
subroutine tooth40a
real*8::pericel,radio2,req1,req2,adjustreq,di
!******* #1 DEFINING SPATIAL DIMENSIONS *******
   
    !epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=6    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer
   ! radius=0    ! no especificar
    !mesenchyme's dimension parameters
    packed=0     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1    !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=2   !number of radial cell layers
    layer=2      !number of planar cell layers
    zmes=1.30  !z-position of uppermost layer
    req=0.15
    !pericel=(pi*radius*radius)/(radicel*2-1)
    !req1=pericel
    !radio2=radius-pericel
   ! req2=(pi*radio2*radio2)/(radicel*2-1)
    radius=(req*(2*radicel-1)/pi) !external
    print*,"radius",radius
    radio2=radius-req! 
    print*,"radio2",radio2
    req2=((pi*radio2)/(2*radicel-1))!*0.9
    print*,"req2",req2
    adjustreq= req2/req
    radius=radius!*1.4
    print*,"radius*1.5",radius
    di=(req+req2)!*1.3

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2    ;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1    !cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
    nd=ndepi+ndmes
    ncels=ncelsepi+ncelsmes
print*,"ncels",ncels


    nda=nd+10
    ncals=ncels+10
print*,"ncals",ncals  
    single=0
    if(radi==1.or.mradi==1) single=1
  !End initializing dimension parameters
nodecel=2
  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals))
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=6
    idum=-8724683
    idumoriginal=idum
    realtime=0

   ! nparam=35       !number of params
    nvarglobal_out=5

    !physical
    temp=1d-2 !low value is low temperature   
    desmax=0.001
    resmax=1d-3
    nuns=100  !for screening
    angleto=90 ; angltormax=cos(angleto*2*pi/360) !maximum angle for epithelial torsion  
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13  
    dmax=2
    khold=1.0d0
    ecmmax=req
    !biological
    !functions used
   
   
    ffu=0
    ffu(1)=1 !spring of the ellipse
    ffu(2)=0 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !eggshell miguel4-1-13

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(17)=1 !conservacion del volumen
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise
    ffu(16)=1
   

!**************Distribution of nodes in space
  
    if(radi>0.and.radicel>0) call epiteli_sphere1(radi,radicel,radius,di,zepi)
    if(mradi>0.and.mradicel>0.and.layer>0) then
                       if(packed==0)then ; call mesenq1(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if
    print*,"node(ndepi+1)%req",node(ndepi+1)%req
    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do

    !do i=1,nd
    !  call random_number(a)
    !  node(i)%x=node(i)%x+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%y=node(i)%y+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%z=node(i)%z+2*a*desmax-desmax
    !end do
 
    !End distribution of nodes in space

   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
           if(node(i)%tipus==2)then  !basal                    !         if(i<=ndepi/2)then  !basal
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%repcel=1d1
          node(i)%req=req*adjustreq ; node(i)%reqcr=node(i)%req
    else                      !apical
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%repcel=1d1
      node(i)%req=req ; node(i)%reqcr=node(i)%req
        end if
!    print*,"node(i)%tipus",node(i)%tipus   
!    if (node(i)%tipus==1) then
!   
!    print*,"node(i)%req1",node(i)%req
!    else       
!   
!    print*,"node(i)%req2",node(i)%req       
!    end if
   
    node(i)%reqs=req
        node(i)%da=node(i)%req*1.80;
        node(i)%ke=1d1
        node(i)%tor=43d0
        node(i)%stor=53d0
        node(i)%mo=temp
        node(i)%dmo=desmax
        node(i)%diffe=0d0
        node(i)%khold=1d1
    node(i)%kvol=1d0
    node(i)%kplast=1.9d-1
      end do
    end if
    print*,"req",req
    print*,"adjustreq",adjustreq

    !Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=5d0 ; node(i)%adh=5d0
        node(i)%rep=1d1 ; node(i)%repcel=1d1
        node(i)%req=req ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.70
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=1d0  !only for epithelium
        node(i)%stor=4.3d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=desmax
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)
    maxdidare=node(1)%req*0.1d0


    !Distribution of nodes in space
!    call epiteli_sphere(radi,radicel,radius)
!    do i=1,nd
!      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
!    end do
    !do i=1,nd
    !  call random_number(a)
    !  node(i)%x=node(i)%x+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%y=node(i)%y+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%z=node(i)%z+2*a*desmax-desmax
    !end do
   
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%hold=0
    !!let's do it for the most external layer of cells
    j=0
    !print*,"radicel", radicel
      do i=1,(radicel)-2
        j=j+i
      end do
    !print*,"j",j
      j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    print*,"j2",j     
    do i=j+1,ndepi
        node(i)%hold=1 ;node(i)%repcel=1d1!
    !node(i)%border=1 !borde abierto
      end do
!    print*,"ndepi", ndepi
!    print*,"cels(1)%nunodes",cels(1)%nunodes

    !do i=ndepi+1,nd        !all mesenchyme hold
     !node(i)%hold=1 ;node(i)%repcel=1d2
    !end do
 !     j=0
 !     do i=1,radicel-2
 !       j=j+i
 !     end do
 !     print*,"j",j
 !     j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
!     print*,"j2",j
!     print*,"nunodes", cels(1)%nunodes   
!    do i=j+1,ndepi
!        node(i)%hold=1
    !    !node(i)%da=node(i)%req*2.0
    !    !node(i)%orix=0 ; node(i)%oriy=0
    !    !node(i)%oriz=node(i)%oriz-100
    !    !node(i)%rep=1d1;node(i)%repcel=1d1
    !    !node(i)%ke=1d1
    !    !node(i)%tor=1d1
    !    !node(i)%stor=1d1
!      end do
      !j=0
      !do i=1,mradicel-2
      !  j=j+i
      !end do
      !j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      !do i=ndepi+j+1,ndepi+ndepi+ndmes/2
      !  node(i)%hold=1
      !  !node(i)%da=node(i)%req*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%repcel=1d1
      !end do
      !do i=ndepi+ndmes/2+j+1,ndepi+ndmes
      !  node(i)%hold=1
      !  !node(i)%da=node(i)%req*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%repcel=1d1
      !end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
print*,"ncelsepi",ncelsepi
print*,"ncels",ncels 
    if(radi>0.and.(radicel)>0)then
        
    do i=1,ncelsepi
        cels(i)%fase=0d0
      !  cels(i)%reqmax=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
print*,"ncelsepi2",ncelsepi
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
     !   cels(i)%reqmax=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=8
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu
    gen(1)%diffu=4d0 ; gen(1)%kindof=1 ; gen(1)%mu=1d-1     !adhesion epi., house keeping
    gen(2)%diffu=2d0 ; gen(2)%kindof=1 ; gen(2)%mu=1d-1    !adhesion mesen., house keeping
    gen(8)%diffu=1d0 ; gen(8)%kindof=1 ; gen(8)%mu=1d0
    gen(6)%diffu=1d0 ; gen(6)%kindof=1 ; gen(6)%mu=1d0    !celular div. mesen.
   
    gen(4)%diffu=1d0 ; gen(4)%kindof=2 ; gen(4)%mu=1d0    !pre extracelular signal
    gen(4)%npost=1  ; allocate(gen(4)%post(gen(4)%npost)) ; gen(4)%post(1)=5
   
    gen(5)%diffu=5d1 ; gen(5)%kindof=4 ; gen(5)%mu=1d1
    gen(5)%npre=1  ; allocate(gen(5)%pre(gen(5)%npre)) ; gen(5)%pre(1)=4    !extracel signal, promotes division from epi

    gen(3)%diffu=1d1 ; gen(3)%kindof=2 ; gen(3)%mu=1d0    !pre division epi
    gen(3)%npost=1  ; allocate(gen(3)%post(gen(3)%npost)) ; gen(3)%post(1)=7

    gen(7)%diffu=1d1 ; gen(7)%kindof=3 ; gen(7)%mu=5d1    !promotes epi division (activated by 5)
    gen(7)%npre=1  ; allocate(gen(7)%pre(gen(7)%npre)) ; gen(7)%pre(1)=3   

   
    !Gene-behavior interactions
    gen(1)%wa(1)=1   
    gen(2)%wa(1)=2
    gen(7)%wa(nparam_per_node+2)=2.8999d-1
       
    gen(6)%wa(nparam_per_node+2)=5.9d-4!40********************************************2.4d-4
    gen(8)%wa(nparam_per_node+2)=4d-5!40*********************************************5d-5

   

    !Gene-gene interactions
        gen(1)%w(1)=5.0d1
    gen(3)%w(1)=1.0d0
    gen(8)%w(1)=1.0d2
        gen(2)%w(2)=5.0d1
    gen(6)%w(2)=1.0d2
    gen(4)%w(2)=1d1
           
    gen(2)%nww=1
    gen(2)%ww(1,1)=4
    gen(2)%ww(1,2)=5
    gen(2)%ww(1,3)=1d2 

    gen(5)%nww=1
    gen(5)%ww(1,1)=3
    gen(5)%ww(1,2)=7
    gen(5)%ww(1,3)=1d2 

    !Adhesion molecules

    ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions
      kadh(1,1)=1d0!*****************************************************************orig de 40 3d0
      kadh(1,2)=10d0; kadh(2,1)=kadh(1,2)
      kadh(2,2)=3d0
   
    end if

    !Gene expression on nodes
      j=0
    do ii=1,radicel-2
      j=j+ii
    end do
    j=(6*j+1) !;print*,"jota",j
    do i=1,ndepi
    gex(i,1)=5d0
        !gex(i,3)=5d0
    !if (node(i)%tipus==2) gex(i,7)=2d0
    end do

    do i=ndepi+1,nd
        gex(i,2)=5d0
   
    end do

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%da)*3
    node(:)%diffe=0.0d0
   
 

end subroutine tooth40a
!*************************************************************************************************************************************
!*****************************************************************
!*****************************************************************
!**********************************************************************************************************************************************
subroutine tooth40
real*8::pericel,radio2,req1,req2,adjustreq,di
!******* #1 DEFINING SPATIAL DIMENSIONS *******
   
    !epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=6    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer
   ! radius=0    ! no especificar
    !mesenchyme's dimension parameters
    packed=0     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1    !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=2   !number of radial cell layers
    layer=2      !number of planar cell layers
    zmes=1.30  !z-position of uppermost layer
    req=0.15
    !pericel=(pi*radius*radius)/(radicel*2-1)
    !req1=pericel
    !radio2=radius-pericel
   ! req2=(pi*radio2*radio2)/(radicel*2-1)
    radius=(req*(2*radicel-1)/pi) !external
    print*,"radius",radius
    radio2=radius-req! 
    print*,"radio2",radio2
    req2=((pi*radio2)/(2*radicel-1))!*0.9
    print*,"req2",req2
    adjustreq= req2/req
    radius=radius!*1.4
    print*,"radius*1.5",radius
    di=(req+req2)!*1.3

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2    ;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1    !cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
    nd=ndepi+ndmes
    ncels=ncelsepi+ncelsmes
print*,"ncels",ncels


    nda=nd+10
    ncals=ncels+10
print*,"ncals",ncals  
    single=0
    if(radi==1.or.mradi==1) single=1
  !End initializing dimension parameters
nodecel=2
  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals))
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=6
    idum=-8724683
    idumoriginal=idum
    realtime=0

   ! nparam=35       !number of params
    nvarglobal_out=5

    !physical
    temp=1d-2 !low value is low temperature   
    desmax=0.001
    resmax=1d-3
    nuns=100  !for screening
    !angleto=90 ; angltormax=cos(angleto*2*pi/360) !maximum angle for epithelial torsion  
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13  
    dmax=2
    khold=0d0
    ecmmax=req
    !biological
    !functions used
   
   
    ffu=0
    ffu(1)=1 !spring of the ellipse
    ffu(2)=0 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !eggshell miguel4-1-13

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(17)=1 !conservacion del volumen
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise
    ffu(16)=1
   

!**************Distribution of nodes in space
  
    if(radi>0.and.radicel>0) call epiteli_sphere1(radi,radicel,radius,di,zepi)
    if(mradi>0.and.mradicel>0.and.layer>0) then
                       if(packed==0)then ; call mesenq1(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if
    print*,"node(ndepi+1)%req",node(ndepi+1)%req
    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do

    !do i=1,nd
    !  call random_number(a)
    !  node(i)%x=node(i)%x+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%y=node(i)%y+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%z=node(i)%z+2*a*desmax-desmax
    !end do
 
    !End distribution of nodes in space

   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
           if(node(i)%tipus==2)then  !basal                    !         if(i<=ndepi/2)then  !basal
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%repcel=1d1
          node(i)%req=req*adjustreq ; node(i)%reqcr=node(i)%req
    else                      !apical
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%repcel=1d1
      node(i)%req=req ; node(i)%reqcr=node(i)%req
        end if
!    print*,"node(i)%tipus",node(i)%tipus   
!    if (node(i)%tipus==1) then
!   
!    print*,"node(i)%req1",node(i)%req
!    else       
!   
!    print*,"node(i)%req2",node(i)%req       
!    end if
   
    node(i)%reqs=req
        node(i)%da=node(i)%req*1.80;
        node(i)%ke=1d1
        node(i)%tor=43d0
        node(i)%stor=43d0
        node(i)%mo=temp
        node(i)%dmo=desmax
        node(i)%diffe=0d0
        node(i)%khold=1d1
    node(i)%kvol=1d0
    node(i)%kplast=1.9d-1
      end do
    end if
    print*,"req",req
    print*,"adjustreq",adjustreq

    !Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=5d0 ; node(i)%adh=5d0
        node(i)%rep=1d1 ; node(i)%repcel=1d1
        node(i)%req=req ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.70
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=1d0  !only for epithelium
        node(i)%stor=4.3d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=desmax
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)
    maxdidare=node(1)%req*0.1d0


    !Distribution of nodes in space
!    call epiteli_sphere(radi,radicel,radius)
!    do i=1,nd
!      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
!    end do
    !do i=1,nd
    !  call random_number(a)
    !  node(i)%x=node(i)%x+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%y=node(i)%y+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%z=node(i)%z+2*a*desmax-desmax
    !end do
   
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%hold=0
    !!let's do it for the most external layer of cells
    j=0
    !print*,"radicel", radicel
      do i=1,(radicel)-2
        j=j+i
      end do
    !print*,"j",j
      j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    print*,"j2",j     
    do i=j+1,ndepi
        node(i)%hold=1 ;node(i)%repcel=1d1!
    !node(i)%border=1 !borde abierto
      end do
!    print*,"ndepi", ndepi
!    print*,"cels(1)%nunodes",cels(1)%nunodes

    !do i=ndepi+1,nd        !all mesenchyme hold
     !node(i)%hold=1 ;node(i)%repcel=1d2
    !end do
 !     j=0
 !     do i=1,radicel-2
 !       j=j+i
 !     end do
 !     print*,"j",j
 !     j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
!     print*,"j2",j
!     print*,"nunodes", cels(1)%nunodes   
!    do i=j+1,ndepi
!        node(i)%hold=1
    !    !node(i)%da=node(i)%req*2.0
    !    !node(i)%orix=0 ; node(i)%oriy=0
    !    !node(i)%oriz=node(i)%oriz-100
    !    !node(i)%rep=1d1;node(i)%repcel=1d1
    !    !node(i)%ke=1d1
    !    !node(i)%tor=1d1
    !    !node(i)%stor=1d1
!      end do
      !j=0
      !do i=1,mradicel-2
      !  j=j+i
      !end do
      !j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      !do i=ndepi+j+1,ndepi+ndepi+ndmes/2
      !  node(i)%hold=1
      !  !node(i)%da=node(i)%req*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%repcel=1d1
      !end do
      !do i=ndepi+ndmes/2+j+1,ndepi+ndmes
      !  node(i)%hold=1
      !  !node(i)%da=node(i)%req*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%repcel=1d1
      !end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
print*,"ncelsepi",ncelsepi
print*,"ncels",ncels 
    if(radi>0.and.(radicel)>0)then
        
    do i=1,ncelsepi
        cels(i)%fase=0d0
      !  cels(i)%reqmax=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
print*,"ncelsepi2",ncelsepi
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
     !   cels(i)%reqmax=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=9
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu
    gen(1)%diffu=4d0 ; gen(1)%kindof=1 ; gen(1)%mu=1d-1     !adhesion epi., house keeping
    gen(2)%diffu=2d0 ; gen(2)%kindof=1 ; gen(2)%mu=1d-1    !adhesion mesen., house keeping
    gen(8)%diffu=1d0 ; gen(8)%kindof=1 ; gen(8)%mu=1d0
    gen(6)%diffu=1d0 ; gen(6)%kindof=1 ; gen(6)%mu=1d0    !celular div. mesen.
   
    gen(4)%diffu=1d0 ; gen(4)%kindof=2 ; gen(4)%mu=1d0    !pre extracelular signal
    gen(4)%npost=1  ; allocate(gen(4)%post(gen(4)%npost)) ; gen(4)%post(1)=5
   
    gen(5)%diffu=5d1 ; gen(5)%kindof=4 ; gen(5)%mu=1d1
    gen(5)%npre=1  ; allocate(gen(5)%pre(gen(5)%npre)) ; gen(5)%pre(1)=4    !extracel signal, promotes division from epi

    gen(3)%diffu=1d1 ; gen(3)%kindof=2 ; gen(3)%mu=1d0    !pre division epi
    gen(3)%npost=1  ; allocate(gen(3)%post(gen(3)%npost)) ; gen(3)%post(1)=7

    gen(7)%diffu=1d1 ; gen(7)%kindof=3 ; gen(7)%mu=5d1    !promotes epi division (activated by 5)
    gen(7)%npre=1  ; allocate(gen(7)%pre(gen(7)%npre)) ; gen(7)%pre(1)=3   

   

    gen(9)%diffu=5d-2 ; gen(9)%kindof=4 ; gen(9)%mu=5d-1

    !Gene-behavior interactions
    gen(1)%wa(1)=1   
    gen(2)%wa(1)=2
    gen(7)%wa(nparam_per_node+2)=2.8999d-1
       
    gen(6)%wa(nparam_per_node+2)=2.4d-4!40********************************************2.4d-4
    gen(8)%wa(nparam_per_node+2)=5d-5!40*********************************************5d-5

   

    !Gene-gene interactions
        gen(1)%w(1)=5.0d1
    gen(3)%w(1)=5.0d1
    gen(8)%w(1)=1.0d2
        gen(2)%w(2)=5.0d1
    gen(6)%w(2)=1.0d2
    gen(4)%w(2)=1d1
           
    gen(2)%nww=1
    gen(2)%ww(1,1)=4
    gen(2)%ww(1,2)=5
    gen(2)%ww(1,3)=1d2 

    gen(5)%nww=1
    gen(5)%ww(1,1)=3
    gen(5)%ww(1,2)=7
    gen(5)%ww(1,3)=1d2 

    gen(9)%nww=1
    gen(9)%ww(1,1)=3
    gen(9)%ww(1,2)=7
    gen(9)%ww(1,3)=1d2 

    !Adhesion molecules

    ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions
      kadh(1,1)=3d0!*****************************************************************orig de 40 3d0
      kadh(1,2)=10d0; kadh(2,1)=kadh(1,2)
      kadh(2,2)=3d0
   
    end if

    !Gene expression on nodes
      j=0
    do ii=1,radicel-2
      j=j+ii
    end do
    j=(6*j+1) !;print*,"jota",j
    do i=1,ndepi
    gex(i,1)=5d0
        !gex(i,3)=5d0
    !if (node(i)%tipus==2) gex(i,7)=2d0
    end do

    do i=ndepi+1,nd
        gex(i,2)=5d0
   
    end do

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%da)*3
    node(:)%diffe=0.0d0
   
 

end subroutine tooth40

subroutine tooth40_ligand
real*8::pericel,radio2,req1,req2,adjustreq,di
!******* #1 DEFINING SPATIAL DIMENSIONS *******
   
    !epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=6    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer
   ! radius=0    ! no especificar
    !mesenchyme's dimension parameters
    packed=0     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1    !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=2   !number of radial cell layers
    layer=2      !number of planar cell layers
    zmes=1.30  !z-position of uppermost layer
    req=0.15
    !pericel=(pi*radius*radius)/(radicel*2-1)
    !req1=pericel
    !radio2=radius-pericel
   ! req2=(pi*radio2*radio2)/(radicel*2-1)
    radius=(req*(2*radicel-1)/pi) !external
    print*,"radius",radius
    radio2=radius-req! 
    print*,"radio2",radio2
    req2=((pi*radio2)/(2*radicel-1))!*0.9
    print*,"req2",req2
    adjustreq= req2/req
    radius=radius!*1.4
    print*,"radius*1.5",radius
    di=(req+req2)!*1.3

!************************************************


    !Initializing dimension parameters
    if(radi>0.and.radicel>0)then
      j=0
      do i=1,radi-1
        j=j+i
      end do
      nodecel=(6*j+1)*2    ;if(radi==1) nodecel=2!number of nodes per cell
      nodecela=2*nodecel+1    !cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13

      j=0
      do i=1,radicel-1
        j=j+i
      end do
      ncelsepi=(6*j+1)
      ndepi=nodecel*ncelsepi
    else
      ncelsepi=0
      ndepi=0
      nodecel=0
      nodecela=0
    end if

    if(mradi>0.and.mradicel>0.and.layer>0)then
      if(packed==1)then !if mesenchyme is packed, we use mradi differently
        j=0
        do i=1,mradi-1
          j=j+i
        end do
        nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
        nodecela=2*nodecel+1
      else
        nodecel=mradi
        nodecela=2*nodecel+1
      end if

      j=0
      do i=1,mradicel-1
        j=j+i
      end do
      ncelsmes=(6*j+1)*layer  !number of mesenchymal cells
      ndmes=nodecel*ncelsmes    !number of mesenchymal nodes
      if(radi==0.and.radicel==0)then
        nodecel=mradi
        nodecela=2*nodecel+1
      else if(nodecel<mradi)then
        nodecel=radi
        nodecela=nodecel*2+1
      end if
    else
      ndmes=0
      ncelsmes=0
    end if
    nd=ndepi+ndmes
    ncels=ncelsepi+ncelsmes
print*,"ncels",ncels


    nda=nd+10
    ncals=ncels+10
print*,"ncals",ncals  
    single=0
    if(radi==1.or.mradi==1) single=1
  !End initializing dimension parameters
nodecel=2
  print*,"nd",nd,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals))
    call iniarrays
    !End Allocatations


   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=6
    idum=-8724683
    idumoriginal=idum
    realtime=0

   ! nparam=35       !number of params
    nvarglobal_out=5

    !physical
    temp=1d-2 !low value is low temperature   
    desmax=0.001
    resmax=1d-3
    nuns=100  !for screening
    !angleto=90 ; angltormax=cos(angleto*2*pi/360) !maximum angle for epithelial torsion  
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13  
    dmax=2
    khold=0d0
    ecmmax=req
    mnn=1000
    !biological
    !functions used
   
   
    ffu=0
    ffu(1)=1 !spring of the ellipse
    ffu(2)=0 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !eggshell miguel4-1-13

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(17)=1 !conservacion del volumen
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise
    ffu(16)=1
   

!**************Distribution of nodes in space
  
    if(radi>0.and.radicel>0) call epiteli_sphere1(radi,radicel,radius,di,zepi)
    if(mradi>0.and.mradicel>0.and.layer>0) then
                       if(packed==0)then ; call mesenq1(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if
    print*,"node(ndepi+1)%req",node(ndepi+1)%req
    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do

    !do i=1,nd
    !  call random_number(a)
    !  node(i)%x=node(i)%x+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%y=node(i)%y+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%z=node(i)%z+2*a*desmax-desmax
    !end do
 
    !End distribution of nodes in space

   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    if(radi>0.and.radicel>0)then
      do i=1,ndepi
           if(node(i)%tipus==2)then  !basal                    !         if(i<=ndepi/2)then  !basal
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%repcel=1d1
          node(i)%req=req*adjustreq ; node(i)%reqcr=node(i)%req
    else                      !apical
          node(i)%you=1d1 ; node(i)%adh=0d0
          node(i)%rep=1d1 ; node(i)%repcel=1d1
      node(i)%req=req ; node(i)%reqcr=node(i)%req
        end if
!    print*,"node(i)%tipus",node(i)%tipus   
!    if (node(i)%tipus==1) then
!   
!    print*,"node(i)%req1",node(i)%req
!    else       
!   
!    print*,"node(i)%req2",node(i)%req       
!    end if
   
    node(i)%reqs=req
        node(i)%da=node(i)%req*1.80;
        node(i)%ke=1d1
        node(i)%tor=43d0
        node(i)%stor=43d0
        node(i)%mo=temp
        node(i)%dmo=desmax
        node(i)%diffe=0d0
        node(i)%khold=1d1
    node(i)%kvol=1d0
    node(i)%kplast=1.9d-1
      end do
    end if
    print*,"req",req
    print*,"adjustreq",adjustreq

    !Mesenchyme
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ndepi+1,nd
        node(i)%you=5d0 ; node(i)%adh=5d0
        node(i)%rep=1d1 ; node(i)%repcel=1d1
        node(i)%req=req ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.70
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=1d0  !only for epithelium
        node(i)%stor=4.3d1 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=desmax
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)
    maxdidare=node(1)%req*0.1d0


    !Distribution of nodes in space
!    call epiteli_sphere(radi,radicel,radius)
!    do i=1,nd
!      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
!    end do
    !do i=1,nd
    !  call random_number(a)
    !  node(i)%x=node(i)%x+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%y=node(i)%y+2*a*desmax-desmax
    !  call random_number(a)
    !  node(i)%z=node(i)%z+2*a*desmax-desmax
    !end do
   
    !End distribution of nodes in space

    !Setting boundary nodes (they will be subject to an external force keeping them in their initial position)
    node(:)%hold=0
    !!let's do it for the most external layer of cells
    j=0
    !print*,"radicel", radicel
      do i=1,(radicel)-2
        j=j+i
      end do
    !print*,"j",j
      j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    print*,"j2",j     
    do i=j+1,ndepi
        node(i)%hold=1 ;node(i)%repcel=1d1!
    !node(i)%border=1 !borde abierto
      end do
!    print*,"ndepi", ndepi
!    print*,"cels(1)%nunodes",cels(1)%nunodes

    !do i=ndepi+1,nd        !all mesenchyme hold
     !node(i)%hold=1 ;node(i)%repcel=1d2
    !end do
 !     j=0
 !     do i=1,radicel-2
 !       j=j+i
 !     end do
 !     print*,"j",j
 !     j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
!     print*,"j2",j
!     print*,"nunodes", cels(1)%nunodes   
!    do i=j+1,ndepi
!        node(i)%hold=1
    !    !node(i)%da=node(i)%req*2.0
    !    !node(i)%orix=0 ; node(i)%oriy=0
    !    !node(i)%oriz=node(i)%oriz-100
    !    !node(i)%rep=1d1;node(i)%repcel=1d1
    !    !node(i)%ke=1d1
    !    !node(i)%tor=1d1
    !    !node(i)%stor=1d1
!      end do
      !j=0
      !do i=1,mradicel-2
      !  j=j+i
      !end do
      !j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
      !do i=ndepi+j+1,ndepi+ndepi+ndmes/2
      !  node(i)%hold=1
      !  !node(i)%da=node(i)%req*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%repcel=1d1
      !end do
      !do i=ndepi+ndmes/2+j+1,ndepi+ndmes
      !  node(i)%hold=1
      !  !node(i)%da=node(i)%req*2.0
      !  !node(i)%orix=0 ; node(i)%oriy=0
      !  !node(i)%rep=1d1;node(i)%repcel=1d1
      !end do

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
print*,"ncelsepi",ncelsepi
print*,"ncels",ncels 
    if(radi>0.and.(radicel)>0)then
        
    do i=1,ncelsepi
        cels(i)%fase=0d0
      !  cels(i)%reqmax=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if
    !Mesenchymal
print*,"ncelsepi2",ncelsepi
    if(mradi>0.and.mradicel>0.and.layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
     !   cels(i)%reqmax=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do

   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=10
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu
    gen(1)%diffu=4d0 ; gen(1)%kindof=1 ; gen(1)%mu=1d-1     !adhesion epi., house keeping
    gen(2)%diffu=2d0 ; gen(2)%kindof=1 ; gen(2)%mu=1d-1    !adhesion mesen., house keeping
    gen(8)%diffu=1d0 ; gen(8)%kindof=1 ; gen(8)%mu=1d0
    gen(6)%diffu=1d0 ; gen(6)%kindof=1 ; gen(6)%mu=1d0    !celular div. mesen.
   
    gen(4)%diffu=1d0 ; gen(4)%kindof=2 ; gen(4)%mu=1d0    !pre extracelular signal
    gen(4)%npost=1  ; allocate(gen(4)%post(gen(4)%npost)) ; gen(4)%post(1)=5
   
    gen(5)%diffu=5d1 ; gen(5)%kindof=4 ; gen(5)%mu=1d1
    gen(5)%npre=1  ; allocate(gen(5)%pre(gen(5)%npre)) ; gen(5)%pre(1)=4    !extracel signal, promotes division from epi

    gen(3)%diffu=1d1 ; gen(3)%kindof=2 ; gen(3)%mu=1d0    !pre division epi
    !gen(3)%npost=1  ; allocate(gen(3)%post(gen(3)%npost)) ; gen(3)%post(1)=7
    !gen(3)%npre=1  ; allocate(gen(3)%pre(gen(3)%npre)) ; gen(3)%pre(1)=7


    gen(7)%diffu=0d0 ; gen(7)%kindof=8 ; gen(7)%mu=0d0    !receptor activated by 5 promotes epi division (activated by 5)
    gen(7)%npre=2  ; allocate(gen(7)%pre(gen(7)%npre)) ; gen(7)%pre(1)=3 ; gen(7)%pre(2)=5
    gen(7)%npost=2  ; allocate(gen(7)%post(gen(7)%npost)) ; gen(7)%post(1)=3 ; gen(7)%post(2)=5

   

    gen(9)%diffu=5d-2 ; gen(9)%kindof=4 ; gen(9)%mu=5d-1  !signal from bead

    gen(10)%diffu=0d0 ; gen(10)%kindof=8 ; gen(10)%mu=0d0    !receptor activated by 5 promotes epi division (activated by 5)
    gen(10)%npre=2  ; allocate(gen(10)%pre(gen(10)%npre)) ; gen(10)%pre(1)=3 ; gen(10)%pre(2)=9
    gen(10)%npost=2  ; allocate(gen(10)%post(gen(10)%npost)) ; gen(10)%post(1)=3 ; gen(10)%post(2)=9

    !Gene-behavior interactions
    gen(1)%wa(1)=1   
    gen(2)%wa(1)=2
    gen(7)%wa(nparam_per_node+2)=2.8999d-1
       
    gen(6)%wa(nparam_per_node+2)=2.4d-4!40********************************************2.4d-4
    gen(8)%wa(nparam_per_node+2)=5d-5!40*********************************************5d-5

   

    !Gene-gene interactions
        gen(1)%w(1)=5.0d1
    gen(3)%w(1)=5.0d1
    gen(8)%w(1)=1.0d2
        gen(2)%w(2)=5.0d1
    gen(6)%w(2)=1.0d2
    gen(4)%w(2)=1d1
           
    gen(2)%nww=1
    gen(2)%ww(1,1)=4
    gen(2)%ww(1,2)=5
    gen(2)%ww(1,3)=1d2 

    !gen(5)%nww=1
    !gen(5)%ww(1,1)=3
    !gen(5)%ww(1,2)=7
    !gen(5)%ww(1,3)=1d2 

    !gen(9)%nww=1
    !gen(9)%ww(1,1)=3
    !gen(9)%ww(1,2)=7
    !gen(9)%ww(1,3)=1d2 


    gen(7)%nww=4
    gen(7)%ww(1,1)=3
    gen(7)%ww(1,2)=7
    gen(7)%ww(1,3)=1d2 
    gen(7)%ww(2,1)=7
    gen(7)%ww(2,2)=3
    gen(7)%ww(2,3)=1d2
    gen(7)%ww(3,1)=5
    gen(7)%ww(3,2)=7
    gen(7)%ww(3,3)=1d2 
    gen(7)%ww(4,1)=7
    gen(7)%ww(4,2)=5
    gen(7)%ww(4,3)=1d2 


    gen(10)%nww=4
    gen(10)%ww(1,1)=3
    gen(10)%ww(1,2)=10
    gen(10)%ww(1,3)=1d2 
    gen(10)%ww(2,1)=10
    gen(10)%ww(2,2)=3
    gen(10)%ww(2,3)=1d2
    gen(10)%ww(3,1)=9
    gen(10)%ww(3,2)=10
    gen(10)%ww(3,3)=1d2 
    gen(10)%ww(4,1)=10
    gen(10)%ww(4,2)=9
    gen(10)%ww(4,3)=1d2 


    !Adhesion molecules

    ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions
      kadh(1,1)=3d0!*****************************************************************orig de 40 3d0
      kadh(1,2)=10d0; kadh(2,1)=kadh(1,2)
      kadh(2,2)=3d0
   
    end if

    !Gene expression on nodes
      j=0
    do ii=1,radicel-2
      j=j+ii
    end do
    j=(6*j+1) !;print*,"jota",j
    do i=1,ndepi
    gex(i,1)=5d0
    gex(i,3)=5d-1

        !gex(i,3)=5d0
    !if (node(i)%tipus==2) gex(i,7)=2d0
    end do

    do i=ndepi+1,nd
        gex(i,2)=5d0
   
    end do

    call update_npag
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%da)*3
    node(:)%diffe=0.0d0
   
 

end subroutine tooth40_ligand




!*************************************************************************************************************************************
!*****************************************************************


subroutine epiteli_sphere1(radi,radicel,radius,di,zepi)        !radi is the number of cell rows from one pole to the equator
integer            ::radi,radicel,valcel

integer::trobat,cont,rep,n,i,j,k,ii,jj,kk,lev,numdos,ncels2,val,val2,ll,val3,val4,val5,val6,iiii,jjjj,kkkk,suma,cocels
real*8::modul,pescab,pescac,pescbc,modula,modulb,modulc,aax,aay,aaz,bbx,bby,bbz,ccx,ccy,ccz,costat,rhex,lad,ax,ay,az,di
real*8::alf,bet,gam,a,b,c,l,aa,bb,cc,aaa,bbb,ccc,angle,rpent,modmin,pesc,radius,radius2,bx,by,bz,cost,ucost,sint,thet,radicelepisphere
real*8,dimension(:)::minu(3),vec1(3),vec2(3),vec3(3),vec4(3),vec5(3)
real*8,allocatable::cmalla(:,:),cveci(:,:),primers(:)
real*8::zepi




    di=di*0.80
    de=0.5d0

    radicelepisphere=(radicel-1)*2

    !SPHERE CODE, FIRST HEMISPHERE
    cont=0 ; cocels=0
    ii=radicelepisphere/2
!    radius=5d0  !radius of the sphere
    radius2=radius+di

    alf=pi*0.5d0/radicelepisphere

!!!pole cell
    cont=cont+1 ; cocels=cocels+1
    node(1)%x=0d0 ; node(1)%y=0d0 ; node(1)%z=radius
    node(1)%altre=2 ; node(1)%tipus=2 ; node(1)%icel=cocels
    cont=cont+1
    node(2)%x=0d0 ; node(2)%y=0d0 ; node(2)%z=radius2
    node(2)%altre=1 ; node(2)%tipus=1 ; node(2)%icel=cocels

    !el radi 2 de la cèlula!
    if (radi>1) then !111
    cont=cont+1    !pone un nodo adyacente al inicial
    node(cont)%x=0d0 ; node(cont)%y=de ; node(cont)%z=radius
    node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
    cont=cont+1
    node(cont)%x=0d0 ; node(cont)%y=de ; node(cont)%z=radius2
    node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels

    gam=2*pi/6d0
    do i=1,5  !pone resto de nodos
      cont=cont+1
      node(cont)%x=de*cos(i*gam+pi/2) ; node(cont)%y=de*sin(i*gam+pi/2) ; node(cont)%z=radius
      node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
      cont=cont+1
      node(cont)%x=de*cos(i*gam+pi/2) ; node(cont)%y=de*sin(i*gam+pi/2) ; node(cont)%z=radius2
      node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels
    end do
    end if !111
!!!!! end pole cell

!!!!! spine cell
    do i=1,ii
       
!  print*,"*********espina******",i
        cont=cont+1 ; cocels=cocels+1    !hace columna de celulas, 1. nodo
        suma=0
        jj=cont
        angle=alf*i
!        print*,"cont",cont
        vec1(1)=radius*sin(angle)
        vec1(2)=0
        vec1(3)=radius*cos(angle)
        node(cont)%x=vec1(1) ; node(cont)%y=vec1(2) ; node(cont)%z=vec1(3)
        node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
!        print*,"altre",node(cont)%altre

        cont=cont+1
        vec2(1)=radius2*sin(angle)
        vec2(2)=0
        vec2(3)=radius2*cos(angle)
        node(cont)%x=vec2(1) ; node(cont)%y=vec2(2) ; node(cont)%z=vec2(3)
        node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels

        ux=vec2(1)-vec1(1) ; uy=vec2(2)-vec1(2) ; uz=vec2(3)-vec1(3)
        a=1d0/sqrt(ux**2+uy**2+uz**2) ; ux=ux*a ; uy=uy*a ; uz=uz*a

        ax=-de*sin(0d0) ; ay=de*cos(0d0) ; az=0d0
!        print*,"vora"," a",ax,ay,az

    if (radi>1) then !111
        do j=1,6    !pone nodos restantes alrededor 1. nodo
          cont=cont+1
          thet=j*gam
   print*,"thet",thet
          cost=cos(thet); ucost=1-cost ; sint=sin(thet)
  print*,"cos",cost,"ucost",ucost,"sint",sint
          bx=(cost+ux**2*ucost)*ax+(ux*uy*ucost-uz*sint)*ay+(ux*uz*ucost+uy*sint)*az
          by=(uy*ux*ucost+uz*sint)*ax+(cost+uy**2*ucost)*ay+(uy*uz*ucost-ux*sint)*az
          bz=(ux*uz*ucost-uy*sint)*ax+(uy*uz*ucost+ux*sint)*ay+(cost+uz**2*ucost)*az
          node(cont)%x=vec1(1)+bx ; node(cont)%y=vec1(2)+by ; node(cont)%z=vec1(3)+bz
          node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
!   print*,"b",bx,by,bz
!  print*,"vec1",vec1
!  print*,"nodecont",node(cont)%x,node(cont)%y,node(cont)%z
          cont=cont+1

          node(cont)%x=vec2(1)+bx ; node(cont)%y=vec2(2)+by ; node(cont)%z=vec2(3)+bz
          node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels
        end do
    end if !111
        !!!!!!!!! end spine cell



!        print*,"altre",node(cont)%altre
      if (radi>1) then  !111
        if(i+1<=ii+1)then
            kk=6+6*(i-1)
        else
            kk=6+6*(radicelepisphere-i-1)
    end if !111
      else    !111
    kk=6+6*(i-1)    !111
      end if
   
        bet=2*pi/kk

!        print*,"i",i,"nombre de cels del paralel",kk

        do j=1,kk-1

            !!!!!!!!!!! rib cell

            !if(i>15 .and. i<20)then
            !  !if(j==1.or.mod(j,kk/6)==0.or.mod(j+1,kk/6)==0.or.mod(j-1,kk/6)==0)then
            !  if(j==1 .or. mod(j+1,kk/6)==0 .or. mod(j-1,kk/6)==0)then
            !    nd=nd-2 ; ndepi=ndepi-2
            !    ncels=ncels-1 ; ncelsepi=ncelsepi-1
            !    cycle
            !  end if
            !end if
            angle=alf

            cont=cont+1 ; cocels=cocels+1
            a=node(1)%z*sin(i*angle)!+node(1)%x*cos(i*angle)!fem els 2 girs respecte al node 1 (pol nord)
            b=node(1)%y
            c=node(1)%z*cos(i*angle)!-node(1)%x*sin(i*angle)
            node(cont)%x=a
            node(cont)%y=b
            node(cont)%z=c
            a=node(cont)%x*cos(j*bet)-node(cont)%y*sin(j*bet)
            b=node(cont)%x*sin(j*bet)+node(cont)%y*cos(j*bet)
            c=node(cont)%z
            node(cont)%x=a
            node(cont)%y=b
            node(cont)%z=c
            node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels
            vec1(1)=node(cont)%x ; vec1(2)=node(cont)%y ; vec1(3)=node(cont)%z

            cont=cont+1
            a=node(2)%z*sin(i*angle)!+node(1)%x*cos(i*angle)!fem els 2 girs respecte al node 1 (pol nord)
            b=node(2)%y
            c=node(2)%z*cos(i*angle)!-node(1)%x*sin(i*angle)
            node(cont)%x=a
            node(cont)%y=b
            node(cont)%z=c
            a=node(cont)%x*cos(j*bet)-node(cont)%y*sin(j*bet)
            b=node(cont)%x*sin(j*bet)+node(cont)%y*cos(j*bet)
            c=node(cont)%z
            node(cont)%x=a
            node(cont)%y=b
            node(cont)%z=c
            node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels

            ux=node(cont)%x-node(cont-1)%x
            uy=node(cont)%y-node(cont-1)%y
            uz=node(cont)%z-node(cont-1)%z
            a=1d0/sqrt(ux**2+uy**2+uz**2) ; ux=ux*a ; uy=uy*a ; uz=uz*a !;print*,"spring",a
            vec2(1)=node(cont)%x ; vec2(2)=node(cont)%y ; vec2(3)=node(cont)%z

            ax=-de*sin(j*bet) ; ay=de*cos(j*bet) ; az=0d0

!        ax=de*sin(0d0) ; ay=de*cos(0d0) ; az=0d0
!        print*,"bet",j*bet,"a",ax,ay,"modul",sqrt(ax**2+ay**2)
    if (radi>1) then           
    do k=1,6
              cont=cont+1
              thet=k*gam
              cost=cos(thet); ucost=1-cost ; sint=sin(thet)
              bx=(cost+ux**2*ucost)*ax+(ux*uy*ucost-uz*sint)*ay+(ux*uz*ucost+uy*sint)*az
              by=(uy*ux*ucost+uz*sint)*ax+(cost+uy**2*ucost)*ay+(uy*uz*ucost-ux*sint)*az
              bz=(ux*uz*ucost-uy*sint)*ax+(uy*uz*ucost+ux*sint)*ay+(cost+uz**2*ucost)*az
              node(cont)%x=vec1(1)+bx ; node(cont)%y=vec1(2)+by ; node(cont)%z=vec1(3)+bz
              node(cont)%altre=cont+1 ; node(cont)%tipus=2 ; node(cont)%icel=cocels

!              print*,"b rib",bx,by,bz,"modul",sqrt(bx**2+by**2+bz**2)

              cont=cont+1

              node(cont)%x=vec2(1)+bx ; node(cont)%y=vec2(2)+by ; node(cont)%z=vec2(3)+bz
              node(cont)%altre=cont-1 ; node(cont)%tipus=1 ; node(cont)%icel=cocels
           end do
    end if


            !!!!!!!!!!! end rib cell



        end do
    end do

    
    do i=1,ndepi
      node(i)%z=node(i)%z+zepi
    end do
    
    
!print*,"cont",cont,"cocels",cocels


!the 2ond hemisphere

!    j=0
!    do i=1,(radicel/2+1)-2
!      j=j+i
!    end do
!    ii=(6*j+1)*2*7   !number of nodes on the 2ond hemisphere, si se quita sale la mitad de la esfera

!    print*,"ii",ii
!    do i=1,ii
!      cont=cont+1
!      node(cont)%x=node(i)%x
!      node(cont)%y=node(i)%y
!      node(cont)%z=-node(i)%z
!      node(cont)%tipus=node(i)%tipus
!      if(node(i)%tipus==2)then
!        node(cont)%altre=cont+1
!      else
!        node(cont)%altre=cont-1
!      end if
!      node(cont)%icel=node(i)%icel+cocels
!    end do

!print*,"cont complet",cont,"cocels",cocels


!    cont=cont+1
!    malla(cont,1)=0;malla(cont,2)=0;malla(cont,3)=-radius
!    primers(radi+1)=cont
!    print*,i,"malla",malla(cont,:)
!    print*,"nombre de cels creades",cont
!    print*,"nombre de paral·lels",lev





    !define the cell's centroid

    do i=1,ncelsepi
      kk=0
      cels(i)%ctipus=1    !tipo celular 1, epitelio
      cels(i)%nunodes=nodecel !
      cels(i)%nodela=nodecela
      cels(i)%cex=0d0 ; cels(i)%cey=0d0 ; cels(i)%cez=0d0
      cels(i)%polx=0d0 ; cels(i)%poly=0d0 ; cels(i)%polz=0d0
      allocate(cels(i)%node(nodecela))
      do j=1,ndepi
        k=node(j)%icel
        if(k==i)then
          kk=kk+1
          cels(i)%node(kk)=j
          if(node(j)%tipus==1)then
            cels(i)%cex=cels(i)%cex+node(j)%x
            cels(i)%cey=cels(i)%cey+node(j)%y
            cels(i)%cez=cels(i)%cez+node(j)%z
          end if
        end if
      end do
      cels(i)%cex=2*cels(i)%cex/cels(i)%nunodes
      cels(i)%cey=2*cels(i)%cey/cels(i)%nunodes
      cels(i)%cez=2*cels(i)%cez/cels(i)%nunodes
    end do

    do i=1,ndepi
      if(node(i)%tipus==2) node(i)%marge=0
      if(node(i)%tipus==1) node(i)%marge=1
    end do


!    do i=1,ncelsepi
!        cels(i)%ctipus=1
!        cels(i)%cex=0;cels(i)%cey=0;cels(i)%cez=0
!        do j=1,cels(i)%nunodes
!            k=cels(i)%node(j)
!            if(node(k)%tipus==1)then
!                 cels(i)%cex=cels(i)%cex+node(k)%x
!                 cels(i)%cey=cels(i)%cey+node(k)%y
!                cels(i)%cez=cels(i)%cez+node(k)%z
!            end if
!        end do
!        cels(i)%cex=2*cels(i)%cex/real(cels(i)%nunodes)
!        cels(i)%cey=2*cels(i)%cey/real(cels(i)%nunodes)
!        cels(i)%cez=2*cels(i)%cez/real(cels(i)%nunodes)
!    end do


 


    realtime=0
    maxdidare=node(1)%req*0.1d0
    dmax=1
    screen_radius=1.0d0
    node(:)%talone=0.0d0
    ramax=maxval(node(:)%da)*3
    node(:)%diffe=0.0d0
end subroutine epiteli_sphere1
!*************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mesenq1(radi,radicel,layer,zmes,dreq)                                    !!! miguel 20.11
integer  ::i,j,k,ii,jj,radi,radicel,layer,signo,ent,ont       ! number of layers and concentric hexagons
real*8   ::rad,der,zmes,de,di,dreq
real*8   :: xx,yy,zz                                                        ! miguel 4-6-13
rad=pi/180d0
dreq=0.05
de=dreq*2                    ! call radius
!di=2.0d0*de                 ! distance between cells and layers (it has to be >2 to avoid cell contacts)
di=2*de+2*de*cos(60d0*rad)
print*,"dreq",dreq
    do i=ncelsepi+1,ncels              
!        cels(i)%nunodes=radi+1
        cels(i)%nunodes=radi                !>>>>>>>>>>>>>>Miquel 21-3-13
        cels(i)%nodela=nodecela
!        allocate(cels(i)%node(radi+1))
        allocate(cels(i)%node(nodecela))        !>>>>>>>>>>>>>>Miquel 21-3-13
        cels(i)%node=0
    end do

node(ndepi+1:)%marge=1

!    print*,"node",size(node),"cels",size(cels),"celsnode",size(cels(1)%node)


    kkk=ndmes/layer
! radi=radi+1
    cont=ncelsepi+1

    ii=ndepi        !node counter
    jj=ncelsepi        !cel counter


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! central nodes (=nucleus) !!!!

   do l=1,layer 
!    jj=(kkk*(l-1))+radi+ndepi+1
    ii=ii+1
    jj=jj+1
!    node(ii)%x=(sqrt(di/2d0))*mod(l+1,2) ; node(ii)%y=(sqrt(di/2d0))*mod(l+1,2) ! Euler packing algorhytm
    node(ii)%x=0d0 ; node(ii)%y=0d0 ! Euler packing algorhytm
    node(ii)%z=zmes-di*(l-1)    !zmes marks the z-position of the uppermost layer
!    print*,"ii",ii,node(ii)%x,node(ii)%y,node(ii)%z
    node(ii)%tipus=3 ; node(ii)%icel=jj                        ! origin   
    cels(jj)%node(1)=ii
!    print*,"cel central, node central. ii:",ii
    ent=ii
    ont=ii
    !fill with the external nodes
    if(radi>1)then
        signo=1                                                              ! miguel 4-6-13
        do k=2,radi
            ii=ii+1
            !print*,"cel perif. ent:",ent,ii,k-1,
             call random_number(a)
            der=(sqrt(de)*sqrt(de*a))                               ! miguel 4-6-13
121         call random_number(a)
            xx=der*(1d0-2*a) ;
            call random_number(a)
            yy=der*(1d0-2*a)            ! miguel 4-6-13
            zz=(der**2)-(xx**2)-(yy**2)                                      ! miguel 4-6-13
            if(zz.lt.0)then;goto 121;endif                                   ! miguel 4-6-13
            zz=signo*sqrt(zz) ; signo=-1*signo                               ! miguel 4-6-13 
            if(mod(k-1,3).eq.0)then      ; b=xx ; c=yy ; a=zz                ! miguel 4-6-13
            else if(mod(k-1,3).eq.1)then ; a=xx ; c=yy ; b=zz                ! miguel 4-6-13
              else if(mod(k-1,3).eq.2)then ; a=xx ; b=yy ; c=zz ; end if       ! miguel 4-6-13      
            !a=(2*ran2(idum)-1d0)*de;b=(2*ran2(idum)-1d0)*de;c=(2*ran2(idum)-1d0)*de   ! miguel 4-6-13    
            node(ii)%x=node(ent)%x+a
            node(ii)%y=node(ent)%y+b
            node(ii)%z=node(ent)%z+c
            node(ii)%tipus=3 ; node(ii)%icel=jj
            cels(jj)%node(k)=ii
        end do
    end if

    do j=2,radicel                                                        ! "cicle"
        dx1=0.0 ; dx2=0.0 ; dx3=0.0 ; dy1=0.0 ; dy2=0.0 ; dy3=0.0
        do i=1,6                                                        ! "sector" of the hexagon
            if(i.eq.1)then  

                dx2=0!*cos(real(i)*60d0*rad)
                dy2=di*(real(j)-1d0)!*sin(real(i)*60d0*rad);print*,"dx2",dx2,"dy2",dy2
                dx1=di*(real(j)-1d0)*sin(-60d0*rad) ; dy1=di*(real(j)-1d0)*cos(-60d0*rad)!;print*,"dx1",dx1

            else
                hip=di*(real(j)-1d0)                                    ! hipotenusa
                dx1=dx2 ; dy1=dy2
                dx2=hip*sin(real(i-1)*60d0*rad) ; dy2=hip*cos(real(i-1)*60d0*rad)         
            end if           
            ii=ii+1
            jj=jj+1
!            print*,"cel perif, node central. ent:",ent,ii
            node(ii)%x=node(ont)%x+dx2 ; node(ii)%y=node(ont)%y+dy2 ; node(ii)%z=node(ont)%z
!            print*,"vertex i",i,"ii",ii,node(ii)%x,node(ii)%y,node(ii)%z
!            print*,i,"ii",node(ii)%x,"ent",node(ent)%x
            node(ii)%icel=jj ; node(ii)%tipus=3
            cels(jj)%node(1)=ii
            ent=ii
            if(radi>1)then
                signo=1                                                                  ! miguel 4-6-13
                do k=2,radi
                    ii=ii+1
!                            print*,"cel perif, node extern. ent:",ent,ii
                            call random_number(a)
                            der=de!(sqrt(de)*sqrt(de*a))                               ! miguel 4-6-13
122                         call random_number(a)
                            xx=der*(1d0-2*a) ;
                            call random_number(a)
                            yy=der*(1d0-2*a)            ! miguel 4-6-13
                            zz=(der**2)-(xx**2)-(yy**2)                                      ! miguel 4-6-13
                            if(zz.lt.0)then;goto 122;endif                                   ! miguel 4-6-13
                            zz=signo*sqrt(zz) ; signo=-1*signo                               ! miguel 4-6-13 
                            if(mod(k-1,3).eq.0)then      ; b=xx ; c=yy ; a=zz               ! miguel 4-6-13
                            else if(mod(k-1,3).eq.1)then ; a=xx ; c=yy ; b=zz               ! miguel 4-6-13
                              else if(mod(k-1,3).eq.2)then ; a=xx ; b=yy ; c=zz ; end if      ! miguel 4-6-13

!                    a=(2*ran2(idum)-1d0)*de;b=(2*ran2(idum)-1d0)*de;c=(2*ran2(idum)-1d0)*de    ! miguel 4-6-13    
                    node(ii)%x=node(ent)%x+a
                    node(ii)%y=node(ent)%y+b
                    node(ii)%z=node(ent)%z+c
                    node(ii)%tipus=3 ; node(ii)%icel=jj
                    cels(jj)%node(k)=ii
                end do
            end if

            if(j.gt.2)then                                              ! intermediate points
                dx3=dx2-dx1       ; dy3=dy2-dy1                         ! vectors which link "extreme" points
                dx3=dx3/real(j-1) ; dy3=dy3/real(j-1)                   ! sub-vectors                   
                do k=1,j-2                              
                    ii=ii+1
                    jj=jj+1
                    node(ii)%x=node(ont)%x+(dx1+(real(k)*dx3))
                    node(ii)%y=node(ont)%y+(dy1+(real(k)*dy3))   
                    node(ii)%z=node(ont)%z                                               
!                    print*,"intra i",i,"ii",ii,node(ii)%x,node(ii)%y,node(ii)%z
                    node(ii)%icel=jj ; node(ii)%tipus=3 ; cels(jj)%node(1)=ii
                    ent=ii
                    if(radi>1)then
                        signo=1                                                        ! miguel 4-6-13
                        do kk=2,radi
                            ii=ii+1
!                            print*,"cel perif, node extern. ent:",ent,ii
                            call random_number(a)
                            der=de!(sqrt(de)*sqrt(de*a))                               ! miguel 4-6-13
123                         call random_number(a)
                            xx=der*(1d0-2*a) ;
                            call random_number(a)
                            yy=der*(1d0-2*a)            ! miguel 4-6-13
                            zz=(der**2)-(xx**2)-(yy**2)                                      ! miguel 4-6-13
                            if(zz.lt.0)then;goto 123;endif                                   ! miguel 4-6-13
                            zz=signo*sqrt(zz) ; signo=-1*signo                               ! miguel 4-6-13 
                            if(mod(ii,3).eq.0)then      ; b=xx ; c=yy ; a=zz               ! miguel 4-6-13
                            else if(mod(ii,3).eq.1)then ; a=xx ; c=yy ; b=zz               ! miguel 4-6-13
                              else if(mod(ii,3).eq.2)then ; a=xx ; b=yy ; c=zz ; end if      ! miguel 4-6-13
                       
                            !a=(2*ran2(idum)-1d0)*de;b=(2*ran2(idum)-1d0)*de;c=(2*ran2(idum)-1d0)*de    ! miguel 4-6-13    
                            node(ii)%x=node(ent)%x+a
                            node(ii)%y=node(ent)%y+b
                            node(ii)%z=node(ent)%z+c
                            node(ii)%tipus=3 ; node(ii)%icel=jj
                            cels(jj)%node(kk)=ii
                        end do
                    end if
                end do
            end if

        end do
    end do
  end do




    !define the cell's centroid and nucleus (margin)            !>>>>>>>>>>>>>> Miquel 14-4-13

    do i=ncelsepi+1,ncels                                            !>>>>>>>>>>>>>> Miquel 14-4-13
        cels(i)%ctipus=3                                    !>>>>>>>>>>>>>> Miquel 14-4-13
        cels(i)%cex=0;cels(i)%cey=0;cels(i)%cez=0            !>>>>>>>>>>>>>> Miquel 14-4-13
!        print*,"i",i,"node",cels(i)%node(:),"ncels",ncels
        do j=1,cels(i)%nunodes                                !>>>>>>>>>>>>>> Miquel 14-4-13
            k=cels(i)%node(j)                                !>>>>>>>>>>>>>> Miquel 14-4-13
!            print*,"k",k,"nunodes",cels(i)%nunodes
            cels(i)%cex=cels(i)%cex+node(k)%x                !>>>>>>>>>>>>>> Miquel 14-4-13
            cels(i)%cey=cels(i)%cey+node(k)%y                !>>>>>>>>>>>>>> Miquel 14-4-13
            cels(i)%cez=cels(i)%cez+node(k)%z                !>>>>>>>>>>>>>> Miquel 14-4-13
        end do                                                !>>>>>>>>>>>>>> Miquel 14-4-13
        cels(i)%cex=cels(i)%cex/real(cels(i)%nunodes)        !>>>>>>>>>>>>>> Miquel 14-4-13
        cels(i)%cey=cels(i)%cey/real(cels(i)%nunodes)        !>>>>>>>>>>>>>> Miquel 14-4-13
        cels(i)%cez=cels(i)%cez/real(cels(i)%nunodes)        !>>>>>>>>>>>>>> Miquel 14-4-13
    end do                                                    !>>>>>>>>>>>>>> Miquel 14-4-13

    do i=ncelsepi+1,ncels
          b=1.0d8                                            !>>>>>>>>>>>>>> Is 14-9-13
         do j=1,cels(i)%nunodes                                !>>>>>>>>>>>>>> Is 14-9-13
            k=cels(i)%node(j)                                !>>>>>>>>>>>>>> Is 14-9-13
           a=sqrt((cels(i)%cex-node(k)%x)**2+(cels(i)%cey-node(k)%y)**2+(cels(i)%cez-node(k)%z)**2)
            if (b>a) then ; b=a ; ii=k ; jj=j ; end if
      end do                                                !>>>>>>>>>>>>>> Is 14-9-13
          node(ii)%marge=0

          !now we want the nucleus to be the first node in cels%node: this allows diffusion to be faster 
          jjj=cels(i)%node(1)
          cels(i)%node(1)=ii
          cels(i)%node(jj)=jjj
    end do                                                    !>>>>>>>>>>>>>> Is 14-9-13

    node(:)%talone=0.0d0

end subroutine mesenq1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111


!*********************************************************************************************************



subroutine epidermal_organ

!real*8::pericel,radio2,req1,req2,adjustreq,di
!integer::kk2

integer:: lx,ly
integer:: x1,x2,y1,y2,ii1,ii2,jj1,jj2
integer:: x3,x4,y3,y4,ii3,ii4,jj3,jj4,jj5,y5
integer, allocatable::cell_grid_epi(:,:),cell_grid_mes(:,:,:)
!integer:: periodic_borders,nborders
integer::geometry
integer::xlayer

!******* #1 DEFINING SPATIAL DIMENSIONS *******

  geometry=1 !geometry=1 means hexagonal epithelium ; geometry=2 means rectangular epithelium
  
  if (geometry==1)then
	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=10    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=10   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=-0.5d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    do i=1,layer
      mesradicel(i)=mradicel
    end do
    
    xlayer=0 !this will make a layer of ECM beneath the layers of mesenchymal cells
    
  elseif(geometry==2)then
    layer=1      !number of planar cell layers
    lx=32
    ly=32
    
    allocate(cell_grid_epi(lx,ly),cell_grid_mes(lx,ly,layer))

  end if  
    !periodic_borders=1
    


    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(geometry==1)then
      if(radi>0.and.radicel>0)then
        j=0
        do i=1,radi-1
          j=j+i
        end do
        nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
        nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13
  
        j=0
        do i=1,radicel-1
          j=j+i
        end do
        ncelsepi=(6*j+1)
        ndepi=nodecel*ncelsepi
      else
        ncelsepi=0
        ndepi=0
        nodecel=0
        nodecela=0
      end if
  
      if(mradi>0.and.mradicel>0.and.layer>0)then
        if(packed==1)then !if mesenchyme is packed, we use mradi differently
          j=0
          do i=1,mradi-1
            j=j+i
          end do
          nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
          nodecela=2*nodecel+1
        else
          nodecel=mradi
          nodecela=2*nodecel+1
        end if
  
        ncelsmes=0
        do k=1,layer
          j=0 !;print*,"mesradicel",mesradicel(k)
          do i=1,mesradicel(k)-1
            j=j+i
          end do
          ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
        end do
        ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
        if(radi==0.and.radicel==0)then
          nodecel=mradi
          nodecela=2*nodecel+1
        else if(nodecel<mradi)then
          nodecel=radi
          nodecela=nodecel*2+1
        end if
      else
        ndmes=0
        ncelsmes=0
      end if
      !ndx=(6*j+1)*xlayer
      nd=ndepi+ndmes!+ndx
      ncels=ncelsepi+ncelsmes
      nda=nd+10
      ncals=ncels+10
    elseif(geometry==2)then

      ncelsepi=lx*ly
      ndepi=ncelsepi*2
      ncelsmes=ncelsepi*layer
      ndmes=ncelsmes
      nd=ndepi+ndmes
      ncels=ncelsepi+ncelsmes
      nodecel=2 ; nodecela=5
      nda=nd+10
      ncals=ncels+10
    end if
    

  !End initializing dimension parameters

  print*,"inicial nd",nd,"ncels",ncels,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=1
    !End Allocatations


    

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    deltamin=1d-2
    khold=1d0
    angletor=0.05
    k_bu=5d0
    ramax=0.35d0
    k_press=5d-1
    
    mi_xwall=-0.5d0
    
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=1
    ffu(2)=0 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !buoyancy

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(14)=1 !physical boundaries (walls)
    ffu(17)=1 !volume conservation (epithelial)
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise
    ffu(24)=0


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    !if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(mod(i,2)/=0)then  !basal
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
          !node(i)%reqp=-0.05d0
          !node(i)%req=0.20d0
        !node(i)%kplast=1d0
        else                      !apical
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !node(i)%kplast=1d1
        !  !node(i)%reqp=0.05d0
        !  !node(i)%req=0.30d0
        end if
        node(i)%req=0.25d0
        node(i)%reqcr=0.25d0 !; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.25
        node(i)%da=node(i)%req*1.60; 
        node(i)%ke=5d-1
        node(i)%tor=3d0
        node(i)%stor=1d0
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=3d0
        node(i)%kvol=3d1
      end do
    !end if

	!Mesenchyme
    if(layer>0)then
      do i=ndepi+1,ndepi+ndmes
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=1d-1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.60
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

	!ECM
    if(xlayer>0)then
      do i=ndepi+ndmes+1,nd
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.60
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if
    
    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)

    node(:)%hold=0


  if(geometry==1)then
    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    !if(radi>0.and.radicel>0)               call epiteli_sphere1(radi,radicel,radius,di,zepi)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if


    !ECM basal layer
    !do i=ndepi+ndmes+1,nd
    !  j=i-ndx
    !  node(i)%x=node(j)%x ; node(i)%y=node(j)%y ; node(i)%z=node(j)%z-0.5d0
    !  node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    !  node(i)%icel=-i ; node(i)%tipus=4
    !  !node(i)%border=1
    !end do

      
    j=0
    do i=1,radicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    do i=j+1,ndepi
      node(i)%hold=2 ;node(i)%repcel=1d0
    end do
    j=0
    do i=1,mradicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    k=0
    do i=1,mradicel-1
      k=k+i
    end do
    k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
      node(i)%hold=2;node(i)%repcel=1d0
      !node(i)%border=1
      !node(i)%da=node(i)%req*2.0
      !node(i)%orix=0 ; node(i)%oriy=0
      !node(i)%rep=1d1;node(i)%repcel=1d1
    end do
    !do i=ndepi+ndmes/2+j+1,nd
    !  node(i)%hold=1;node(i)%repcel=1d2
    !  !node(i)%border=1
    !end do
    do i=ndepi+k+1,nd
      node(i)%hold=2;node(i)%repcel=1d0
      node(i)%border=0
    end do

    
  elseif(geometry==2)then
	do i=1,ncels
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do


    !epithelial cells
    !origin of rectangle
    aa=node(1)%req*2
    bb=sqrt(aa**2-(aa*0.5)**2)
    a=-real(lx)*0.5*aa ; b=-real(ly)*0.5*bb
    !print*,"a",a,"b",b,"aa",aa,"bb",bb,"lx",lx,"ly",ly
    ii=0 ; jj=0
    do i=1,ly
      do j=1,lx
        ii=ii+1 ;jj=jj+1 ; cell_grid_epi(j,i)=jj
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii+1 ;node(ii)%tipus=2
        ii=ii+1
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii-1 ;node(ii)%tipus=1
        cels(jj)%node(1)=ii-1 ; cels(jj)%node(2)=ii ;cels(jj)%nunodes=2; cels(jj)%ctipus=1
        if(i==1.or.i==ly.or.j==1.or.j==lx)then
          node(ii)%hold=2;node(ii-1)%hold=2 ;
          node(ii)%repcel=1d1;node(ii-1)%repcel=1d1;
          !node(ii)%border=1;node(ii-1)%border=1;
        end if
      end do
    end do

    !mesenchymals
    do k=1,layer
      !print*,"k",k
      do i=1,ly
        do j=1,lx
          ii=ii+1 ; jj=jj+1 ; cell_grid_mes(j,i,k)=jj
          if(mod(i,2)==0)then
            node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          else
            node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          end if
          node(ii)%icel=jj ; node(ii)%altre=0 ;node(ii)%tipus=3
          cels(jj)%node(1)=ii ;cels(jj)%nunodes=1 ; cels(jj)%ctipus=3
          if(i==1 .or. i==ly .or. j==1 .or. j==lx)then
            node(ii)%hold=2 ;node(ii)%repcel=1d1 ; !node(ii)%border=2 ;
          end if
          !if(i==1 .or. i==ly .or. j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
          !if(j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
        end do
      end do
    end do
  end if

  !m_xwall=maxval(node(1:nd)%x)+0.01
  !mi_xwall=minval(node(1:nd)%x)-0.01
  !m_ywall=maxval(node(1:nd)%y)+0.01
  !mi_ywall=minval(node(1:nd)%y)-0.01
  


    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
  !node(:)%hold=2
  !node(:)%hold=0  

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    !if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    !end if
    !Mesenchymal
    if(layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=20
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=2.5d0 ; gen(1)%mu=5d-1 ;gen(1)%name="epithelial suprabasal marker"
      gen(1)%npre=1 ; allocate(gen(1)%pre(gen(1)%npre))
      gen(1)%pre(1)=7

      
      gen(2)%kindof=9 ; gen(2)%diffu=4d-2 ; gen(2)%mu=2.0d-2 ;gen(2)%name="primary activator signal"
      gen(2)%npre=1 ; allocate(gen(2)%pre(gen(2)%npre))
      gen(2)%pre(1)=3 ;! gen(2)%pre(2)=5

      
      gen(3)%kindof=1 ; gen(3)%diffu=1d-1 ; gen(3)%mu=5d-1 ;gen(3)%name="epithelial basal marker"
      !gen(3)%npost=1 ; allocate(gen(3)%post(gen(3)%npost))
      !gen(3)%post(1)=4

      gen(4)%kindof=1 ; gen(4)%diffu=0d0 ; gen(4)%mu=5d-1 ;gen(4)%name="epithelial knot marker"
      !gen(4)%kindof=8 ; gen(4)%diffu=0d0 ; gen(4)%mu=5d-1 !
      !gen(4)%npre=2 ; allocate(gen(4)%pre(gen(4)%npre))
      !gen(4)%pre(1)=2 ; gen(4)%pre(2)=3
      !gen(4)%npost=2 ; allocate(gen(4)%post(gen(4)%npost))
      !gen(4)%post(1)=2 ; gen(4)%post(2)=3

      gen(5)%kindof=1 ; gen(5)%diffu=0d0 ; gen(5)%mu=5d-1 ;gen(5)%name="mesenchymal marker"
      !gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost))
      !gen(5)%post(1)=6

      gen(6)%kindof=1 ; gen(6)%diffu=0d0 ; gen(6)%mu=5.0d-1 ;gen(6)%name="mesenchymal growth zone marker"
      !gen(6)%kindof=8 ; gen(6)%diffu=0d0 ; gen(6)%mu=5.0d-1 !
      !gen(6)%npre=2 ; allocate(gen(6)%pre(gen(6)%npre))
      !gen(6)%pre(1)=2 ; gen(6)%pre(2)=5
      !gen(6)%npost=2 ; allocate(gen(6)%post(gen(6)%npost))
      !gen(6)%post(1)=2 ; gen(6)%post(2)=5

      gen(7)%kindof=1 ; gen(7)%diffu=1d0 ; gen(7)%mu=2.5d-1 ;gen(7)%name="epithelial placode marker"

      gen(8)%kindof=1 ; gen(8)%diffu=0d0 ; gen(8)%mu=5d-1 ;gen(8)%name="mesenchymal placode marker"

      gen(9)%kindof=9 ; gen(9)%diffu=2.75d0 ; gen(9)%mu=8.0d-2 ;gen(9)%name="primary inhibitor signal"
      gen(9)%npre=1 ; allocate(gen(9)%pre(gen(9)%npre))
      gen(9)%pre(1)=3 ;! gen(9)%pre(2)=5
      
      gen(10)%kindof=9 ; gen(10)%diffu=1d-1 ; gen(10)%mu=5.0d-1 ;gen(10)%name="secondary activator signal"
      gen(10)%npre=1 ; allocate(gen(10)%pre(gen(10)%npre))
      gen(10)%pre(1)=7 ; !gen(10)%pre(2)=8

      gen(11)%kindof=9 ; gen(11)%diffu=2d0 ; gen(11)%mu=1.0d0 ;gen(11)%name="secondary inhibitor signal"
      gen(11)%npre=1 ; allocate(gen(11)%pre(gen(11)%npre))
      gen(11)%pre(1)=7 ; !gen(11)%pre(2)=8
      !gen(11)%kindof=8 ; gen(11)%diffu=0d0 ; gen(11)%mu=5.0d-1 !
      !gen(11)%npre=2 ; allocate(gen(11)%pre(gen(11)%npre))
      !gen(11)%pre(1)=9 ; gen(11)%pre(2)=10
      !gen(11)%npost=2 ; allocate(gen(11)%post(gen(11)%npost))
      !gen(11)%post(1)=9 ; gen(11)%post(2)=10

      gen(12)%kindof=1 ; gen(12)%diffu=1d0 ; gen(12)%mu=2.5d-1 ;gen(12)%name="epithelial basal adhesion molecule"
      gen(13)%kindof=1 ; gen(13)%diffu=0d0 ; gen(13)%mu=5.0d-1 ;gen(13)%name="mesench adhesion molecule"

      gen(14)%kindof=1 ; gen(14)%diffu=0d0 ; gen(14)%mu=5.0d-1 ;gen(14)%name="epithelial interplacode marker"
      gen(15)%kindof=1 ; gen(15)%diffu=0d0 ; gen(15)%mu=5.0d-1  ;gen(15)%name="mesenchymal interplacode marker"
      gen(16)%kindof=1 ; gen(16)%diffu=0d0 ; gen(16)%mu=0.0d0  ;gen(16)%name="epithelial suprabasal adhesion molecule"
      
      gen(17)%kindof=9 ; gen(17)%diffu=1d-1 ; gen(17)%mu=5.0d-1  ;gen(17)%name="2ary activator-like, signal for mesenchyme"      
      gen(17)%npre=1 ; allocate(gen(17)%pre(gen(17)%npre))
      gen(17)%pre(1)=8 ; !gen(11)%pre(2)=8

      gen(18)%kindof=9 ; gen(18)%diffu=2d-1 ; gen(18)%mu=5.0d-1  ;gen(18)%name="1ary activator-like, signal for mesenchyme"      
      gen(18)%npre=1 ; allocate(gen(18)%pre(gen(18)%npre))
      gen(18)%pre(1)=5 ; !gen(11)%pre(2)=8
      
      gen(19)%kindof=1 ; gen(19)%diffu=0d0 ; gen(19)%mu=5.0d-1  ;gen(19)%name="epithelial growth zone marker"      

      gen(20)%kindof=4 ; gen(20)%diffu=1d-2 ; gen(20)%mu=5.0d-2  ;gen(20)%name="external inhibitor signal"      
      gen(20)%npre=2 ; allocate(gen(20)%pre(gen(20)%npre))
      gen(20)%pre(1)=3 ; gen(20)%pre(2)=5


    !Gene-behavior interactions
    

    !!!!STANDARDIZED CELL BEHAVIORS
    
      gen(12)%wa(1)=1 !adhesion molecules
      gen(13)%wa(1)=2
      gen(16)%wa(1)=3
      
      !gen(2)%wa(nparam_per_node+8)=1d0     !1ary activator as a polarizing cue
      gen(7)%wa(nparam_per_node+16)=0d0   !epithelial placode marker triggers polar cell migration

      gen(10)%wa(nparam_per_node+8)=1d0     !1ary activator as a polarizing cue
      gen(6)%wa(nparam_per_node+11)=1d2     !1ary activator as a polarizing cue
      gen(19)%wa(nparam_per_node+11)=1d2     !1ary activator as a polarizing cue
      
      !regulation of cell proliferation

      !gen(7)%wa(nparam_per_node+2)=-5d-2  !placode and interplacode markers form an inhibitor threshold to cell proliferation
      !gen(8)%wa(nparam_per_node+2)=-5d-2
      !gen(14)%wa(nparam_per_node+2)=-1d1
      !gen(15)%wa(nparam_per_node+2)=-1d1
    
      !gen(4)%wa(nparam_per_node+2)=5d-2
      gen(19)%wa(nparam_per_node+2)=2.0d-2
      gen(6)%wa(nparam_per_node+2)=1.0d-1
      gen(8)%wa(nparam_per_node+2)=1.0d-3
    !!!!!!!!!!! 
    
      !gen(12)%wa(1)=1
      !gen(13)%wa(1)=2
      !!gen(7)%wa(1)=3
      !!gen(8)%wa(1)=4
      !
      !gen(16)%wa(1)=3
      !
      !gen(2)%wa(nparam_per_node+8)=1d0
      !gen(15)%wa(nparam_per_node+16)=0d0
      !!gen(7)%wa(13)=0.0d0
      !gen(8)%wa(nparam_per_node+2)=2d-2
      !gen(14)%wa(nparam_per_node+2)=8d-3
      

      
    !Gene-gene interactions

     !wavefront setting

    !!!!STANDARDIZED GENE NETWORK

      gen(2)%w(2)=8d0  !Primary RD network module
      gen(9)%w(2)=1d0
      gen(2)%w(9)=-5d0

      gen(18)%w(7)=1d0 !Primary activator activates a signal that only communicates to the mesenchyme
      gen(8)%w(18)=5d0
      
      gen(7)%w(2)=1d0  !placode markers activated by primary RD
      
      !gen(2)%w(7)=1d0 !placode marker activates activator?
      !gen(8)%w(2)=1d0

      gen(7)%w(7)=1d1  !autoactivation of placode markers
      gen(8)%w(8)=1d1


      
      gen(7)%w(3)=-1d1  !extended inhibition allows for threshold-based expression of placode markers
      gen(8)%w(3)=-1d5
      gen(7)%w(5)=-1d5
      gen(8)%w(5)=-5d0
      gen(14)%w(7)=-1d2
      gen(15)%w(8)=-1d2
      gen(15)%w(3)=-1d5
      gen(14)%w(5)=-1d5
      
      !gen(7)%w(14)=-1d2
      !gen(8)%w(15)=-1d2


      gen(3)%w(3)=1d0  !autoactivation of tissue markers
      gen(5)%w(5)=1d0
      
      gen(12)%w(12)=1d0 !autoactivation of adhesion molecules
      gen(13)%w(13)=1d0

      gen(14)%w(3)=1d1 !autoactivation of interplacode markers
      gen(15)%w(5)=1d1

      !gen(12)%w(7)=-2d0
      !gen(13)%w(8)=-2d0
      
      gen(10)%w(8)=1d-1 !placode marker activates secondary RD module

      gen(10)%w(10)=2d1 !Secondary RD network module
      gen(11)%w(10)=1d0
      gen(10)%w(11)=-5d0
      
      !gen(20)%w(14)=1d0
      !gen(10)%w(20)=-5d0

      gen(20)%w(20)=1d0
      gen(7)%w(20)=-1d1
      gen(8)%w(20)=-1d1


      gen(10)%w(14)=-1d5 !secondary RD signals not allowed outside the placode
      gen(11)%w(14)=-1d5
      gen(10)%w(15)=-1d5
      gen(11)%w(15)=-1d5


      gen(17)%w(4)=1d0 !Secondary activator activates a signal that only communicates to the mesenchyme

      
      gen(4)%w(10)=4d0 !epithelial knot marker activated by 2ari activator
      gen(4)%w(7)=-1d1 !basal inhibitor threshold
      
      gen(19)%w(11)=1d0 !epithelial GZ marker activated by 2ari inhibitor
      gen(19)%w(4)=-1d5 !epithelial knot marker excludes GZ marker
      gen(19)%w(3)=-3d1 !epithelial GZ marker inhibition threshold
 
      !gen(4)%w(7)=-1d-1 !basal inhibitor threshold
      
      
      
      gen(6)%w(17)=3d1 !mesench GZ marker activated by 2ari activator
      

      gen(4)%w(7)=-1d1 !GZ markers threshold-based inhibition
      gen(6)%w(8)=-5d0 !
      
      gen(4)%w(8)=-1d5 !epi GZ marker not allowed on mesenchyme
      gen(6)%w(7)=-1d5 !mesench GZ marker not allowed on epithelium
      gen(19)%w(8)=-1d5
      
      
      gen(4)%w(14)=-1d5 !GZ markers not allowed outside the placode
      gen(4)%w(15)=-1d5 
      gen(6)%w(14)=-1d5 
      gen(6)%w(15)=-1d5 



      
      !gen(10)%w(10)=5d0 !Secondary RD network module
      !gen(1)%w(10)=1d0
      !gen(1)%w(1)=-1d1
      !gen(1)%w(8)=-1d5
      !gen(10)%w(7)=-1d5
      
      !gen(11)%w(10)=1d0
      !gen(11)%w(11)=-5d0
      !gen(11)%w(8)=-1d5


      !gen(10)%w(1)=-1d2
      !gen(11)%w(1)=-1d2
      
      
      
      !gen(1)%w(14)=1d-1 !external inhibitor of the 2ari RD module (something like EDA?)
      !gen(10)%w(1)=-1d1
      !gen(11)%w(1)=-1d1
      
      
      !gen(4)%w(1)=1d0  !growth zone markers (tipus knot)
      !gen(6)%w(10)=1d0

      !gen(4)%w(10)=-1d0
      !gen(6)%w(10)=-1d0

      !gen(6)%w(7)=-1d2
      !gen(4)%w(8)=-1d2
      !gen(6)%w(14)=-1d2
      !gen(4)%w(15)=-1d2      
      !gen(6)%w(15)=-1d2
      !gen(4)%w(14)=-1d2      

      


    !Adhesion molecules

	ntipusadh=3
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=5d-1
      kadh(1,2)=1d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=5d-1
      kadh(1,3)=5d-1 ; kadh(3,1)=kadh(1,3)
      kadh(3,3)=5d-1
      !kadh(3,4)=1d0 ; kadh(4,3)=kadh(3,4)
      !kadh(4,2)=5d-1 ; kadh(4,2)=kadh(2,4)
      !kadh(4,4)=3d0
      !kadh(5,5)=5d0
      !kadh(5,2)=5d-1 ; kadh(2,5)=kadh(5,2)
      !kadh(5,4)=5d-1 ; kadh(4,5)=kadh(5,4)

    end if

    !Gene expression on nodes

    
    if(geometry==1)then
      gex(1,2)=1d0
      !gex(ndepi+1,2)=1d0
    end if
    
    if(geometry==2)then
    !draw rectangles
    !!!!!!!!!!!!!!!!!!!!!!!!!
    jj1=1 ; jj2=ly
    ii1=1 ; ii2=lx
    do i=ii1,ii2
      do j=jj1,jj2
        k=cels(cell_grid_epi(i,j))%node(1)
        if(i<=2.or.i>=lx-1.or.j<=2.or.j>=ly-1)then
          gex(k,1)=1d0
          node(k)%hold=3
          node(node(k)%altre)%hold=3
          !gex(k,2)=0d0 ; gex(k,9)=0d0
        end if
        do ii=1,layer
          k=cels(cell_grid_mes(i,j,ii))%node(1)
          if(i<=2.or.i>=lx-1.or.j<=2.or.j>=ly-1)then
            gex(k,1)=1d0
            node(k)%hold=3
            !gex(k,2)=0d0 ; gex(k,9)=0d0
          end if
        end do
      end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!
    end if

   
    do i=1,nd
    
      if(node(i)%tipus<3)then
        if(node(i)%tipus==2)then
    
          !if(node(i)%x<-0.2d0) gex(i,1)=1d0 !AP regionalization
    
          gex(i,3)=1d0  !receptors
          !gex(i,10)=1d0
          gex(i,14)=1d0
          !call random_number(a)
          !!!d=sqrt(node(i)%x**2+node(i)%y**2)
          !gex(i,2)=a!/(1+d)!*0.1d0
          !gex(i,9)=a!/(1+d)!*0.1d0
        end if
        !if(node(i)%tipus==1) gex(i,14)=1d0  !receptors
       
        gex(i,12)=1d0
        !gex(i,10)=1d0

      if(node(i)%hold>0) gex(i,20)=1d0

      end if
      if(node(i)%tipus==3)then
        !if(node(i)%x<-0.2d0) gex(i,1)=1d0 !AP regionalization
        gex(i,5)=1d0
        gex(i,13)=1d0
        gex(i,15)=1d0
          !call random_number(a)
          !!!d=sqrt(node(i)%x**2+node(i)%y**2)
          !gex(i,2)=a!/(1+d)!*0.1d0
          !gex(i,9)=a!/(1+d)!*0.1d0

      if(node(i)%hold>0) gex(i,20)=1d0
 
     end if
      if(node(i)%tipus==4)then
        gex(i,16)=1d0
        !print*,"node border", node(i)%border
      end if
      !if(node(i)%hold>0)then
      !  gex(i,:)=0d0 !clear the hold nodes
      !  if(node(i)%tipus<3)then
      !    gex(i,12)=1d0
      !  end if
      !  if(node(i)%tipus==3)then
      !    gex(i,13)=1d0
      !  end if
      !end if

    end do


    call update_npag

    node(:)%talone=0.0d0

    !ramax=maxval(node(:)%da)*3

    
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine epidermal_organ_in

!real*8::pericel,radio2,req1,req2,adjustreq,di
!integer::kk2

integer:: lx,ly
integer:: x1,x2,y1,y2,ii1,ii2,jj1,jj2
integer:: x3,x4,y3,y4,ii3,ii4,jj3,jj4,jj5,y5
integer, allocatable::cell_grid_epi(:,:),cell_grid_mes(:,:,:)
!integer:: periodic_borders,nborders
integer::geometry
integer::xlayer

!******* #1 DEFINING SPATIAL DIMENSIONS *******

  geometry=1 !geometry=1 means hexagonal epithelium ; geometry=2 means rectangular epithelium
  
  if (geometry==1)then
	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=10    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=10   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=1.0d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    do i=1,layer
      mesradicel(i)=mradicel
    end do
    
    xlayer=1 !this will make a layer of ECM beneath the layers of mesenchymal cells
    
  elseif(geometry==2)then
    layer=1      !number of planar cell layers
    lx=32
    ly=32
    
    allocate(cell_grid_epi(lx,ly),cell_grid_mes(lx,ly,layer))

  end if  
    !periodic_borders=1
    


    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(geometry==1)then
      if(radi>0.and.radicel>0)then
        j=0
        do i=1,radi-1
          j=j+i
        end do
        nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
        nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13
  
        j=0
        do i=1,radicel-1
          j=j+i
        end do
        ncelsepi=(6*j+1)
        ndepi=nodecel*ncelsepi
      else
        ncelsepi=0
        ndepi=0
        nodecel=0
        nodecela=0
      end if
  
      if(mradi>0.and.mradicel>0.and.layer>0)then
        if(packed==1)then !if mesenchyme is packed, we use mradi differently
          j=0
          do i=1,mradi-1
            j=j+i
          end do
          nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
          nodecela=2*nodecel+1
        else
          nodecel=mradi
          nodecela=2*nodecel+1
        end if
  
        ncelsmes=0
        do k=1,layer
          j=0 !;print*,"mesradicel",mesradicel(k)
          do i=1,mesradicel(k)-1
            j=j+i
          end do
          ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
        end do
        ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
        if(radi==0.and.radicel==0)then
          nodecel=mradi
          nodecela=2*nodecel+1
        else if(nodecel<mradi)then
          nodecel=radi
          nodecela=nodecel*2+1
        end if
      else
        ndmes=0
        ncelsmes=0
      end if
      ndx=(6*j+1)*xlayer
      nd=ndepi+ndmes+ndx
      !ncels=ncelsepi+ncelsmes
      ncels=ncelsepi+ncelsmes+ndx
   print*,"ncels",ncels,ncelsepi,ncelsmes,ndx
      nda=nd+10
      ncals=ncels+10
    elseif(geometry==2)then

      ncelsepi=lx*ly
      ndepi=ncelsepi*2
      ncelsmes=ncelsepi*layer
      ndmes=ncelsmes
      nd=ndepi+ndmes
      ncels=ncelsepi+ncelsmes
      nodecel=2 ; nodecela=5
      nda=nd+10
      ncals=ncels+10
    end if
    

  !End initializing dimension parameters

  print*,"inicial nd",nd,"ncels",ncels,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=1
    !End Allocatations


    

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    deltamin=1d-2
    khold=1d0
    angletor=0.05
    k_bu=5d0
    ramax=0.35d0
    k_press=5d-1
    m_xwall=1d0
  
    
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=1
    ffu(2)=0 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !buoyancy

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(14)=1 !physical boundaries (walls)
    ffu(17)=1 !volume conservation (epithelial)
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise
    ffu(24)=0


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    !if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(mod(i,2)/=0)then  !basal
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
          !node(i)%reqp=-0.05d0
          !node(i)%req=0.20d0
        !node(i)%kplast=1d0
        else                      !apical
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !node(i)%kplast=1d1
        !  !node(i)%reqp=0.05d0
        !  !node(i)%req=0.30d0
        end if
        node(i)%req=0.25d0
        node(i)%reqcr=0.25d0 !; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.25
        node(i)%da=node(i)%req*1.60; 
        node(i)%ke=5d-1
        node(i)%tor=3d0
        node(i)%stor=1d0
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=3d0
        node(i)%kvol=3d1
      end do
    !end if

	!Mesenchyme
    if(layer>0)then
      do i=ndepi+1,ndepi+ndmes
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.60
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

	!ECM
    if(xlayer>0)then
      do i=ndepi+ndmes+1,nd
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.60
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if
    
    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)

    node(:)%hold=0


  if(geometry==1)then
    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    !if(radi>0.and.radicel>0)               call epiteli_sphere1(radi,radicel,radius,di,zepi)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if
!print*,"tall0"

    !Mesenchymal layer
    ii=0
    do i=ndepi+ndmes+1,nd
      ii=ii+1
      j=i-ndx
      node(i)%x=node(j)%x ; node(i)%y=node(j)%y ; node(i)%z=-0.5d0
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
      node(i)%icel=node(j)%icel+ndx ; node(i)%tipus=3
      !node(i)%icel=ncelsepi+ncelsmes+ii ; node(i)%tipus=3

      node(i)%repcel=5d-1
      !allocate(cels(node(i)%icel)%node(nodela))
      cels(node(i)%icel)%nunodes=1
      cels(node(i)%icel)%node(1)=i
      node(i)%marge=0
      !print*,"i",i,"j",j,"icel i",node(i)%icel,"icel j",node(j)%icel
      !node(i)%border=1
    end do
    ncelsmes=ncelsmes+ndx ; ndx=0
    ndmes=ncelsmes
!print*,"tall1"
      
    j=0
    do i=1,radicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    do i=j+1,ndepi
      node(i)%hold=2 ;node(i)%repcel=1d0
    end do
    j=0
    do i=1,mradicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    k=0
    do i=1,mradicel-1
      k=k+i
    end do
    k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
      node(i)%hold=2;node(i)%repcel=1d0
      !node(i)%border=1
      !node(i)%da=node(i)%req*2.0
      !node(i)%orix=0 ; node(i)%oriy=0
      !node(i)%rep=1d1;node(i)%repcel=1d1
    end do
    do i=ndepi+ndmes/2+j+1,nd
      node(i)%hold=2;node(i)%repcel=1d0
      !node(i)%border=1
    end do
    !do i=ndepi+k+1,nd
    !  node(i)%hold=2;node(i)%repcel=1d0
    !  node(i)%border=1
    !end do

    
  elseif(geometry==2)then
	do i=1,ncels
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do


    !epithelial cells
    !origin of rectangle
    aa=node(1)%req*2
    bb=sqrt(aa**2-(aa*0.5)**2)
    a=-real(lx)*0.5*aa ; b=-real(ly)*0.5*bb
    !print*,"a",a,"b",b,"aa",aa,"bb",bb,"lx",lx,"ly",ly
    ii=0 ; jj=0
    do i=1,ly
      do j=1,lx
        ii=ii+1 ;jj=jj+1 ; cell_grid_epi(j,i)=jj
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii+1 ;node(ii)%tipus=2
        ii=ii+1
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii-1 ;node(ii)%tipus=1
        cels(jj)%node(1)=ii-1 ; cels(jj)%node(2)=ii ;cels(jj)%nunodes=2; cels(jj)%ctipus=1
        if(i==1.or.i==ly.or.j==1.or.j==lx)then
          node(ii)%hold=2;node(ii-1)%hold=2 ;
          node(ii)%repcel=1d1;node(ii-1)%repcel=1d1;
          !node(ii)%border=1;node(ii-1)%border=1;
        end if
      end do
    end do

    !mesenchymals
    do k=1,layer
      !print*,"k",k
      do i=1,ly
        do j=1,lx
          ii=ii+1 ; jj=jj+1 ; cell_grid_mes(j,i,k)=jj
          if(mod(i,2)==0)then
            node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          else
            node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          end if
          node(ii)%icel=jj ; node(ii)%altre=0 ;node(ii)%tipus=3
          cels(jj)%node(1)=ii ;cels(jj)%nunodes=1 ; cels(jj)%ctipus=3
          if(i==1 .or. i==ly .or. j==1 .or. j==lx)then
            node(ii)%hold=2 ;node(ii)%repcel=1d1 ; !node(ii)%border=2 ;
          end if
          !if(i==1 .or. i==ly .or. j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
          !if(j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
        end do
      end do
    end do
  end if

  !m_xwall=maxval(node(1:nd)%x)+0.01
  !mi_xwall=minval(node(1:nd)%x)-0.01
  !m_ywall=maxval(node(1:nd)%y)+0.01
  !mi_ywall=minval(node(1:nd)%y)-0.01
  
!print*,"tall2"

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
  !node(:)%hold=2
  !node(:)%hold=0  

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    !if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    !end if
    !Mesenchymal
    if(layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=21
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=2.5d0 ; gen(1)%mu=5d-1 ;gen(1)%name="epithelial suprabasal marker"
      gen(1)%npre=1 ; allocate(gen(1)%pre(gen(1)%npre))
      gen(1)%pre(1)=7

      
      gen(2)%kindof=9 ; gen(2)%diffu=4d-2 ; gen(2)%mu=2.0d-2 ;gen(2)%name="primary activator signal"
      gen(2)%npre=2 ; allocate(gen(2)%pre(gen(2)%npre))
      gen(2)%pre(1)=3 ; gen(2)%pre(2)=1

      
      gen(3)%kindof=1 ; gen(3)%diffu=1d-1 ; gen(3)%mu=5d-1 ;gen(3)%name="epithelial basal marker"
      !gen(3)%npost=1 ; allocate(gen(3)%post(gen(3)%npost))
      !gen(3)%post(1)=4

      gen(4)%kindof=1 ; gen(4)%diffu=0d0 ; gen(4)%mu=5d-1 ;gen(4)%name="epithelial knot marker"
      !gen(4)%kindof=8 ; gen(4)%diffu=0d0 ; gen(4)%mu=5d-1 !
      !gen(4)%npre=2 ; allocate(gen(4)%pre(gen(4)%npre))
      !gen(4)%pre(1)=2 ; gen(4)%pre(2)=3
      !gen(4)%npost=2 ; allocate(gen(4)%post(gen(4)%npost))
      !gen(4)%post(1)=2 ; gen(4)%post(2)=3

      gen(5)%kindof=1 ; gen(5)%diffu=0d0 ; gen(5)%mu=5d-1 ;gen(5)%name="mesenchymal marker"
      !gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost))
      !gen(5)%post(1)=6

      gen(6)%kindof=1 ; gen(6)%diffu=0d0 ; gen(6)%mu=5.0d-1 ;gen(6)%name="mesenchymal growth zone marker"
      !gen(6)%kindof=8 ; gen(6)%diffu=0d0 ; gen(6)%mu=5.0d-1 !
      !gen(6)%npre=2 ; allocate(gen(6)%pre(gen(6)%npre))
      !gen(6)%pre(1)=2 ; gen(6)%pre(2)=5
      !gen(6)%npost=2 ; allocate(gen(6)%post(gen(6)%npost))
      !gen(6)%post(1)=2 ; gen(6)%post(2)=5

      gen(7)%kindof=1 ; gen(7)%diffu=1d0 ; gen(7)%mu=2.5d-1 ;gen(7)%name="epithelial placode marker"

      gen(8)%kindof=1 ; gen(8)%diffu=0d0 ; gen(8)%mu=5d-1 ;gen(8)%name="mesenchymal placode marker"

      gen(9)%kindof=9 ; gen(9)%diffu=2.75d0 ; gen(9)%mu=8.0d-2 ;gen(9)%name="primary inhibitor signal"
      gen(9)%npre=2 ; allocate(gen(9)%pre(gen(9)%npre))
      gen(9)%pre(1)=3 ; gen(9)%pre(2)=1
      
      gen(10)%kindof=9 ; gen(10)%diffu=1d-1 ; gen(10)%mu=5.0d-1 ;gen(10)%name="secondary activator signal"
      gen(10)%npre=1 ; allocate(gen(10)%pre(gen(10)%npre))
      gen(10)%pre(1)=7 ; !gen(10)%pre(2)=8

      gen(11)%kindof=9 ; gen(11)%diffu=2d0 ; gen(11)%mu=1.0d0 ;gen(11)%name="secondary inhibitor signal"
      gen(11)%npre=1 ; allocate(gen(11)%pre(gen(11)%npre))
      gen(11)%pre(1)=7 ; !gen(11)%pre(2)=8
      !gen(11)%kindof=8 ; gen(11)%diffu=0d0 ; gen(11)%mu=5.0d-1 !
      !gen(11)%npre=2 ; allocate(gen(11)%pre(gen(11)%npre))
      !gen(11)%pre(1)=9 ; gen(11)%pre(2)=10
      !gen(11)%npost=2 ; allocate(gen(11)%post(gen(11)%npost))
      !gen(11)%post(1)=9 ; gen(11)%post(2)=10

      gen(12)%kindof=1 ; gen(12)%diffu=1d0 ; gen(12)%mu=2.5d-1 ;gen(12)%name="epithelial basal adhesion molecule"
      gen(13)%kindof=1 ; gen(13)%diffu=0d0 ; gen(13)%mu=5.0d-1 ;gen(13)%name="mesench adhesion molecule"

      gen(14)%kindof=1 ; gen(14)%diffu=0d0 ; gen(14)%mu=5.0d-1 ;gen(14)%name="epithelial interplacode marker"
      gen(15)%kindof=1 ; gen(15)%diffu=0d0 ; gen(15)%mu=5.0d-1  ;gen(15)%name="mesenchymal interplacode marker"
      gen(16)%kindof=1 ; gen(16)%diffu=0d0 ; gen(16)%mu=0.0d0  ;gen(16)%name="epithelial suprabasal adhesion molecule"
      
      gen(17)%kindof=9 ; gen(17)%diffu=1d-1 ; gen(17)%mu=5.0d-1  ;gen(17)%name="2ary activator-like, signal for mesenchyme"      
      gen(17)%npre=1 ; allocate(gen(17)%pre(gen(17)%npre))
      gen(17)%pre(1)=8 ; !gen(11)%pre(2)=8

      gen(18)%kindof=9 ; gen(18)%diffu=2d-1 ; gen(18)%mu=5.0d-1  ;gen(18)%name="1ary activator-like, signal for mesenchyme"      
      gen(18)%npre=1 ; allocate(gen(18)%pre(gen(18)%npre))
      gen(18)%pre(1)=5 ; !gen(11)%pre(2)=8
      
      gen(19)%kindof=1 ; gen(19)%diffu=0d0 ; gen(19)%mu=5.0d-1  ;gen(19)%name="epithelial growth zone marker"      

      gen(20)%kindof=4 ; gen(20)%diffu=1d-2 ; gen(20)%mu=5.0d-2  ;gen(20)%name="external inhibitor signal"      
      gen(20)%npre=3 ; allocate(gen(20)%pre(gen(20)%npre))
      gen(20)%pre(1)=3 ; gen(20)%pre(2)=5 ; gen(20)%pre(3)=1

      gen(21)%kindof=1 ; gen(21)%diffu=0d0 ; gen(21)%mu=5.0d-1  ;gen(21)%name="placode basal side molecule"      

    !Gene-behavior interactions
    

    !!!!STANDARDIZED CELL BEHAVIORS
    
      gen(12)%wa(1)=1 !adhesion molecules
      gen(13)%wa(1)=2
      gen(16)%wa(1)=3
      
      gen(2)%wa(nparam_per_node+8)=1d0     !1ary activator as a polarizing cue
      gen(7)%wa(nparam_per_node+16)=5d-2   !epithelial placode marker triggers polar cell migration

      gen(10)%wa(nparam_per_node+8)=1d0     !1ary activator as a polarizing cue
      gen(6)%wa(nparam_per_node+11)=1d2     !1ary activator as a polarizing cue
      gen(19)%wa(nparam_per_node+11)=1d2     !1ary activator as a polarizing cue
      
      !regulation of cell proliferation

      !gen(7)%wa(nparam_per_node+2)=-5d-2  !placode and interplacode markers form an inhibitor threshold to cell proliferation
      !gen(8)%wa(nparam_per_node+2)=-5d-2
      !gen(14)%wa(nparam_per_node+2)=-1d1
      !gen(15)%wa(nparam_per_node+2)=-1d1
    
      !gen(4)%wa(nparam_per_node+2)=5d-2
      !gen(19)%wa(nparam_per_node+2)=2.0d-2
      !gen(6)%wa(nparam_per_node+2)=1.0d-1
      !gen(8)%wa(nparam_per_node+2)=1.0d-3
      gen(14)%wa(nparam_per_node+2)=5d-3
      
      
      gen(21)%wa(27)=-1d2
      gen(21)%wa(28)=-1d2
      gen(21)%wa(21)=5d-2
    !!!!!!!!!!! 
    
      !gen(12)%wa(1)=1
      !gen(13)%wa(1)=2
      !!gen(7)%wa(1)=3
      !!gen(8)%wa(1)=4
      !
      !gen(16)%wa(1)=3
      !
      !gen(2)%wa(nparam_per_node+8)=1d0
      !gen(15)%wa(nparam_per_node+16)=0d0
      !!gen(7)%wa(13)=0.0d0
      !gen(8)%wa(nparam_per_node+2)=2d-2
      !gen(14)%wa(nparam_per_node+2)=8d-3
      

      
    !Gene-gene interactions

     !wavefront setting

    !!!!STANDARDIZED GENE NETWORK

      gen(2)%w(2)=8d0  !Primary RD network module
      gen(9)%w(2)=1d0
      gen(2)%w(9)=-5d0

      gen(18)%w(7)=1d0 !Primary activator activates a signal that only communicates to the mesenchyme
      gen(8)%w(18)=5d0
      
      gen(7)%w(2)=1d0  !placode markers activated by primary RD
      
      !gen(2)%w(7)=1d0 !placode marker activates activator?
      !gen(8)%w(2)=1d0

      gen(7)%w(7)=1d1  !autoactivation of placode markers
      gen(8)%w(8)=1d1

      gen(21)%w(7)=1d0

      
      gen(7)%w(3)=-1d1  !extended inhibition allows for threshold-based expression of placode markers
      gen(7)%w(1)=-1d1  !extended inhibition allows for threshold-based expression of placode markers
      gen(8)%w(3)=-1d5
      gen(8)%w(1)=-1d5
      gen(7)%w(5)=-1d5
      gen(8)%w(5)=-5d0
      gen(14)%w(7)=-1d2
      gen(15)%w(8)=-1d2
      gen(15)%w(3)=-1d5
      gen(15)%w(1)=-1d5
      gen(14)%w(5)=-1d5
      
      !gen(7)%w(14)=-1d2
      !gen(8)%w(15)=-1d2

      gen(1)%w(1)=1d0  !autoactivation of tissue markers
      gen(3)%w(3)=1d0  !autoactivation of tissue markers
      gen(5)%w(5)=1d0
      
      gen(12)%w(12)=1d0 !autoactivation of adhesion molecules
      gen(13)%w(13)=1d0

      gen(14)%w(3)=1d1 !autoactivation of interplacode markers
      gen(14)%w(1)=1d1 !autoactivation of interplacode markers
      gen(15)%w(5)=1d1

      !gen(12)%w(7)=-2d0
      !gen(13)%w(8)=-2d0
      
      gen(10)%w(8)=1d-1 !placode marker activates secondary RD module

      gen(10)%w(10)=2d1 !Secondary RD network module
      gen(11)%w(10)=1d0
      gen(10)%w(11)=-5d0
      
      !gen(20)%w(14)=1d0
      !gen(10)%w(20)=-5d0

      gen(20)%w(20)=1d0
      gen(7)%w(20)=-1d1
      gen(8)%w(20)=-1d1


      gen(10)%w(14)=-1d5 !secondary RD signals not allowed outside the placode
      gen(11)%w(14)=-1d5
      gen(10)%w(15)=-1d5
      gen(11)%w(15)=-1d5


      gen(17)%w(4)=1d0 !Secondary activator activates a signal that only communicates to the mesenchyme

      
      gen(4)%w(10)=4d0 !epithelial knot marker activated by 2ari activator
      gen(4)%w(7)=-1d1 !basal inhibitor threshold
      
      gen(19)%w(11)=1d0 !epithelial GZ marker activated by 2ari inhibitor
      gen(19)%w(4)=-1d5 !epithelial knot marker excludes GZ marker
      gen(19)%w(3)=-3d1 !epithelial GZ marker inhibition threshold
      gen(19)%w(1)=-3d1 !epithelial GZ marker inhibition threshold
 
      !gen(4)%w(7)=-1d-1 !basal inhibitor threshold
      
      
      
      gen(6)%w(17)=3d1 !mesench GZ marker activated by 2ari activator
      

      gen(4)%w(7)=-1d1 !GZ markers threshold-based inhibition
      gen(6)%w(8)=-5d0 !
      
      gen(4)%w(8)=-1d5 !epi GZ marker not allowed on mesenchyme
      gen(6)%w(7)=-1d5 !mesench GZ marker not allowed on epithelium
      gen(19)%w(8)=-1d5
      
      
      gen(4)%w(14)=-1d5 !GZ markers not allowed outside the placode
      gen(4)%w(15)=-1d5 
      gen(6)%w(14)=-1d5 
      gen(6)%w(15)=-1d5 



      !!gen(1)%w(1)=1d0
      !!gen(2)%w(1)=-1d5
      !
      !gen(2)%w(2)=8d0  !Primary RD network module
      !gen(9)%w(2)=1d0
      !gen(2)%w(9)=-5d0
      !
      !gen(7)%w(2)=7d-1
      !gen(8)%w(2)=1d0
      !
      !gen(7)%w(3)=-1d1
      !gen(8)%w(3)=-1d3
      !gen(7)%w(5)=-1d3
      !gen(8)%w(5)=-1d1
      !gen(14)%w(7)=-1d2
      !gen(15)%w(8)=-1d2
      !
      !gen(1)%w(1)=1d0
      !gen(3)%w(3)=1d0  !autoactivation of tissue markers
      !gen(5)%w(5)=1d0
      !!gen(10)%w(10)=1d0
      !
      !gen(12)%w(12)=1d0 !autoactivation of adhesion molecules
      !gen(13)%w(13)=1d0
      !gen(16)%w(16)=1d0
      !
      !gen(14)%w(3)=1d0
      !gen(15)%w(5)=1d0
      !
      !!gen(4)%w(1)=-1d3
      !!gen(4)%w(7)=1d0
      !
      !!gen(6)%w(1)=-1d3
      !!gen(6)%w(8)=1d0
      !
      !gen(12)%w(7)=-2d0
      !gen(13)%w(8)=-2d0
      !
      !gen(10)%w(7)=1d-3
      !!gen(10)%w(4)=5d0
      !gen(10)%w(10)=2d1 !Secondary RD network module
      !gen(11)%w(10)=1d0
      !gen(10)%w(11)=-5d0
      !
      !
      !!gen(10)%w(10)=5d0 !Secondary RD network module
      !!gen(1)%w(10)=1d0
      !!gen(1)%w(1)=-1d1
      !!gen(1)%w(8)=-1d5
      !!gen(10)%w(7)=-1d5
      !
      !!gen(11)%w(10)=1d0
      !!gen(11)%w(11)=-5d0
      !!gen(11)%w(8)=-1d5
      !
      !
      !!gen(10)%w(1)=-1d2
      !!gen(11)%w(1)=-1d2
      !
      !gen(10)%w(14)=-1d2 !secondary RD signals not allowed outside the placode
      !gen(11)%w(14)=-1d2
      !gen(10)%w(15)=-1d2
      !gen(11)%w(15)=-1d2
      !
      !
      !!gen(1)%w(14)=1d-1 !external inhibitor of the 2ari RD module (something like EDA?)
      !!gen(10)%w(1)=-1d1
      !!gen(11)%w(1)=-1d1
      !
      !
      !gen(4)%w(1)=1d0  !growth zone markers (tipus knot)
      !gen(6)%w(10)=1d0
      !
      !!gen(4)%w(10)=-1d0
      !!gen(6)%w(10)=-1d0
      !
      !gen(6)%w(7)=-1d2
      !gen(4)%w(8)=-1d2
      !gen(6)%w(14)=-1d2
      !gen(4)%w(15)=-1d2      
      !gen(6)%w(15)=-1d2
      !gen(4)%w(14)=-1d2      

      


    !Adhesion molecules

	ntipusadh=3
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=1d0
      kadh(1,2)=5d-1 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=5d-1
      kadh(1,3)=5d-1 ; kadh(3,1)=kadh(1,3)
      kadh(3,3)=1d0
      !kadh(3,4)=1d0 ; kadh(4,3)=kadh(3,4)
      !kadh(4,2)=5d-1 ; kadh(4,2)=kadh(2,4)
      !kadh(4,4)=3d0
      !kadh(5,5)=5d0
      !kadh(5,2)=5d-1 ; kadh(2,5)=kadh(5,2)
      !kadh(5,4)=5d-1 ; kadh(4,5)=kadh(5,4)

    end if

    !Gene expression on nodes

    
    if(geometry==1)then
      gex(1,2)=1d0
      gex(ndepi+1,2)=1d0
      !gex(ndepi+1,2)=1d0
    end if
    
    if(geometry==2)then
    !draw rectangles
    !!!!!!!!!!!!!!!!!!!!!!!!!
    jj1=1 ; jj2=ly
    ii1=1 ; ii2=lx
    do i=ii1,ii2
      do j=jj1,jj2
        k=cels(cell_grid_epi(i,j))%node(1)
        if(i<=2.or.i>=lx-1.or.j<=2.or.j>=ly-1)then
          gex(k,1)=1d0
          node(k)%hold=3
          node(node(k)%altre)%hold=3
          !gex(k,2)=0d0 ; gex(k,9)=0d0
        end if
        do ii=1,layer
          k=cels(cell_grid_mes(i,j,ii))%node(1)
          if(i<=2.or.i>=lx-1.or.j<=2.or.j>=ly-1)then
            gex(k,1)=1d0
            node(k)%hold=3
            !gex(k,2)=0d0 ; gex(k,9)=0d0
          end if
        end do
      end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!
    end if

   
    do i=1,nd
    
      if(node(i)%tipus<3)then
        if(node(i)%tipus==2)then
    
          !if(node(i)%x<-0.2d0) gex(i,1)=1d0 !AP regionalization
    
          gex(i,3)=1d0  !receptors
          !gex(i,10)=1d0
          gex(i,14)=1d0
          !call random_number(a)
          !!!d=sqrt(node(i)%x**2+node(i)%y**2)
          !gex(i,2)=a!/(1+d)!*0.1d0
          !gex(i,9)=a!/(1+d)!*0.1d0
        end if
        !if(node(i)%tipus==1) gex(i,14)=1d0  !receptors
       
        gex(i,12)=1d0
        !gex(i,10)=1d0
      end if
      if(node(i)%tipus==3)then
        !if(node(i)%x<-0.2d0) gex(i,1)=1d0 !AP regionalization

        if(node(i)%z<=0)then
          gex(i,5)=1d0
          gex(i,13)=1d0
          gex(i,15)=1d0
        else
          gex(i,1)=1d0
          gex(i,16)=1d0
          gex(i,14)=1d0
        end if
          !call random_number(a)
          !!!d=sqrt(node(i)%x**2+node(i)%y**2)
          !gex(i,2)=a!/(1+d)!*0.1d0
          !gex(i,9)=a!/(1+d)!*0.1d0
      end if
      if(node(i)%tipus==4)then
        gex(i,16)=1d0
        !print*,"node border", node(i)%border
      end if
      !if(node(i)%hold>0)then
      !  gex(i,:)=0d0 !clear the hold nodes
      !  if(node(i)%tipus<3)then
      !    gex(i,12)=1d0
      !  end if
      !  if(node(i)%tipus==3)then
      !    gex(i,13)=1d0
      !  end if
      !end if

    end do


    call update_npag

    node(:)%talone=0.0d0

    !ramax=maxval(node(:)%da)*3

    
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1





subroutine epidermal_organ_reduc

!real*8::pericel,radio2,req1,req2,adjustreq,di
!integer::kk2

integer:: lx,ly
integer:: x1,x2,y1,y2,ii1,ii2,jj1,jj2
integer:: x3,x4,y3,y4,ii3,ii4,jj3,jj4,jj5,y5
integer, allocatable::cell_grid_epi(:,:),cell_grid_mes(:,:,:)
!integer:: periodic_borders,nborders
integer::geometry
integer::xlayer,ring,ring2

!******* #1 DEFINING SPATIAL DIMENSIONS *******

  geometry=1 !geometry=1 means hexagonal epithelium ; geometry=2 means rectangular epithelium
  
  if (geometry==1)then
	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=10    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=10   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=-0.5d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    do i=1,layer
      mesradicel(i)=mradicel
    end do
    
    xlayer=0 !this will make a layer of ECM beneath the layers of mesenchymal cells
    
  elseif(geometry==2)then
    layer=1      !number of planar cell layers
    lx=32
    ly=32
    
    allocate(cell_grid_epi(lx,ly),cell_grid_mes(lx,ly,layer))

  end if  
    !periodic_borders=1
    


    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(geometry==1)then
      if(radi>0.and.radicel>0)then
        j=0
        do i=1,radi-1
          j=j+i
        end do
        nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
        nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13
  
        j=0
        do i=1,radicel-1
          j=j+i
        end do
        ncelsepi=(6*j+1)
        ndepi=nodecel*ncelsepi
      else
        ncelsepi=0
        ndepi=0
        nodecel=0
        nodecela=0
      end if
  
      if(mradi>0.and.mradicel>0.and.layer>0)then
        if(packed==1)then !if mesenchyme is packed, we use mradi differently
          j=0
          do i=1,mradi-1
            j=j+i
          end do
          nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
          nodecela=2*nodecel+1
        else
          nodecel=mradi
          nodecela=2*nodecel+1
        end if
  
        ncelsmes=0
        do k=1,layer
          j=0 !;print*,"mesradicel",mesradicel(k)
          do i=1,mesradicel(k)-1
            j=j+i
          end do
          ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
        end do
        ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
        if(radi==0.and.radicel==0)then
          nodecel=mradi
          nodecela=2*nodecel+1
        else if(nodecel<mradi)then
          nodecel=radi
          nodecela=nodecel*2+1
        end if
      else
        ndmes=0
        ncelsmes=0
      end if
      !ndx=(6*j+1)*xlayer
      nd=ndepi+ndmes!+ndx
      ncels=ncelsepi+ncelsmes
      nda=nd+10
      ncals=ncels+10
    elseif(geometry==2)then

      ncelsepi=lx*ly
      ndepi=ncelsepi*2
      ncelsmes=ncelsepi*layer
      ndmes=ncelsmes
      nd=ndepi+ndmes
      ncels=ncelsepi+ncelsmes
      nodecel=2 ; nodecela=5
      nda=nd+10
      ncals=ncels+10
    end if
    

  !End initializing dimension parameters

  print*,"inicial nd",nd,"ncels",ncels,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=1
    !End Allocatations


    

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    deltamin=1d-2
    khold=1d0
    angletor=0.05
    k_bu=5d0
    ramax=0.35d0
    k_press=5d-1
    
    mi_xwall=-0.5d0
    
    ndmax=3d3
    
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=1
    ffu(2)=1 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !buoyancy

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(14)=1 !physical boundaries (walls)
    ffu(17)=1 !volume conservation (epithelial)
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise
    ffu(24)=0


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    !if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(mod(i,2)/=0)then  !basal
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
          !node(i)%reqp=-0.05d0
          !node(i)%req=0.20d0
        !node(i)%kplast=1d0
        else                      !apical
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !node(i)%kplast=1d1
        !  !node(i)%reqp=0.05d0
        !  !node(i)%req=0.30d0
        end if
        node(i)%req=0.25d0
        node(i)%reqcr=0.25d0 !; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.25
        node(i)%da=node(i)%req*1.60; 
        node(i)%ke=5d-1
        node(i)%tor=3d0
        node(i)%stor=1d0
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=3d0
        node(i)%kvol=3d1
      end do
    !end if

	!Mesenchyme
    if(layer>0)then
      do i=ndepi+1,ndepi+ndmes
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.60
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

	!ECM
    if(xlayer>0)then
      do i=ndepi+ndmes+1,nd
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.60
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if
    
    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)

    node(:)%hold=0


  if(geometry==1)then
    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    !if(radi>0.and.radicel>0)               call epiteli_sphere1(radi,radicel,radius,di,zepi)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if


    !ECM basal layer
    !do i=ndepi+ndmes+1,nd
    !  j=i-ndx
    !  node(i)%x=node(j)%x ; node(i)%y=node(j)%y ; node(i)%z=node(j)%z-0.5d0
    !  node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    !  node(i)%icel=-i ; node(i)%tipus=4
    !  !node(i)%border=1
    !end do

      
    j=0
    do i=1,radicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    do i=j+1,ndepi
      node(i)%hold=2 ;node(i)%repcel=1d0
    end do
    j=0
    do i=1,mradicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    k=0
    do i=1,mradicel-1
      k=k+i
    end do
    k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
      node(i)%hold=2;node(i)%repcel=1d0
      !node(i)%border=1
      !node(i)%da=node(i)%req*2.0
      !node(i)%orix=0 ; node(i)%oriy=0
      !node(i)%rep=1d1;node(i)%repcel=1d1
    end do
    !do i=ndepi+ndmes/2+j+1,nd
    !  node(i)%hold=1;node(i)%repcel=1d2
    !  !node(i)%border=1
    !end do
    do i=ndepi+k+1,nd
      node(i)%hold=2;node(i)%repcel=1d0
      node(i)%border=0
    end do

    
  elseif(geometry==2)then
	do i=1,ncels
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do


    !epithelial cells
    !origin of rectangle
    aa=node(1)%req*2
    bb=sqrt(aa**2-(aa*0.5)**2)
    a=-real(lx)*0.5*aa ; b=-real(ly)*0.5*bb
    !print*,"a",a,"b",b,"aa",aa,"bb",bb,"lx",lx,"ly",ly
    ii=0 ; jj=0
    do i=1,ly
      do j=1,lx
        ii=ii+1 ;jj=jj+1 ; cell_grid_epi(j,i)=jj
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii+1 ;node(ii)%tipus=2
        ii=ii+1
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii-1 ;node(ii)%tipus=1
        cels(jj)%node(1)=ii-1 ; cels(jj)%node(2)=ii ;cels(jj)%nunodes=2; cels(jj)%ctipus=1
        if(i==1.or.i==ly.or.j==1.or.j==lx)then
          node(ii)%hold=2;node(ii-1)%hold=2 ;
          node(ii)%repcel=1d1;node(ii-1)%repcel=1d1;
          !node(ii)%border=1;node(ii-1)%border=1;
        end if
      end do
    end do

    !mesenchymals
    do k=1,layer
      !print*,"k",k
      do i=1,ly
        do j=1,lx
          ii=ii+1 ; jj=jj+1 ; cell_grid_mes(j,i,k)=jj
          if(mod(i,2)==0)then
            node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          else
            node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          end if
          node(ii)%icel=jj ; node(ii)%altre=0 ;node(ii)%tipus=3
          cels(jj)%node(1)=ii ;cels(jj)%nunodes=1 ; cels(jj)%ctipus=3
          if(i==1 .or. i==ly .or. j==1 .or. j==lx)then
            node(ii)%hold=2 ;node(ii)%repcel=1d1 ; !node(ii)%border=2 ;
          end if
          !if(i==1 .or. i==ly .or. j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
          !if(j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
        end do
      end do
    end do
  end if

  !m_xwall=maxval(node(1:nd)%x)+0.01
  !mi_xwall=minval(node(1:nd)%x)-0.01
  !m_ywall=maxval(node(1:nd)%y)+0.01
  !mi_ywall=minval(node(1:nd)%y)-0.01
  


    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
  !node(:)%hold=2
  !node(:)%hold=0  

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    !if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    !end if
    !Mesenchymal
    if(layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=21
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=2.5d0 ; gen(1)%mu=5d-1 ;gen(1)%name="epithelial suprabasal marker"
      gen(1)%npre=1 ; allocate(gen(1)%pre(gen(1)%npre))
      gen(1)%pre(1)=7

      
      gen(2)%kindof=4 ; gen(2)%diffu=1d-2 ; gen(2)%mu=2.0d-2 ;gen(2)%name="placode expansion signal"
      gen(2)%npre=2 ; allocate(gen(2)%pre(gen(2)%npre))
      gen(2)%pre(1)=3 ; gen(2)%pre(2)=5

      
      gen(3)%kindof=1 ; gen(3)%diffu=1d-1 ; gen(3)%mu=5d-1 ;gen(3)%name="epithelial basal marker"
      !gen(3)%npost=1 ; allocate(gen(3)%post(gen(3)%npost))
      !gen(3)%post(1)=4

      gen(4)%kindof=1 ; gen(4)%diffu=0d0 ; gen(4)%mu=5d-1 ;gen(4)%name="epithelial knot marker"
      !gen(4)%kindof=8 ; gen(4)%diffu=0d0 ; gen(4)%mu=5d-1 !
      !gen(4)%npre=2 ; allocate(gen(4)%pre(gen(4)%npre))
      !gen(4)%pre(1)=2 ; gen(4)%pre(2)=3
      !gen(4)%npost=2 ; allocate(gen(4)%post(gen(4)%npost))
      !gen(4)%post(1)=2 ; gen(4)%post(2)=3

      gen(5)%kindof=1 ; gen(5)%diffu=0d0 ; gen(5)%mu=5d-1 ;gen(5)%name="mesenchymal marker"
      !gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost))
      !gen(5)%post(1)=6

      gen(6)%kindof=1 ; gen(6)%diffu=0d0 ; gen(6)%mu=5.0d-1 ;gen(6)%name="mesenchymal growth zone marker"
      !gen(6)%kindof=8 ; gen(6)%diffu=0d0 ; gen(6)%mu=5.0d-1 !
      !gen(6)%npre=2 ; allocate(gen(6)%pre(gen(6)%npre))
      !gen(6)%pre(1)=2 ; gen(6)%pre(2)=5
      !gen(6)%npost=2 ; allocate(gen(6)%post(gen(6)%npost))
      !gen(6)%post(1)=2 ; gen(6)%post(2)=5

      gen(7)%kindof=1 ; gen(7)%diffu=1d0 ; gen(7)%mu=2.5d-1 ;gen(7)%name="epithelial placode marker"

      gen(8)%kindof=1 ; gen(8)%diffu=0d0 ; gen(8)%mu=5d-1 ;gen(8)%name="mesenchymal placode marker"

      gen(9)%kindof=9 ; gen(9)%diffu=2.75d0 ; gen(9)%mu=8.0d-2 ;gen(9)%name="undef"
      gen(9)%npre=1 ; allocate(gen(9)%pre(gen(9)%npre))
      gen(9)%pre(1)=3 ;! gen(9)%pre(2)=5
      
      gen(10)%kindof=9 ; gen(10)%diffu=1d-1 ; gen(10)%mu=5.0d-1 ;gen(10)%name="activator signal"
      gen(10)%npre=1 ; allocate(gen(10)%pre(gen(10)%npre))
      gen(10)%pre(1)=7 ; !gen(10)%pre(2)=8

      gen(11)%kindof=9 ; gen(11)%diffu=2d0 ; gen(11)%mu=1.0d0 ;gen(11)%name="inhibitor signal"
      gen(11)%npre=1 ; allocate(gen(11)%pre(gen(11)%npre))
      gen(11)%pre(1)=7 ; !gen(11)%pre(2)=8
      !gen(11)%kindof=8 ; gen(11)%diffu=0d0 ; gen(11)%mu=5.0d-1 !
      !gen(11)%npre=2 ; allocate(gen(11)%pre(gen(11)%npre))
      !gen(11)%pre(1)=9 ; gen(11)%pre(2)=10
      !gen(11)%npost=2 ; allocate(gen(11)%post(gen(11)%npost))
      !gen(11)%post(1)=9 ; gen(11)%post(2)=10

      gen(12)%kindof=1 ; gen(12)%diffu=1d0 ; gen(12)%mu=2.5d-1 ;gen(12)%name="epithelial basal adhesion molecule"
      gen(13)%kindof=1 ; gen(13)%diffu=0d0 ; gen(13)%mu=5.0d-1 ;gen(13)%name="mesench adhesion molecule"

      gen(14)%kindof=1 ; gen(14)%diffu=0d0 ; gen(14)%mu=5.0d-1 ;gen(14)%name="epithelial interplacode marker"
      gen(15)%kindof=1 ; gen(15)%diffu=0d0 ; gen(15)%mu=5.0d-1  ;gen(15)%name="mesenchymal interplacode marker"
      gen(16)%kindof=1 ; gen(16)%diffu=0d0 ; gen(16)%mu=0.0d0  ;gen(16)%name="epithelial suprabasal adhesion molecule"
      
      gen(17)%kindof=9 ; gen(17)%diffu=1d-1 ; gen(17)%mu=5.0d-1  ;gen(17)%name="2ary activator-like, signal for mesenchyme"      
      gen(17)%npre=1 ; allocate(gen(17)%pre(gen(17)%npre))
      gen(17)%pre(1)=8 ; !gen(11)%pre(2)=8

      gen(18)%kindof=9 ; gen(18)%diffu=2d-1 ; gen(18)%mu=5.0d-1  ;gen(18)%name="undef"      
      gen(18)%npre=1 ; allocate(gen(18)%pre(gen(18)%npre))
      gen(18)%pre(1)=5 ; !gen(11)%pre(2)=8
      
      gen(19)%kindof=1 ; gen(19)%diffu=0d0 ; gen(19)%mu=5.0d-1  ;gen(19)%name="epithelial growth zone marker"      

      gen(20)%kindof=4 ; gen(20)%diffu=1d-2 ; gen(20)%mu=5.0d-2  ;gen(20)%name="undef"      
      gen(20)%npre=2 ; allocate(gen(20)%pre(gen(20)%npre))
      gen(20)%pre(1)=3 ; gen(20)%pre(2)=5

      gen(21)%kindof=1 ; gen(21)%diffu=0d0 ; gen(21)%mu=5.0d-1  ;gen(21)%name="undef"      

    !Gene-behavior interactions
    

    !!!!STANDARDIZED CELL BEHAVIORS
    
      gen(12)%wa(1)=1 !adhesion molecules
      gen(13)%wa(1)=2
      gen(16)%wa(1)=3
      
      !gen(2)%wa(nparam_per_node+8)=1d0     !1ary activator as a polarizing cue
      gen(7)%wa(nparam_per_node+16)=0d0   !epithelial placode marker triggers polar cell migration

      gen(10)%wa(nparam_per_node+8)=1d0     !1ary activator as a polarizing cue
      gen(6)%wa(nparam_per_node+11)=1d2     !1ary activator as a polarizing cue
      gen(19)%wa(nparam_per_node+11)=1d2     !1ary activator as a polarizing cue
      
      !regulation of cell proliferation

      !gen(7)%wa(nparam_per_node+2)=-5d-2  !placode and interplacode markers form an inhibitor threshold to cell proliferation
      !gen(8)%wa(nparam_per_node+2)=-5d-2
      !gen(14)%wa(nparam_per_node+2)=-1d1
      !gen(15)%wa(nparam_per_node+2)=-1d1
    
      !gen(4)%wa(nparam_per_node+2)=5d-2
      !gen(19)%wa(nparam_per_node+2)=2.0d-2
      !gen(6)%wa(nparam_per_node+2)=1.0d-1
      !gen(8)%wa(nparam_per_node+2)=1.0d-3


      gen(7)%wa(nparam_per_node+2)=1d-2  !placode and interplacode markers form an inhibitor threshold to cell proliferation
      !gen(8)%wa(nparam_per_node+2)=8d-2
      gen(6)%wa(nparam_per_node+2)=1d-1
      
      
    !!!!!!!!!!! 
    
      !gen(12)%wa(1)=1
      !gen(13)%wa(1)=2
      !!gen(7)%wa(1)=3
      !!gen(8)%wa(1)=4
      !
      !gen(16)%wa(1)=3
      !
      !gen(2)%wa(nparam_per_node+8)=1d0
      !gen(15)%wa(nparam_per_node+16)=0d0
      !!gen(7)%wa(13)=0.0d0
      !gen(8)%wa(nparam_per_node+2)=2d-2
      !gen(14)%wa(nparam_per_node+2)=8d-3
      

      
    !Gene-gene interactions

     !wavefront setting

    !!!!STANDARDIZED GENE NETWORK

      !gen(2)%w(2)=8d0  !Primary RD network module
      !gen(9)%w(2)=1d0
      !gen(2)%w(9)=-5d0

      !gen(18)%w(7)=1d0 !Primary activator activates a signal that only communicates to the mesenchyme
      !gen(8)%w(18)=5d0
      
      !gen(7)%w(2)=1d0  !placode markers activated by primary RD
      
      gen(7)%w(7)=1d1  !autoactivation of placode markers
      gen(8)%w(8)=1d1


      gen(7)%w(2)=1d1 !placode marker activates activator?
      gen(8)%w(2)=1d0

      gen(2)%w(7)=1d0  !placode markers activated by primary RD
      gen(2)%w(8)=1d0  !placode markers activated by primary RD
      
      gen(7)%w(3)=-3d0  !extended inhibition allows for threshold-based expression of placode markers
      gen(8)%w(3)=-1d5
      gen(7)%w(5)=-1d5
      gen(8)%w(5)=-2d0
      
      !gen(14)%w(7)=-1d5
      !gen(15)%w(8)=-1d5
      gen(7)%w(14)=-1d10
      gen(8)%w(15)=-1d10
      gen(15)%w(3)=-1d5
      gen(14)%w(5)=-1d5
      
      gen(15)%w(15)=1d1
      gen(14)%w(14)=1d1
      
      !gen(7)%w(14)=-1d2
      !gen(8)%w(15)=-1d2


      gen(3)%w(3)=1d0  !autoactivation of tissue markers
      gen(5)%w(5)=1d0
      
      gen(12)%w(12)=1d0 !autoactivation of adhesion molecules
      gen(13)%w(13)=1d0

      !gen(14)%w(3)=1d1 !autoactivation of interplacode markers
      !gen(15)%w(5)=1d1

      !gen(12)%w(7)=-2d0
      !gen(13)%w(8)=-2d0
      
      gen(10)%w(8)=1d-1 !placode marker activates secondary RD module

      gen(10)%w(10)=2d1 !Secondary RD network module
      gen(11)%w(10)=1d0
      gen(10)%w(11)=-5d0
      
      !gen(20)%w(14)=1d0
      !gen(10)%w(20)=-5d0

      !gen(20)%w(20)=1d0
      !gen(7)%w(20)=-1d1
      !gen(8)%w(20)=-1d1


      gen(10)%w(14)=-1d5 !secondary RD signals not allowed outside the placode
      gen(11)%w(14)=-1d5
      gen(10)%w(15)=-1d5
      gen(11)%w(15)=-1d5


      gen(17)%w(4)=1d0 !Secondary activator activates a signal that only communicates to the mesenchyme

      
      gen(4)%w(10)=4d0 !epithelial knot marker activated by 2ari activator
      gen(4)%w(7)=-1d1 !basal inhibitor threshold
      
      gen(19)%w(11)=1d0 !epithelial GZ marker activated by 2ari inhibitor
      gen(19)%w(4)=-1d5 !epithelial knot marker excludes GZ marker
      gen(19)%w(3)=-3d1 !epithelial GZ marker inhibition threshold
 
      !gen(4)%w(7)=-1d-1 !basal inhibitor threshold
      
      
      
      gen(6)%w(17)=3d1 !mesench GZ marker activated by 2ari activator
      

      gen(4)%w(7)=-1d1 !GZ markers threshold-based inhibition
      gen(6)%w(8)=-5d0 !
      
      gen(4)%w(8)=-1d5 !epi GZ marker not allowed on mesenchyme
      gen(6)%w(7)=-1d5 !mesench GZ marker not allowed on epithelium
      gen(19)%w(8)=-1d5
      
      
      gen(4)%w(14)=-1d5 !GZ markers not allowed outside the placode
      gen(4)%w(15)=-1d5 
      gen(6)%w(14)=-1d5 
      gen(6)%w(15)=-1d5 



      
      !gen(10)%w(10)=5d0 !Secondary RD network module
      !gen(1)%w(10)=1d0
      !gen(1)%w(1)=-1d1
      !gen(1)%w(8)=-1d5
      !gen(10)%w(7)=-1d5
      
      !gen(11)%w(10)=1d0
      !gen(11)%w(11)=-5d0
      !gen(11)%w(8)=-1d5


      !gen(10)%w(1)=-1d2
      !gen(11)%w(1)=-1d2
      
      
      
      !gen(1)%w(14)=1d-1 !external inhibitor of the 2ari RD module (something like EDA?)
      !gen(10)%w(1)=-1d1
      !gen(11)%w(1)=-1d1
      
      
      !gen(4)%w(1)=1d0  !growth zone markers (tipus knot)
      !gen(6)%w(10)=1d0

      !gen(4)%w(10)=-1d0
      !gen(6)%w(10)=-1d0

      !gen(6)%w(7)=-1d2
      !gen(4)%w(8)=-1d2
      !gen(6)%w(14)=-1d2
      !gen(4)%w(15)=-1d2      
      !gen(6)%w(15)=-1d2
      !gen(4)%w(14)=-1d2      

      


    !Adhesion molecules

	ntipusadh=3
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=5d-1
      kadh(1,2)=1d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=5d-1
      kadh(1,3)=5d-1 ; kadh(3,1)=kadh(1,3)
      kadh(3,3)=5d-1
      !kadh(3,4)=1d0 ; kadh(4,3)=kadh(3,4)
      !kadh(4,2)=5d-1 ; kadh(4,2)=kadh(2,4)
      !kadh(4,4)=3d0
      !kadh(5,5)=5d0
      !kadh(5,2)=5d-1 ; kadh(2,5)=kadh(5,2)
      !kadh(5,4)=5d-1 ; kadh(4,5)=kadh(5,4)

    end if

    !Gene expression on nodes

    
    !if(geometry==1)then
    !  gex(1,2)=1d0
    !  !gex(ndepi+1,2)=1d0
    !end if
    
    if(geometry==2)then
    !draw rectangles
    !!!!!!!!!!!!!!!!!!!!!!!!!
    jj1=1 ; jj2=ly
    ii1=1 ; ii2=lx
    do i=ii1,ii2
      do j=jj1,jj2
        k=cels(cell_grid_epi(i,j))%node(1)
        if(i<=2.or.i>=lx-1.or.j<=2.or.j>=ly-1)then
          gex(k,1)=1d0
          node(k)%hold=3
          node(node(k)%altre)%hold=3
          !gex(k,2)=0d0 ; gex(k,9)=0d0
        end if
        do ii=1,layer
          k=cels(cell_grid_mes(i,j,ii))%node(1)
          if(i<=2.or.i>=lx-1.or.j<=2.or.j>=ly-1)then
            gex(k,1)=1d0
            node(k)%hold=3
            !gex(k,2)=0d0 ; gex(k,9)=0d0
          end if
        end do
      end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!
    end if


    !determine the ring of expression

    ring=4
    ring2=4

    !j=0
    !do i=1,ring-1
    !  j=j+i
    !end do
    !ii=(6*j+1)
   
    do i=1,nd
    
      if(node(i)%tipus<3)then
        if(node(i)%tipus==2)then
    
          !if(node(i)%x<-0.2d0) gex(i,1)=1d0 !AP regionalization
    
          gex(i,3)=1d0  !receptors
          !gex(i,10)=1d0
          !gex(i,14)=1d0
          !call random_number(a)
          !!!d=sqrt(node(i)%x**2+node(i)%y**2)
          !gex(i,2)=a!/(1+d)!*0.1d0
          !gex(i,9)=a!/(1+d)!*0.1d0
        end if
        !if(node(i)%tipus==1) gex(i,14)=1d0  !receptors
       
        gex(i,12)=1d0
        !gex(i,10)=1d0

        
        !if(node(i)%icel<=ii) gex(i,7)=1d0
        
        d=sqrt(node(i)%x**2+node(i)%y**2)
        if(d-2*node(i)%req*(ring-1)<epsilod)then
          gex(i,7)=1d0
          !gex(i,14)=0d0
        end if
        
        if(d-2*node(i)%req*(ring2-1)>=epsilod)then
          gex(i,14)=1d0
          !gex(i,14)=0d0
        end if

        
      !if(node(i)%hold>0) gex(i,20)=1d0

      end if
      if(node(i)%tipus==3)then
        !if(node(i)%x<-0.2d0) gex(i,1)=1d0 !AP regionalization
        gex(i,5)=1d0
        gex(i,13)=1d0
        !gex(i,15)=1d0
          !call random_number(a)
          !!!d=sqrt(node(i)%x**2+node(i)%y**2)
          !gex(i,2)=a!/(1+d)!*0.1d0
          !gex(i,9)=a!/(1+d)!*0.1d0

      !if(node(i)%hold>0) gex(i,20)=1d0
 
        !if(node(i)%icel<=ncelsepi+ii) gex(i,8)=1d0
        
        d=sqrt(node(i)%x**2+node(i)%y**2)
        if(d-2*node(i)%req*(ring-1)<epsilod)then
          gex(i,8)=1d0
          !gex(i,15)=0d0
        end if
        
        if(d-2*node(i)%req*(ring2-1)>=epsilod)then
          !gex(i,8)=1d0
          gex(i,15)=1d0
        end if
        !print*,"d",d,"radius",2*node(i)%req*(ring-1)
      end if
      !if(node(i)%tipus==4)then
      !  gex(i,16)=1d0
      !  !print*,"node border", node(i)%border
      !end if
      !if(node(i)%hold>0)then
      !  gex(i,:)=0d0 !clear the hold nodes
      !  if(node(i)%tipus<3)then
      !    gex(i,12)=1d0
      !  end if
      !  if(node(i)%tipus==3)then
      !    gex(i,13)=1d0
      !  end if
      !end if

    end do


    call update_npag

    node(:)%talone=0.0d0

    !ramax=maxval(node(:)%da)*3

    
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





















subroutine epidermal_organ_2

!real*8::pericel,radio2,req1,req2,adjustreq,di
!integer::kk2

integer:: lx,ly
integer:: x1,x2,y1,y2,ii1,ii2,jj1,jj2
integer:: x3,x4,y3,y4,ii3,ii4,jj3,jj4,jj5,y5
integer, allocatable::cell_grid_epi(:,:),cell_grid_mes(:,:,:)
!integer:: periodic_borders,nborders
integer::geometry
integer::xlayer

!******* #1 DEFINING SPATIAL DIMENSIONS *******

  geometry=1 !geometry=1 means hexagonal epithelium ; geometry=2 means rectangular epithelium
  
  if (geometry==1)then
	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=10    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=10   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=-0.5d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    do i=1,layer
      mesradicel(i)=mradicel
    end do
    
    xlayer=1 !this will make a layer of ECM beneath the layers of mesenchymal cells
    
  elseif(geometry==2)then
    layer=1      !number of planar cell layers
    lx=32
    ly=32
    
    allocate(cell_grid_epi(lx,ly),cell_grid_mes(lx,ly,layer))

  end if  
    !periodic_borders=1
    


    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(geometry==1)then
      if(radi>0.and.radicel>0)then
        j=0
        do i=1,radi-1
          j=j+i
        end do
        nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
        nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13
  
        j=0
        do i=1,radicel-1
          j=j+i
        end do
        ncelsepi=(6*j+1)
        ndepi=nodecel*ncelsepi
      else
        ncelsepi=0
        ndepi=0
        nodecel=0
        nodecela=0
      end if
  
      if(mradi>0.and.mradicel>0.and.layer>0)then
        if(packed==1)then !if mesenchyme is packed, we use mradi differently
          j=0
          do i=1,mradi-1
            j=j+i
          end do
          nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
          nodecela=2*nodecel+1
        else
          nodecel=mradi
          nodecela=2*nodecel+1
        end if
  
        ncelsmes=0
        do k=1,layer
          j=0 !;print*,"mesradicel",mesradicel(k)
          do i=1,mesradicel(k)-1
            j=j+i
          end do
          ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
        end do
        ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
        if(radi==0.and.radicel==0)then
          nodecel=mradi
          nodecela=2*nodecel+1
        else if(nodecel<mradi)then
          nodecel=radi
          nodecela=nodecel*2+1
        end if
      else
        ndmes=0
        ncelsmes=0
      end if
      ndx=(6*j+1)*xlayer
      nd=ndepi+ndmes+ndx
      ncels=ncelsepi+ncelsmes
      nda=nd+10
      ncals=ncels+10
    elseif(geometry==2)then

      ncelsepi=lx*ly
      ndepi=ncelsepi*2
      ncelsmes=ncelsepi*layer
      ndmes=ncelsmes
      nd=ndepi+ndmes
      ncels=ncelsepi+ncelsmes
      nodecel=2 ; nodecela=5
      nda=nd+10
      ncals=ncels+10
    end if
    

  !End initializing dimension parameters

  print*,"inicial nd",nd,"ncels",ncels,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=1
    !End Allocatations


    

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    deltamin=1d-2
    khold=1d0
    angletor=0.05
    k_bu=5d0
    ramax=0.35d0
    
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=1
    ffu(2)=0 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !buoyancy

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(14)=0 !physical boundaries (walls)
    ffu(17)=1 !volume conservation (epithelial)
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise
    ffu(24)=0


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    !if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(mod(i,2)/=0)then  !basal
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
          !node(i)%reqp=-0.05d0
          !node(i)%req=0.20d0
        !node(i)%kplast=1d0
        else                      !apical
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !node(i)%kplast=1d1
        !  !node(i)%reqp=0.05d0
        !  !node(i)%req=0.30d0
        end if
        node(i)%req=0.25d0
        node(i)%reqcr=0.25d0 !; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.25
        node(i)%da=node(i)%req*1.60; 
        node(i)%ke=5d-1
        node(i)%tor=3d0
        node(i)%stor=1d0
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=1d0
        node(i)%kvol=1d0
      end do
    !end if

	!Mesenchyme
    if(layer>0)then
      do i=ndepi+1,ndepi+ndmes
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.60
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

	!ECM
    if(xlayer>0)then
      do i=ndepi+ndmes+1,nd
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.60
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if
    
    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)

    node(:)%hold=0


  if(geometry==1)then
    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    !if(radi>0.and.radicel>0)               call epiteli_sphere1(radi,radicel,radius,di,zepi)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if



    do i=ndepi+ndmes+1,nd
      j=i-ndx
      node(i)%x=node(j)%x ; node(i)%y=node(j)%y ; node(i)%z=node(j)%z-0.5d0
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
      node(i)%icel=-i ; node(i)%tipus=4
    end do

      
    j=0
    do i=1,radicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    do i=j+1,ndepi
      node(i)%hold=1 ;node(i)%repcel=5d0
    end do
    j=0
    do i=1,mradicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    k=0
    do i=1,mradicel-1
      k=k+i
    end do
    k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
      node(i)%hold=1;node(i)%repcel=5d0
      !node(i)%border=1
      !node(i)%da=node(i)%req*2.0
      !node(i)%orix=0 ; node(i)%oriy=0
      !node(i)%rep=1d1;node(i)%repcel=1d1
    end do
    !do i=ndepi+ndmes/2+j+1,nd
    !  node(i)%hold=1;node(i)%repcel=1d2
    !  !node(i)%border=1
    !end do
    do i=ndepi+k+1,nd
      node(i)%hold=1;node(i)%repcel=5d0
      !node(i)%border=1
    end do

    
  elseif(geometry==2)then
	do i=1,ncels
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do


    !epithelial cells
    !origin of rectangle
    aa=node(1)%req*2
    bb=sqrt(aa**2-(aa*0.5)**2)
    a=-real(lx)*0.5*aa ; b=-real(ly)*0.5*bb
    !print*,"a",a,"b",b,"aa",aa,"bb",bb,"lx",lx,"ly",ly
    ii=0 ; jj=0
    do i=1,ly
      do j=1,lx
        ii=ii+1 ;jj=jj+1 ; cell_grid_epi(j,i)=jj
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii+1 ;node(ii)%tipus=2
        ii=ii+1
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii-1 ;node(ii)%tipus=1
        cels(jj)%node(1)=ii-1 ; cels(jj)%node(2)=ii ;cels(jj)%nunodes=2; cels(jj)%ctipus=1
        if(i==1.or.i==ly.or.j==1.or.j==lx)then
          node(ii)%hold=2;node(ii-1)%hold=2 ;
          node(ii)%repcel=1d2;node(ii-1)%repcel=1d1;
          !node(ii)%border=1;node(ii-1)%border=1;
        end if
      end do
    end do

    !mesenchymals
    do k=1,layer
      !print*,"k",k
      do i=1,ly
        do j=1,lx
          ii=ii+1 ; jj=jj+1 ; cell_grid_mes(j,i,k)=jj
          if(mod(i,2)==0)then
            node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          else
            node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          end if
          node(ii)%icel=jj ; node(ii)%altre=0 ;node(ii)%tipus=3
          cels(jj)%node(1)=ii ;cels(jj)%nunodes=1 ; cels(jj)%ctipus=3
          if(i==1 .or. i==ly .or. j==1 .or. j==lx)then
            node(ii)%hold=2 ;node(ii)%repcel=1d1 ; !node(ii)%border=2 ;
          end if
          !if(i==1 .or. i==ly .or. j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
          !if(j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
        end do
      end do
    end do
  end if

  !m_xwall=maxval(node(1:nd)%x)+0.01
  !mi_xwall=minval(node(1:nd)%x)-0.01
  !m_ywall=maxval(node(1:nd)%y)+0.01
  !mi_ywall=minval(node(1:nd)%y)-0.01
  


    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
  !node(:)%hold=2
  !node(:)%hold=0  

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    !if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    !end if
    !Mesenchymal
    if(layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=16
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=4 ; gen(1)%diffu=1d-1 ; gen(1)%mu=5d-1 ;gen(1)%name="initial dermal signal"
      gen(1)%npre=1 ; allocate(gen(1)%pre(gen(1)%npre))
      gen(1)%pre(1)=3
      
      gen(2)%kindof=4 ; gen(2)%diffu=4d-2 ; gen(2)%mu=2.0d-2 ;gen(2)%name="primary activator signal"
      gen(2)%npre=1 ; allocate(gen(2)%pre(gen(2)%npre))
      gen(2)%pre(1)=3

      !gen(2)%npost=1 ; allocate(gen(2)%post(gen(2)%npost))
      !gen(2)%post(1)=4
      
      gen(3)%kindof=1 ; gen(3)%diffu=1d-1 ; gen(3)%mu=5d-1 ;gen(3)%name="epithelial marker"
      gen(3)%npost=1 ; allocate(gen(3)%post(gen(3)%npost))
      gen(3)%post(1)=1

      gen(4)%kindof=1 ; gen(4)%diffu=0d0 ; gen(4)%mu=5d-1 ;gen(4)%name="epithelial posterior placode marker"
      !gen(4)%kindof=8 ; gen(4)%diffu=0d0 ; gen(4)%mu=5d-1 !
      !gen(4)%npre=2 ; allocate(gen(4)%pre(gen(4)%npre))
      !gen(4)%pre(1)=2 ; gen(4)%pre(2)=3
      !gen(4)%npost=2 ; allocate(gen(4)%post(gen(4)%npost))
      !gen(4)%post(1)=2 ; gen(4)%post(2)=3

      gen(5)%kindof=1 ; gen(5)%diffu=0d0 ; gen(5)%mu=5d-1 ;gen(5)%name="mesenchymal marker"
      !gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost))
      !gen(5)%post(1)=6

      gen(6)%kindof=1 ; gen(6)%diffu=0d0 ; gen(6)%mu=5.0d-1 ;gen(6)%name="mesenchymal posterior placode marker"
      !gen(6)%kindof=8 ; gen(6)%diffu=0d0 ; gen(6)%mu=5.0d-1 !
      !gen(6)%npre=2 ; allocate(gen(6)%pre(gen(6)%npre))
      !gen(6)%pre(1)=2 ; gen(6)%pre(2)=5
      !gen(6)%npost=2 ; allocate(gen(6)%post(gen(6)%npost))
      !gen(6)%post(1)=2 ; gen(6)%post(2)=5

      gen(7)%kindof=9 ; gen(7)%diffu=1d0 ; gen(7)%mu=2.5d-1 ;gen(7)%name="epithelial placode marker"
      !gen(7)%npost=2 ; allocate(gen(7)%post(gen(7)%npost))
      !gen(7)%post(1)=2 ; gen(7)%post(2)=9

      gen(8)%kindof=9 ; gen(8)%diffu=0d0 ; gen(8)%mu=5d-1 ;gen(8)%name="mesenchymal placode marker"

      gen(9)%kindof=4 ; gen(9)%diffu=5d-1 ; gen(9)%mu=1.0d-1 ;gen(9)%name="primary inhibitor signal"
      gen(9)%npre=1 ; allocate(gen(9)%pre(gen(9)%npre))
      gen(9)%pre(1)=3
      !gen(9)%npost=1 ; allocate(gen(9)%post(gen(9)%npost))
      !gen(9)%post(1)=11
      
      gen(10)%kindof=4 ; gen(10)%diffu=1d-2 ; gen(10)%mu=5.0d-1 ;gen(10)%name="secondary activator signal"
      !gen(10)%npost=1 ; allocate(gen(10)%post(gen(10)%npost))
      !gen(10)%post(1)=11

      gen(11)%kindof=4 ; gen(11)%diffu=1d-1 ; gen(11)%mu=1.0d0 ;gen(11)%name="secondary inhibitor signal"
      !gen(11)%kindof=8 ; gen(11)%diffu=0d0 ; gen(11)%mu=5.0d-1 !
      !gen(11)%npre=2 ; allocate(gen(11)%pre(gen(11)%npre))
      !gen(11)%pre(1)=9 ; gen(11)%pre(2)=10
      !gen(11)%npost=2 ; allocate(gen(11)%post(gen(11)%npost))
      !gen(11)%post(1)=9 ; gen(11)%post(2)=10

      gen(12)%kindof=1 ; gen(12)%diffu=1d0 ; gen(12)%mu=2.5d-1 ;gen(12)%name="epithelial adhesion molecule"
      gen(13)%kindof=1 ; gen(13)%diffu=0d0 ; gen(13)%mu=5.0d-1 ;gen(13)%name="mesench adhesion molecule"

      gen(14)%kindof=1 ; gen(14)%diffu=0d0 ; gen(14)%mu=5.0d-1 ;gen(14)%name="epithelial interplacode marker"
      gen(15)%kindof=1 ; gen(15)%diffu=0d0 ; gen(15)%mu=5.0d-1  ;gen(15)%name="mesenchymal interplacode marker"
      gen(16)%kindof=1 ; gen(16)%diffu=0d0 ; gen(16)%mu=0.0d0  ;gen(16)%name="ECM adhesion molecule"
      
      
      
    !Gene-behavior interactions
    

      gen(12)%wa(1)=1
      gen(13)%wa(1)=2
      !gen(7)%wa(1)=3
      !gen(8)%wa(1)=4

      gen(16)%wa(1)=3
      
      gen(2)%wa(nparam_per_node+8)=1d0
      gen(15)%wa(nparam_per_node+16)=0d0
      gen(7)%wa(13)=0.0d0
      gen(15)%wa(nparam_per_node+2)=0d0

      

    !Gene-gene interactions

     !wavefront setting
      !gen(2)%w(1)=1d2  !this will spontaneously generate activator
      !gen(15)%w(1)=-1d1  !this will spontaneously generate activator
      !gen(1)%w(14)=1d0  !this will produce a ware
      !!gen(1)%w(15)=-1d0  !this is everywhere, inhibiting production of activator (in a threshold manner)
      !gen(14)%w(14)=1d0  !this is everywhere, inhibiting production of activator (in a threshold manner)
      !gen(15)%w(15)=1d0  !this is everywhere, inhibiting production of activator (in a threshold manner)
      !gen(2)%w(15)=-1d5  !this is everywhere, inhibiting production of activator (in a threshold manner)




      gen(1)%w(1)=1d0
      !gen(2)%w(1)=-1d5

      gen(2)%w(2)=8d0  !the activated receptor directly transcripts signal
      gen(9)%w(2)=1d0
      gen(2)%w(9)=-1d1

      !gen(7)%w(2)=1d0
      !gen(8)%w(2)=1d0
      gen(7)%w(3)=-3d0
      gen(8)%w(3)=-1d3
      gen(7)%w(5)=-1d3
      gen(8)%w(5)=-3d0
      gen(14)%w(7)=-1d2
      gen(15)%w(8)=-1d2
      
      
      gen(3)%w(3)=1d0  !autoactivation of receptors
      gen(5)%w(5)=1d0
      !gen(10)%w(10)=1d0
      
      gen(12)%w(12)=1d0 !autoactivation of adhesion molecules
      gen(13)%w(13)=1d0
      gen(14)%w(3)=1d0
      gen(15)%w(5)=1d0


      !gen(4)%w(1)=-1d3
      !gen(4)%w(7)=1d0
      
      !gen(6)%w(1)=-1d3
      !gen(6)%w(8)=1d0
      
      gen(12)%w(7)=-2d0
      gen(13)%w(8)=-2d0

      gen(1)%w(5)=1d0 !mesenchymal marker secretes initial inductive signal
      
      gen(7)%w(1)=1d-2
      
      gen(2)%w(7)=1d0
      gen(9)%w(7)=1d0
      
      gen(7)%w(2)=2d1
      gen(7)%w(9)=-5d0
      
      gen(7)%w(2)=8d0
      gen(8)%w(9)=-5d0
      
      !gen(5)%w(7)=-1d3
      !gen(3)%w(8)=-1d3
      


    !Adhesion molecules

	ntipusadh=3
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=5d-1
      kadh(1,2)=1d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=7d-1
      !kadh(1,3)=3d-1 ; kadh(3,1)=kadh(1,3)
      !kadh(3,3)=1d0
      !kadh(3,4)=1d0 ; kadh(4,3)=kadh(3,4)
      !kadh(4,2)=5d-1 ; kadh(4,2)=kadh(2,4)
      !kadh(4,4)=3d0
      kadh(3,3)=5d0
      kadh(3,2)=5d-1 ; kadh(2,3)=kadh(3,2)
      !kadh(5,4)=5d-1 ; kadh(4,5)=kadh(5,4)

    end if

    !Gene expression on nodes

    
    !if(geometry==1)then
    !  gex(1,2)=1d0
    !end if
    
    if(geometry==2)then
    !draw rectangles
    !!!!!!!!!!!!!!!!!!!!!!!!!
    jj1=1 ; jj2=ly
    ii1=1 ; ii2=lx
    do i=ii1,ii2
      do j=jj1,jj2
        k=cels(cell_grid_epi(i,j))%node(1)
        if(i<=2.or.i>=lx-1.or.j<=2.or.j>=ly-1)then
          gex(k,1)=1d0
          node(k)%hold=3
          node(node(k)%altre)%hold=3
          !gex(k,2)=0d0 ; gex(k,9)=0d0
        end if
        do ii=1,layer
          k=cels(cell_grid_mes(i,j,ii))%node(1)
          if(i<=2.or.i>=lx-1.or.j<=2.or.j>=ly-1)then
            gex(k,1)=1d0
            node(k)%hold=3
            !gex(k,2)=0d0 ; gex(k,9)=0d0
          end if
        end do
      end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!
    end if

   
    do i=1,nd
    
      if(node(i)%tipus<3)then
        if(node(i)%tipus==2)then
    
          !if(node(i)%x<-0.2d0) gex(i,1)=1d0 !AP regionalization
    
          gex(i,3)=1d0  !receptors
          gex(i,10)=1d0
          gex(i,14)=1d0
          !call random_number(a)
          !d=sqrt(node(i)%x**2+node(i)%y**2)
          !gex(i,2)=a/(1+d)!*0.1d0
          !gex(i,9)=a/(1+d)!*0.1d0
        end if
        !if(node(i)%tipus==1) gex(i,14)=1d0  !receptors
       
        gex(i,12)=1d0
        !gex(i,10)=1d0
      end if
      if(node(i)%tipus==3)then
        !if(node(i)%x<-0.2d0) gex(i,1)=1d0 !AP regionalization
        gex(i,5)=1d0
        gex(i,13)=1d0
        gex(i,15)=1d0

        call random_number(a)
        !d=sqrt(node(i)%x**2+node(i)%y**2)
        gex(i,5)=a!/(1+d)!*0.1d0
        gex(i,5)=a!/(1+d)!*0.1d0
        
      end if
      if(node(i)%tipus==4)then
        gex(i,16)=1d0
      end if
      !if(node(i)%hold>0)then
      !  gex(i,:)=0d0 !clear the hold nodes
      !  if(node(i)%tipus<3)then
      !    gex(i,12)=1d0
      !  end if
      !  if(node(i)%tipus==3)then
      !    gex(i,13)=1d0
      !  end if
      !end if

    end do


    call update_npag

    node(:)%talone=0.0d0

    !ramax=maxval(node(:)%da)*3

    
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine torque_test

!real*8::pericel,radio2,req1,req2,adjustreq,di
!integer::kk2

integer:: lx,ly
integer:: x1,x2,y1,y2,ii1,ii2,jj1,jj2
integer:: x3,x4,y3,y4,ii3,ii4,jj3,jj4,jj5,y5
integer, allocatable::cell_grid_epi(:,:),cell_grid_mes(:,:,:)
!integer:: periodic_borders,nborders
integer::geometry

!******* #1 DEFINING SPATIAL DIMENSIONS *******

  geometry=1 !geometry=1 means hexagonal epithelium ; geometry=2 means rectangular epithelium
  
  if (geometry==1)then
	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=2    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=2   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=-0.5d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    do i=1,layer
      mesradicel(i)=mradicel
    end do
    
  elseif(geometry==2)then
    layer=1      !number of planar cell layers
    lx=32
    ly=32
    
    allocate(cell_grid_epi(lx,ly),cell_grid_mes(lx,ly,layer))

  end if  
    !periodic_borders=1
    


    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(geometry==1)then
      if(radi>0.and.radicel>0)then
        j=0
        do i=1,radi-1
          j=j+i
        end do
        nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
        nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13
  
        j=0
        do i=1,radicel-1
          j=j+i
        end do
        ncelsepi=(6*j+1)
        ndepi=nodecel*ncelsepi
      else
        ncelsepi=0
        ndepi=0
        nodecel=0
        nodecela=0
      end if
  
      if(mradi>0.and.mradicel>0.and.layer>0)then
        if(packed==1)then !if mesenchyme is packed, we use mradi differently
          j=0
          do i=1,mradi-1
            j=j+i
          end do
          nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
          nodecela=2*nodecel+1
        else
          nodecel=mradi
          nodecela=2*nodecel+1
        end if
  
        ncelsmes=0
        do k=1,layer
          j=0 !;print*,"mesradicel",mesradicel(k)
          do i=1,mesradicel(k)-1
            j=j+i
          end do
          ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
        end do
        ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
        if(radi==0.and.radicel==0)then
          nodecel=mradi
          nodecela=2*nodecel+1
        else if(nodecel<mradi)then
          nodecel=radi
          nodecela=nodecel*2+1
        end if
      else
        ndmes=0
        ncelsmes=0
      end if
      nd=ndepi+ndmes
      ncels=ncelsepi+ncelsmes
      nda=nd+10
      ncals=ncels+10
    elseif(geometry==2)then

      ncelsepi=lx*ly
      ndepi=ncelsepi*2
      ncelsmes=ncelsepi*layer
      ndmes=ncelsmes
      nd=ndepi+ndmes
      ncels=ncelsepi+ncelsmes
      nodecel=2 ; nodecela=5
      nda=nd+10
      ncals=ncels+10
    end if
    

  !End initializing dimension parameters

  print*,"inicial nd",nd,"ncels",ncels,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=1
    !End Allocatations


    

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    deltamin=1d-3
    khold=1d1
    angletor=0.05
    k_bu=5d0
    ramax=0.35d0
    
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=1
    ffu(2)=0 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=2 !external signal source
    ffu(6)=0 !buoyancy

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(14)=0 !physical boundaries (walls)
    ffu(17)=1 !volume conservation (epithelial)
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise
    ffu(24)=0


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    !if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(mod(i,2)/=0)then  !basal
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d0
          !node(i)%reqp=-0.05d0
          !node(i)%req=0.20d0
        !node(i)%kplast=1d0
        else                      !apical
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d0
        !node(i)%kplast=1d1
        !  !node(i)%reqp=0.05d0
        !  !node(i)%req=0.30d0
        end if
        node(i)%req=0.25d0
        node(i)%reqcr=0.25d0 !; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.25
        node(i)%da=node(i)%req*1.60; 
        node(i)%ke=5d0
        node(i)%tor=3d1
        node(i)%stor=1d1
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=1d1
        node(i)%kvol=1d1
      end do
    !end if

	!Mesenchyme
    if(layer>0)then
      do i=ndepi+1,nd
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=5d0 ; node(i)%repcel=5d0
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.60
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)

    node(:)%hold=0


  if(geometry==1)then
    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    !if(radi>0.and.radicel>0)               call epiteli_sphere1(radi,radicel,radius,di,zepi)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if

    !j=0
    !do i=1,radicel-2
    !  j=j+i
    !end do
    !j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    !do i=j+1,ndepi
    !  node(i)%hold=1 ;node(i)%repcel=1d2
    !end do
    !j=0
    !do i=1,mradicel-2
    !  j=j+i
    !end do
    !j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    !k=0
    !do i=1,mradicel-1
    !  k=k+i
    !end do
    !k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    !do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
    !  node(i)%hold=1;node(i)%repcel=1d2
    !  !node(i)%border=1
    !  !node(i)%da=node(i)%req*2.0
    !  !node(i)%orix=0 ; node(i)%oriy=0
    !  !node(i)%rep=1d1;node(i)%repcel=1d1
    !end do
    !do i=ndepi+k+1,nd
    !  node(i)%hold=1;node(i)%repcel=1d2
    !  !node(i)%border=1
    !end do

  elseif(geometry==2)then
	do i=1,ncels
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do


    !epithelial cells
    !origin of rectangle
    aa=node(1)%req*2
    bb=sqrt(aa**2-(aa*0.5)**2)
    a=-real(lx)*0.5*aa ; b=-real(ly)*0.5*bb
    !print*,"a",a,"b",b,"aa",aa,"bb",bb,"lx",lx,"ly",ly
    ii=0 ; jj=0
    do i=1,ly
      do j=1,lx
        ii=ii+1 ;jj=jj+1 ; cell_grid_epi(j,i)=jj
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii+1 ;node(ii)%tipus=2
        ii=ii+1
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii-1 ;node(ii)%tipus=1
        cels(jj)%node(1)=ii-1 ; cels(jj)%node(2)=ii ;cels(jj)%nunodes=2; cels(jj)%ctipus=1
        if(i==1.or.i==ly.or.j==1.or.j==lx)then
          node(ii)%hold=2;node(ii-1)%hold=2 ;
          node(ii)%repcel=1d2;node(ii-1)%repcel=1d2;
          !node(ii)%border=1;node(ii-1)%border=1;
        end if
      end do
    end do

    !mesenchymals
    do k=1,layer
      !print*,"k",k
      do i=1,ly
        do j=1,lx
          ii=ii+1 ; jj=jj+1 ; cell_grid_mes(j,i,k)=jj
          if(mod(i,2)==0)then
            node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          else
            node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          end if
          node(ii)%icel=jj ; node(ii)%altre=0 ;node(ii)%tipus=3
          cels(jj)%node(1)=ii ;cels(jj)%nunodes=1 ; cels(jj)%ctipus=3
          if(i==1 .or. i==ly .or. j==1 .or. j==lx)then
            node(ii)%hold=2 ;node(ii)%repcel=1d2 ; !node(ii)%border=2 ;
          end if
          !if(i==1 .or. i==ly .or. j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
          !if(j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
        end do
      end do
    end do
  end if

  !m_xwall=maxval(node(1:nd)%x)+0.01
  !mi_xwall=minval(node(1:nd)%x)-0.01
  !m_ywall=maxval(node(1:nd)%y)+0.01
  !mi_ywall=minval(node(1:nd)%y)-0.01
  


    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
  !node(:)%hold=2
  !node(:)%hold=0  

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    !if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    !end if
    !Mesenchymal
    if(layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=1
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=0d0 ; gen(1)%mu=0d0 !
      
      
      
      
    !Gene-behavior interactions
    

      gen(1)%wa(nparam_per_node+16)=1d0




    !Adhesion molecules

	ntipusadh=0
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=5d0
      kadh(1,2)=3d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=5d0
      !kadh(3,3)=5d0
      !kadh(3,4)=2d0 ; kadh(4,3)=kadh(3,4)
      !kadh(4,4)=5d0
    end if

    !Gene expression on nodes

    
    if(geometry==1)then
!      gex(1,2)=1d0
      gex(15:21,1)=1d0
      !gex(3,1)=1d0
    end if
    
    !if(geometry==2)then
    !!draw rectangles
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !jj1=1 ; jj2=ly
    !ii1=1 ; ii2=lx
    !do i=ii1,ii2
    !  do j=jj1,jj2
    !    k=cels(cell_grid_epi(i,j))%node(1)
    !    if(i<=2.or.i>=lx-1.or.j<=2.or.j>=ly-1)then
    !      gex(k,1)=1d0
    !      node(k)%hold=3
    !      node(node(k)%altre)%hold=3
    !      !gex(k,2)=0d0 ; gex(k,9)=0d0
    !    end if
    !    do ii=1,layer
    !      k=cels(cell_grid_mes(i,j,ii))%node(1)
    !      if(i<=2.or.i>=lx-1.or.j<=2.or.j>=ly-1)then
    !        gex(k,1)=1d0
    !        node(k)%hold=3
    !        !gex(k,2)=0d0 ; gex(k,9)=0d0
    !      end if
    !    end do
    !  end do
    !end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !end if
    !
    !
    !do i=1,nd
    !
    !  if(node(i)%tipus<3)then
    !    if(node(i)%tipus==2)then
    !
    !      if(node(i)%x<-0.2d0) gex(i,1)=1d0 !AP regionalization
    !
    !      gex(i,3)=1d0  !receptors
    !      gex(i,10)=1d0
    !      gex(i,14)=1d0
    !      !call random_number(a)
    !      !gex(i,2)=a!*0.1d0
    !      !gex(i,9)=a!*0.1d0
    !    end if
    !    !if(node(i)%tipus==1) gex(i,14)=1d0  !receptors
    !   
    !    gex(i,12)=1d0
    !    !gex(i,10)=1d0
    !  end if
    !  if(node(i)%tipus==3)then
    !    if(node(i)%x<-0.2d0) gex(i,1)=1d0 !AP regionalization
    !    gex(i,5)=1d0
    !    gex(i,13)=1d0
    !    gex(i,15)=1d0
    !  end if
    !  !if(node(i)%hold>0)then
    !  !  gex(i,:)=0d0 !clear the hold nodes
    !  !  if(node(i)%tipus<3)then
    !  !    gex(i,12)=1d0
    !  !  end if
    !  !  if(node(i)%tipus==3)then
    !  !    gex(i,13)=1d0
    !  !  end if
    !  !end if
    !
    !end do


    call update_npag

    node(:)%talone=0.0d0

    !ramax=maxval(node(:)%da)*3

    
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine epidermal_rectangle_test

!real*8::pericel,radio2,req1,req2,adjustreq,di
!integer::kk2

integer:: lx,ly
integer:: x1,x2,y1,y2,ii1,ii2,jj1,jj2
integer:: x3,x4,y3,y4,ii3,ii4,jj3,jj4,jj5,y5
integer, allocatable::cell_grid_epi(:,:),cell_grid_mes(:,:,:)
!integer:: periodic_borders,nborders


!******* #1 DEFINING SPATIAL DIMENSIONS *******

	!epithelium's dimension parameters
    radi=0       !number of radial layers of nodes per cell
    radicel=0    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=0     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=0      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=0   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=-0.5d0  !z-position of uppermost layer
    
    lx=10
    ly=10
    
    allocate(cell_grid_epi(lx,ly),cell_grid_mes(lx,ly,layer))

    
    !periodic_borders=1
    


    !ECM dimension parameters?

!************************************************

    !if(mod(ly,2)==0)then;
    !  ncelsepi=(ly/2)*(lx+lx-1)
    !else
    !  ncelsepi=((ly-1)/2)*(lx+lx-1)+lx-1
    !endif
    ncelsepi=lx*ly

    !ncelsepi=(lx-1)*(ly-1)
    !if(mod(ly,2)==0)then; ncelsepi=ncelsepi+lx-2 ;else;ncelsepi=ncelsepi+lx; endif
    ndepi=ncelsepi*2
    ncelsmes=ncelsepi*layer
    ndmes=ncelsmes
    nd=ndepi+ndmes
    ncels=ncelsepi+ncelsmes
    nodecel=2 ; nodecela=5
    

    nda=nd+10
	ncals=ncels+10
    

  !End initializing dimension parameters

  print*,"inicial nd",nd,"ncels",ncels,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=1
    !End Allocatations


    

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.1d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1.0d0
    deltamin=1d-3
    khold=2d0
    angletor=0.05
    k_bu=1d-1
    ramax=0.35d0
    
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=1
    ffu(2)=0 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !buoyancy

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(14)=0 !physical boundaries (walls)
    ffu(17)=1 !volume conservation (epithelial)
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise
    ffu(24)=0


    if(ffu(24)==1)then
      allocate(borders_neigh((2*lx+2*ly-4)*2+(2*lx+2*ly-4)*layer,4))
      borders_neigh=0
      nborders=0
    end if


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    !if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(mod(i,2)/=0)then  !basal
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d0
          !node(i)%reqp=-0.05d0
          !node(i)%req=0.20d0
        !node(i)%kplast=1d0
        else                      !apical
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d0
        !node(i)%kplast=1d1
        !  !node(i)%reqp=0.05d0
        !  !node(i)%req=0.30d0
        end if
        node(i)%req=0.25d0
        node(i)%reqcr=0.25d0 !; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.25
        node(i)%da=node(i)%req*1.60; 
        node(i)%ke=5d0
        node(i)%tor=3d1
        node(i)%stor=1d1
        node(i)%mo=temp
        node(i)%dmo=desmax
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=1d1
        node(i)%kvol=1d1
      end do
    !end if

	!Mesenchyme
    if(layer>0)then
      do i=ndepi+1,nd
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=5d0 ; node(i)%repcel=5d0
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*1.60
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=desmax
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)

    node(:)%hold=0

	do i=1,ncels
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do

    !epithelial cells
    !origin of rectangle
    aa=node(1)%req*2
    bb=sqrt(aa**2-(aa*0.5)**2)
    a=-real(lx)*0.5*aa ; b=-real(ly)*0.5*bb
    !print*,"a",a,"b",b,"aa",aa,"bb",bb,"lx",lx,"ly",ly
    ii=0 ; jj=0
    do i=1,ly
      do j=1,lx
        ii=ii+1 ;jj=jj+1 ; cell_grid_epi(j,i)=jj
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii+1 ;node(ii)%tipus=2
        ii=ii+1
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii-1 ;node(ii)%tipus=1
        cels(jj)%node(1)=ii-1 ; cels(jj)%node(2)=ii ;cels(jj)%nunodes=2; cels(jj)%ctipus=1
        if(i==1.or.i==ly.or.j==1.or.j==lx)then
          node(ii)%hold=2;node(ii-1)%hold=2 ;
          node(ii)%repcel=1d2;node(ii-1)%repcel=1d2;
          !node(ii)%border=1;node(ii-1)%border=1;
          if(ffu(24)==1)then
            nborders=nborders+1
            node(ii-1)%border=nborders
            nborders=nborders+1
            node(ii)%border=nborders
          end if
        end if
      end do
    end do

    !mesenchymals
    do k=1,layer
      !print*,"k",k
      do i=1,ly
        do j=1,lx
          ii=ii+1 ; jj=jj+1 ; cell_grid_mes(j,i,k)=jj
          if(mod(i,2)==0)then
            node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          else
            node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          end if
          node(ii)%icel=jj ; node(ii)%altre=0 ;node(ii)%tipus=3
          cels(jj)%node(1)=ii ;cels(jj)%nunodes=1 ; cels(jj)%ctipus=3
          if(i==1 .or. i==ly .or. j==1 .or. j==lx)then
            node(ii)%hold=2 ;node(ii)%repcel=1d2 ; !node(ii)%border=2 ;
            if(ffu(24)==1)then
              nborders=nborders+1
              node(ii)%border=nborders
            end if
          end if
          !if(i==1 .or. i==ly .or. j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
          !if(j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
        end do
      end do
    end do


  m_xwall=maxval(node(1:nd)%x)+0.01
  mi_xwall=minval(node(1:nd)%x)-0.01
  m_ywall=maxval(node(1:nd)%y)+0.01
  mi_ywall=minval(node(1:nd)%y)-0.01
  
    if(ffu(24)==1)then !set the neighbors of the periodic boundaries
      j=1
      do i=2,lx-1 !i=1,j=1 comença ficat cap endins
        k=cell_grid_epi(i,j) !the cell
        kk=cell_grid_epi(i,ly)!first neighbor
        borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
        borders_neigh(node(cels(k)%node(2))%border,1)=cels(kk)%node(2)

        kk=cell_grid_epi(i+1,ly)!second neighbor
        borders_neigh(node(cels(k)%node(1))%border,2)=cels(kk)%node(1)
        borders_neigh(node(cels(k)%node(2))%border,2)=cels(kk)%node(2)

        do ii=1,layer
          k=cell_grid_mes(i,j,ii) !the cell
          kk=cell_grid_mes(i,ly,ii)!first neighbor
          borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
  
          kk=cell_grid_mes(i+1,ly,ii)!second neighbor
          borders_neigh(node(cels(k)%node(1))%border,2)=cels(kk)%node(1)
        end do
      end do
      j=ly
      do i=2,lx-1 !i=1,j=ly comença ficat cap enfora
        k=cell_grid_epi(i,j) !the cell
        kk=cell_grid_epi(i-1,1)!first neighbor
        borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
        borders_neigh(node(cels(k)%node(2))%border,1)=cels(kk)%node(2)

        kk=cell_grid_epi(i,1)!second neighbor
        borders_neigh(node(cels(k)%node(1))%border,2)=cels(kk)%node(1)
        borders_neigh(node(cels(k)%node(2))%border,2)=cels(kk)%node(2)

        do ii=1,layer
          k=cell_grid_mes(i,j,ii) !the cell
          kk=cell_grid_mes(i-1,1,ii)!first neighbor
          borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
  
          kk=cell_grid_mes(i,1,ii)!second neighbor
          borders_neigh(node(cels(k)%node(1))%border,2)=cels(kk)%node(1)
        end do
      end do
           
      i=1
      do j=2,ly-1 !i=1,j=1 comença ficat cap endins
        if(mod(j,2)==0)then
          k=cell_grid_epi(i,j) !the cell
          kk=cell_grid_epi(lx,j-1)!first neighbor
          borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
          borders_neigh(node(cels(k)%node(2))%border,1)=cels(kk)%node(2)
  
          kk=cell_grid_epi(lx,j)!second neighbor
          borders_neigh(node(cels(k)%node(1))%border,2)=cels(kk)%node(1)
          borders_neigh(node(cels(k)%node(2))%border,2)=cels(kk)%node(2)

          kk=cell_grid_epi(lx,j+1)!third neighbor
          borders_neigh(node(cels(k)%node(1))%border,3)=cels(kk)%node(1)
          borders_neigh(node(cels(k)%node(2))%border,3)=cels(kk)%node(2)
          do ii=1,layer
            k=cell_grid_mes(i,j,ii) !the cell
            kk=cell_grid_mes(lx,j-1,ii)!first neighbor
            borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
  
            kk=cell_grid_mes(lx,j,ii)!second neighbor
            borders_neigh(node(cels(k)%node(1))%border,2)=cels(kk)%node(1)

            kk=cell_grid_mes(lx,j+1,ii)!third neighbor
            borders_neigh(node(cels(k)%node(1))%border,3)=cels(kk)%node(1)
          end do
        else
          k=cell_grid_epi(i,j) !the cell
          kk=cell_grid_epi(lx,j)!first neighbor
          borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
          borders_neigh(node(cels(k)%node(2))%border,1)=cels(kk)%node(2)
          do ii=1,layer
            k=cell_grid_mes(i,j,ii) !the cell
            kk=cell_grid_mes(lx,j,ii)!first neighbor
            borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
          end do
        end if
      end do
  
      i=lx
      do j=2,ly-1 !i=lx,j=1 comença ficat cap endins
        if(mod(j,2)/=0)then
          k=cell_grid_epi(i,j) !the cell
          kk=cell_grid_epi(1,j-1)!first neighbor
          borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
          borders_neigh(node(cels(k)%node(2))%border,1)=cels(kk)%node(2)
  
          kk=cell_grid_epi(1,j)!second neighbor
          borders_neigh(node(cels(k)%node(1))%border,2)=cels(kk)%node(1)
          borders_neigh(node(cels(k)%node(2))%border,2)=cels(kk)%node(2)

          kk=cell_grid_epi(1,j+1)!third neighbor
          borders_neigh(node(cels(k)%node(1))%border,3)=cels(kk)%node(1)
          borders_neigh(node(cels(k)%node(2))%border,3)=cels(kk)%node(2)
          do ii=1,layer
            k=cell_grid_mes(i,j,ii) !the cell
            kk=cell_grid_mes(1,j-1,ii)!first neighbor
            borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
  
            kk=cell_grid_mes(1,j,ii)!second neighbor
            borders_neigh(node(cels(k)%node(1))%border,2)=cels(kk)%node(1)

            kk=cell_grid_mes(1,j+1,ii)!third neighbor
            borders_neigh(node(cels(k)%node(1))%border,3)=cels(kk)%node(1)
          end do
        else
          k=cell_grid_epi(i,j) !the cell
          kk=cell_grid_epi(1,j)!first neighbor
          borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
          borders_neigh(node(cels(k)%node(2))%border,1)=cels(kk)%node(2)
          do ii=1,layer
            k=cell_grid_mes(i,j,ii) !the cell
            kk=cell_grid_mes(1,j,ii)!first neighbor
            borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
          end do
        end if
      end do
      !the corners
      i=1 ; j=1
        k=cell_grid_epi(i,j) !the cell
        kk=cell_grid_epi(lx,1)!first neighbor
        borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
        borders_neigh(node(cels(k)%node(2))%border,1)=cels(kk)%node(2)
  
        kk=cell_grid_epi(1,ly)!second neighbor
        borders_neigh(node(cels(k)%node(1))%border,2)=cels(kk)%node(1)
        borders_neigh(node(cels(k)%node(2))%border,2)=cels(kk)%node(2)

        kk=cell_grid_epi(2,ly)!third neighbor
        borders_neigh(node(cels(k)%node(1))%border,3)=cels(kk)%node(1)
        borders_neigh(node(cels(k)%node(2))%border,3)=cels(kk)%node(2)
        do ii=1,layer
          k=cell_grid_mes(i,j,ii) !the cell
          kk=cell_grid_mes(lx,1,ii)!first neighbor
          borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
  
          kk=cell_grid_mes(1,ly,ii)!second neighbor
          borders_neigh(node(cels(k)%node(1))%border,2)=cels(kk)%node(1)

          kk=cell_grid_mes(2,ly,ii)!third neighbor
          borders_neigh(node(cels(k)%node(1))%border,3)=cels(kk)%node(1)
        end do
      i=lx ; j=1
        k=cell_grid_epi(i,j) !the cell
        kk=cell_grid_epi(1,2)!first neighbor
        borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
        borders_neigh(node(cels(k)%node(2))%border,1)=cels(kk)%node(2)
  
        kk=cell_grid_epi(1,1)!second neighbor
        borders_neigh(node(cels(k)%node(1))%border,2)=cels(kk)%node(1)
        borders_neigh(node(cels(k)%node(2))%border,2)=cels(kk)%node(2)

        kk=cell_grid_epi(1,ly)!third neighbor
        borders_neigh(node(cels(k)%node(1))%border,3)=cels(kk)%node(1)
        borders_neigh(node(cels(k)%node(2))%border,3)=cels(kk)%node(2)

        kk=cell_grid_epi(lx,ly)!fourth neighbor
        borders_neigh(node(cels(k)%node(1))%border,4)=cels(kk)%node(1)
        borders_neigh(node(cels(k)%node(2))%border,4)=cels(kk)%node(2)        
        do ii=1,layer
          k=cell_grid_mes(i,j,ii) !the cell
          kk=cell_grid_mes(1,2,ii)!first neighbor
          borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
  
          kk=cell_grid_mes(1,1,ii)!second neighbor
          borders_neigh(node(cels(k)%node(1))%border,2)=cels(kk)%node(1)

          kk=cell_grid_mes(1,ly,ii)!third neighbor
          borders_neigh(node(cels(k)%node(1))%border,3)=cels(kk)%node(1)

          kk=cell_grid_mes(lx,ly,ii)!third neighbor
          borders_neigh(node(cels(k)%node(1))%border,4)=cels(kk)%node(1)
        end do
        i=1 ; j=ly
        k=cell_grid_epi(i,j) !the cell
        kk=cell_grid_epi(1,1)!first neighbor
        borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
        borders_neigh(node(cels(k)%node(2))%border,1)=cels(kk)%node(2)
  
        kk=cell_grid_epi(lx,1)!second neighbor
        borders_neigh(node(cels(k)%node(1))%border,2)=cels(kk)%node(1)
        borders_neigh(node(cels(k)%node(2))%border,2)=cels(kk)%node(2)

        kk=cell_grid_epi(lx,ly)!third neighbor
        borders_neigh(node(cels(k)%node(1))%border,3)=cels(kk)%node(1)
        borders_neigh(node(cels(k)%node(2))%border,3)=cels(kk)%node(2)

        kk=cell_grid_epi(lx,ly-1)!fourth neighbor
        borders_neigh(node(cels(k)%node(1))%border,4)=cels(kk)%node(1)
        borders_neigh(node(cels(k)%node(2))%border,4)=cels(kk)%node(2)        
        do ii=1,layer
          k=cell_grid_mes(i,j,ii) !the cell
          kk=cell_grid_mes(1,1,ii)!first neighbor
          borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
  
          kk=cell_grid_mes(lx,1,ii)!second neighbor
          borders_neigh(node(cels(k)%node(1))%border,2)=cels(kk)%node(1)

          kk=cell_grid_mes(lx,ly,ii)!third neighbor
          borders_neigh(node(cels(k)%node(1))%border,3)=cels(kk)%node(1)

          kk=cell_grid_mes(lx,ly-1,ii)!third neighbor
          borders_neigh(node(cels(k)%node(1))%border,4)=cels(kk)%node(1)
        end do
  
      i=lx ; j=ly
        k=cell_grid_epi(i,j) !the cell
        kk=cell_grid_epi(lx-1,1)!first neighbor
        borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
        borders_neigh(node(cels(k)%node(2))%border,1)=cels(kk)%node(2)
  
        kk=cell_grid_epi(lx,1)!second neighbor
        borders_neigh(node(cels(k)%node(1))%border,2)=cels(kk)%node(1)
        borders_neigh(node(cels(k)%node(2))%border,2)=cels(kk)%node(2)

        kk=cell_grid_epi(1,ly)!third neighbor
        borders_neigh(node(cels(k)%node(1))%border,3)=cels(kk)%node(1)
        borders_neigh(node(cels(k)%node(2))%border,3)=cels(kk)%node(2)
        do ii=1,layer
          k=cell_grid_mes(i,j,ii) !the cell
          kk=cell_grid_mes(lx-1,1,ii)!first neighbor
          borders_neigh(node(cels(k)%node(1))%border,1)=cels(kk)%node(1)
  
          kk=cell_grid_mes(lx,1,ii)!second neighbor
          borders_neigh(node(cels(k)%node(1))%border,2)=cels(kk)%node(1)

          kk=cell_grid_mes(1,ly,ii)!third neighbor
          borders_neigh(node(cels(k)%node(1))%border,3)=cels(kk)%node(1)
        end do
    end if

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
  !node(:)%hold=2
  
    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    !if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    !end if
    !Mesenchymal
    if(layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=15
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=0d0 ; gen(1)%mu=0d0 !
      
      gen(2)%kindof=4 ; gen(2)%diffu=5d-2 ; gen(2)%mu=1.0d-1 !activator signal
      gen(2)%npost=1 ; allocate(gen(2)%post(gen(2)%npost))
      gen(2)%post(1)=4
      
      gen(3)%kindof=2 ; gen(3)%diffu=1d-1 ; gen(3)%mu=5d-1 !activator receptor inactive epithelial
      gen(3)%npost=1 ; allocate(gen(3)%post(gen(3)%npost))
      gen(3)%post(1)=4

      gen(4)%kindof=8 ; gen(4)%diffu=0d0 ; gen(4)%mu=5d-1 !activator receptor active epithelial
      gen(4)%npre=2 ; allocate(gen(4)%pre(gen(4)%npre))
      gen(4)%pre(1)=2 ; gen(4)%pre(2)=3
      gen(4)%npost=2 ; allocate(gen(4)%post(gen(4)%npost))
      gen(4)%post(1)=2 ; gen(4)%post(2)=3

      gen(5)%kindof=2 ; gen(5)%diffu=0d0 ; gen(5)%mu=5d-1 !activator receptor inactive mesench
      gen(5)%npost=1 ; allocate(gen(5)%post(gen(5)%npost))
      gen(5)%post(1)=6

      gen(6)%kindof=8 ; gen(6)%diffu=0d0 ; gen(6)%mu=5.0d-1 !activator receptor active mesench
      gen(6)%npre=2 ; allocate(gen(6)%pre(gen(6)%npre))
      gen(6)%pre(1)=2 ; gen(6)%pre(2)=5
      gen(6)%npost=2 ; allocate(gen(6)%post(gen(6)%npost))
      gen(6)%post(1)=2 ; gen(6)%post(2)=5

      gen(7)%kindof=1 ; gen(7)%diffu=1d0 ; gen(7)%mu=2.5d-1 !

      gen(8)%kindof=1 ; gen(8)%diffu=0d0 ; gen(8)%mu=5d-1 !

      gen(9)%kindof=4 ; gen(9)%diffu=2d0 ; gen(9)%mu=5.0d-1 !inhibitor signal, epithelial
      gen(9)%npost=1 ; allocate(gen(9)%post(gen(9)%npost))
      gen(9)%post(1)=11
      
      gen(10)%kindof=2 ; gen(10)%diffu=0d0 ; gen(10)%mu=1.0d0 !inhibitor receptor inactive epithelial
      gen(10)%npost=1 ; allocate(gen(10)%post(gen(10)%npost))
      gen(10)%post(1)=11

      gen(11)%kindof=8 ; gen(11)%diffu=0d0 ; gen(11)%mu=5.0d-1 !inhibitor receptor active epithelial
      gen(11)%npre=2 ; allocate(gen(11)%pre(gen(11)%npre))
      gen(11)%pre(1)=9 ; gen(11)%pre(2)=10
      gen(11)%npost=2 ; allocate(gen(11)%post(gen(11)%npost))
      gen(11)%post(1)=9 ; gen(11)%post(2)=10

      gen(12)%kindof=1 ; gen(12)%diffu=1d0 ; gen(12)%mu=2.5d-1 !epithelial adhesion molecule
      gen(13)%kindof=1 ; gen(13)%diffu=0d0 ; gen(13)%mu=5.0d-1 !mesench adhesion molecule

      gen(14)%kindof=1 ; gen(14)%diffu=0d0 ; gen(14)%mu=5.0d-1 !
      gen(15)%kindof=1 ; gen(15)%diffu=0d0 ; gen(15)%mu=5.0d-1  !
      
      
      
    !Gene-behavior interactions
    

      !gen(7)%wa(1)=1
      !gen(8)%wa(1)=2
      gen(12)%wa(1)=1
      gen(13)%wa(1)=2

      gen(2)%wa(nparam_per_node+8)=1d0
      gen(15)%wa(nparam_per_node+16)=1d0

    !Gene-gene interactions

     !wavefront setting
      !gen(2)%w(1)=1d2  !this will spontaneously generate activator
      !gen(15)%w(1)=-1d1  !this will spontaneously generate activator
      !gen(1)%w(14)=1d0  !this will produce a ware
      !!gen(1)%w(15)=-1d0  !this is everywhere, inhibiting production of activator (in a threshold manner)
      !gen(14)%w(14)=1d0  !this is everywhere, inhibiting production of activator (in a threshold manner)
      !gen(15)%w(15)=1d0  !this is everywhere, inhibiting production of activator (in a threshold manner)
      !gen(2)%w(15)=-1d5  !this is everywhere, inhibiting production of activator (in a threshold manner)



      gen(2)%w(1)=-1d5
      gen(2)%w(1)=-1d5

      !gen(2)%w(2)=1d3  !the activated receptor directly transcripts signal
      !gen(9)%w(2)=1d3

      !gen(7)%w(2)=1d0
      !gen(8)%w(2)=1d0
      !gen(7)%w(3)=-6d-2
      !gen(8)%w(3)=-1d3
      !gen(7)%w(5)=-1d3
      !gen(8)%w(5)=-6d-2
      !gen(14)%w(7)=-1d2
      !gen(15)%w(8)=-1d2
      
      !gen(2)%w(9)=-1d2
      
      gen(3)%w(3)=1d0  !autoactivation of receptors
      gen(5)%w(5)=1d0
      gen(10)%w(10)=1d0
      
      gen(12)%w(12)=1d0 !autoactivation of adhesion molecules
      gen(13)%w(13)=1d0
      gen(14)%w(3)=1d0
      gen(15)%w(5)=1d0

      
      !gen(4)%nww=4      !activator signal activates activator receptor epi
      !gen(4)%ww(1,1)=2  
      !gen(4)%ww(1,2)=4
      !gen(4)%ww(1,3)=1d0
      !gen(4)%ww(2,1)=4  
      !gen(4)%ww(2,2)=2
      !gen(4)%ww(2,3)=1d0
      !
      !gen(4)%ww(3,1)=3  
      !gen(4)%ww(3,2)=4
      !gen(4)%ww(3,3)=1d0
      !gen(4)%ww(4,1)=4  
      !gen(4)%ww(4,2)=3
      !gen(4)%ww(4,3)=1d0
      !
      !
      !gen(6)%nww=4      !activator morphogen activates activator receptor mesench
      !gen(6)%ww(1,1)=2  
      !gen(6)%ww(1,2)=6
      !gen(6)%ww(1,3)=1d0
      !gen(6)%ww(2,1)=6  
      !gen(6)%ww(2,2)=2
      !gen(6)%ww(2,3)=1d0
      !
      !gen(6)%ww(3,1)=5
      !gen(6)%ww(3,2)=6
      !gen(6)%ww(3,3)=1d0
      !gen(6)%ww(4,1)=6  
      !gen(6)%ww(4,2)=5
      !gen(6)%ww(4,3)=1d0
      !
      !gen(11)%nww=4      !inhibitor signal activates inhibitor receptor epi
      !gen(11)%ww(1,1)=9  
      !gen(11)%ww(1,2)=11
      !gen(11)%ww(1,3)=1d0
      !gen(11)%ww(2,1)=11  
      !gen(11)%ww(2,2)=9
      !gen(11)%ww(2,3)=1d0
      !
      !gen(11)%ww(3,1)=10  
      !gen(11)%ww(3,2)=11
      !gen(11)%ww(3,3)=1d0
      !gen(11)%ww(4,1)=11  
      !gen(11)%ww(4,2)=10
      !gen(11)%ww(4,3)=1d0


    !Adhesion molecules

	ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      kadh(1,1)=5d0
      kadh(1,2)=3d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=5d0
      !kadh(3,3)=5d0
      !kadh(3,4)=2d0 ; kadh(4,3)=kadh(3,4)
      !kadh(4,4)=5d0
    end if

    !Gene expression on nodes

    
    
    !draw rectangles
    !!!!!!!!!!!!!!!!!!!!!!!!!
    jj1=1 ; jj2=ly
    ii1=1 ; ii2=lx
    do i=ii1,ii2
      do j=jj1,jj2
        k=cels(cell_grid_epi(i,j))%node(1)
        
        if(i==5.and.j==5) gex(k,2)=1d1
        
        if(i<=2.or.i>=lx-1.or.j<=2.or.j>=ly-1)then
          gex(k,1)=1d0
          !gex(k,2)=0d0 ; gex(k,9)=0d0
        end if
        do ii=1,layer
          k=cels(cell_grid_mes(i,j,ii))%node(1)
          if(i<=2.or.i>=lx-1.or.j<=2.or.j>=ly-1)then
            gex(k,1)=1d0
            !gex(k,2)=0d0 ; gex(k,9)=0d0
          end if
        end do
      end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!


   
    do i=1,nd
    
      if(node(i)%tipus<3)then
        if(node(i)%tipus==2)then
          gex(i,3)=1d0  !receptors
          gex(i,10)=1d0
          gex(i,14)=1d0
          !call random_number(a)
          !gex(i,2)=a!*0.1d0
          !gex(i,9)=a!*0.1d0
        end if
        !if(node(i)%tipus==1) gex(i,14)=1d0  !receptors
       
        gex(i,12)=1d0
        !gex(i,10)=1d0
      end if
      if(node(i)%tipus==3)then
        gex(i,5)=1d0
        gex(i,13)=1d0
        gex(i,15)=1d0
      end if
      !if(node(i)%hold>0)then
      !  gex(i,:)=0d0 !clear the hold nodes
      !  if(node(i)%tipus<3)then
      !    gex(i,12)=1d0
      !  end if
      !  if(node(i)%tipus==3)then
      !    gex(i,13)=1d0
      !  end if
      !end if

    end do


    call update_npag

    node(:)%talone=0.0d0

    !ramax=maxval(node(:)%da)*3

    
end subroutine

subroutine epidermal_organ_in_tooth_bud

!real*8::pericel,radio2,req1,req2,adjustreq,di
!integer::kk2

integer:: lx,ly
integer:: x1,x2,y1,y2,ii1,ii2,jj1,jj2
integer:: x3,x4,y3,y4,ii3,ii4,jj3,jj4,jj5,y5
integer, allocatable::cell_grid_epi(:,:),cell_grid_mes(:,:,:)
!integer:: periodic_borders,nborders
integer::geometry
integer::xlayer
integer::ring,ring2,polarized
integer::px1,px2,py1,py2,sx1,sx2,sy1,sy2,apx1,apx2,apy1,apy2

!******* #1 DEFINING SPATIAL DIMENSIONS *******

  geometry=1 !geometry=1 means hexagonal epithelium ; geometry=2 means rectangular epithelium
  
  if (geometry==1)then
	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=10    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=10   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=0.75d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    do i=1,layer
      mesradicel(i)=mradicel
    end do
    
    xlayer=1 !this will make a layer of ECM beneath the layers of mesenchymal cells
    
  elseif(geometry==2)then
    layer=1      !number of planar cell layers
    lx=16
    ly=32
    
    allocate(cell_grid_epi(lx,ly),cell_grid_mes(lx,ly,layer))

  end if  
    !periodic_borders=1
    


    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(geometry==1)then
      if(radi>0.and.radicel>0)then
        j=0
        do i=1,radi-1
          j=j+i
        end do
        nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
        nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13
  
        j=0
        do i=1,radicel-1
          j=j+i
        end do
        ncelsepi=(6*j+1)
        ndepi=nodecel*ncelsepi
      else
        ncelsepi=0
        ndepi=0
        nodecel=0
        nodecela=0
      end if
  
      if(mradi>0.and.mradicel>0.and.layer>0)then
        if(packed==1)then !if mesenchyme is packed, we use mradi differently
          j=0
          do i=1,mradi-1
            j=j+i
          end do
          nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
          nodecela=2*nodecel+1
        else
          nodecel=mradi
          nodecela=2*nodecel+1
        end if
  
        ncelsmes=0
        do k=1,layer
          j=0 !;print*,"mesradicel",mesradicel(k)
          do i=1,mesradicel(k)-1
            j=j+i
          end do
          ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
        end do
        ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
        if(radi==0.and.radicel==0)then
          nodecel=mradi
          nodecela=2*nodecel+1
        else if(nodecel<mradi)then
          nodecel=radi
          nodecela=nodecel*2+1
        end if
      else
        ndmes=0
        ncelsmes=0
      end if
      ndx=(6*j+1)*xlayer
      nd=ndepi+ndmes+ndx
      !ncels=ncelsepi+ncelsmes
      ncels=ncelsepi+ncelsmes+ndx
   print*,"ncels",ncels,ncelsepi,ncelsmes,ndx
      nda=nd+10
      ncals=ncels+10
    elseif(geometry==2)then

      ncelsepi=lx*ly
      ndepi=ncelsepi*2
      ncelsmes=ncelsepi*layer
      ndmes=ncelsmes
      nd=ndepi+ndmes
      ncels=ncelsepi+ncelsmes
      nodecel=2 ; nodecela=5
      nda=nd+10
      ncals=ncels+10
    end if
    

  !End initializing dimension parameters

  print*,"inicial nd",nd,"ncels",ncels,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=1
    !End Allocatations


    

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=0.85d0
    deltamin=1d-2
    khold=1d0
    angletor=0.05



    k_bu=5d0
    ramax=0.35d0
    k_press=5d-1
    m_xwall=zmes !tooth
    mi_xwall=0d0 !tooth

    !m_xwall=0d0  !SCALE
    !mi_xwall=-5d-1 !SCALE

    ndmax=9d4
    
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=1
    ffu(2)=1 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !buoyancy

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(14)=1 !physical boundaries (walls)
    ffu(17)=1 !volume conservation (epithelial)
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise
    ffu(24)=0


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    !if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(mod(i,2)/=0)then  !basal
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
          !node(i)%reqp=-0.05d0
          !node(i)%req=0.20d0
        !node(i)%kplast=1d0
        else                      !apical
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !node(i)%kplast=1d1
        !  !node(i)%reqp=0.05d0
        !  !node(i)%req=0.30d0
        end if
        node(i)%req=0.25d0
        node(i)%reqcr=0.25d0 !; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.125
        node(i)%da=node(i)%req*2.0; 
        node(i)%ke=5d0
        node(i)%tor=3d0
        node(i)%stor=1d0
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=1d1
        node(i)%kvol=5d1
      end do
    !end if

	!Mesenchyme
    if(layer>0)then
      do i=ndepi+1,ndepi+ndmes
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !scale rect
        !node(i)%rep=0d0 ; node(i)%repcel=1d0
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*2.0
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

	!Suprabasal cells
    if(xlayer>0)then
      do i=ndepi+ndmes+1,nd
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*2.0
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if
    
    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)

    node(:)%hold=0


  if(geometry==1)then
    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    !if(radi>0.and.radicel>0)               call epiteli_sphere1(radi,radicel,radius,di,zepi)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if
!print*,"tall0"

    !Mesenchymal layer
    ii=0
    do i=ndepi+ndmes+1,nd
      ii=ii+1
      j=i-ndx
      node(i)%x=node(j)%x ; node(i)%y=node(j)%y ; node(i)%z=-0.5d0
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
      node(i)%icel=node(j)%icel+ndx ; node(i)%tipus=3
      !node(i)%icel=ncelsepi+ncelsmes+ii ; node(i)%tipus=3

      node(i)%repcel=5d-1
      !allocate(cels(node(i)%icel)%node(nodela))
      cels(node(i)%icel)%nunodes=1
      cels(node(i)%icel)%node(1)=i
      node(i)%marge=0
      !print*,"i",i,"j",j,"icel i",node(i)%icel,"icel j",node(j)%icel
      !node(i)%border=1
    end do
    ncelsmes=ncelsmes+ndx ; ndx=0
    ndmes=ncelsmes
!print*,"tall1"
      
    j=0
    do i=1,radicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    do i=j+1,ndepi
      node(i)%hold=2 ;node(i)%repcel=1d0
    end do
    j=0
    do i=1,mradicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    k=0
    do i=1,mradicel-1
      k=k+i
    end do
    k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
      node(i)%hold=2;node(i)%repcel=1d0
      !node(i)%border=1
      !node(i)%da=node(i)%req*2.0
      !node(i)%orix=0 ; node(i)%oriy=0
      !node(i)%rep=1d1;node(i)%repcel=1d1
    end do
    do i=ndepi+ndmes/2+j+1,nd
      node(i)%hold=2;node(i)%repcel=1d0
      !node(i)%border=1
    end do
    !do i=ndepi+k+1,nd
    !  node(i)%hold=2;node(i)%repcel=1d0
    !  node(i)%border=1
    !end do

    
  elseif(geometry==2)then
	do i=1,ncels
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do


    !epithelial cells
    !origin of rectangle
    aa=node(1)%req*2
    bb=sqrt(aa**2-(aa*0.5)**2)
    a=-real(lx)*0.5*aa ; b=-real(ly)*0.5*bb
    !print*,"a",a,"b",b,"aa",aa,"bb",bb,"lx",lx,"ly",ly
    ii=0 ; jj=0
    do i=1,ly
      do j=1,lx
        ii=ii+1 ;jj=jj+1 ; cell_grid_epi(j,i)=jj
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii+1 ;node(ii)%tipus=2
        ii=ii+1
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii-1 ;node(ii)%tipus=1
        cels(jj)%node(1)=ii-1 ; cels(jj)%node(2)=ii ;cels(jj)%nunodes=2; cels(jj)%ctipus=1
        if(i==1.or.i==ly.or.j==1.or.j==lx)then
          node(ii)%hold=2;node(ii-1)%hold=2 ;
          node(ii)%repcel=1d1;node(ii-1)%repcel=1d1;
          !node(ii)%border=1;node(ii-1)%border=1;
        end if
      end do
    end do

    !mesenchymals
    do k=1,layer
      !print*,"k",k
      do i=1,ly
        do j=1,lx
          ii=ii+1 ; jj=jj+1 ; cell_grid_mes(j,i,k)=jj
          if(mod(i,2)==0)then
            node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          else
            node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          end if
          node(ii)%icel=jj ; node(ii)%altre=0 ;node(ii)%tipus=3
          cels(jj)%node(1)=ii ;cels(jj)%nunodes=1 ; cels(jj)%ctipus=3
          if(i==1 .or. i==ly .or. j==1 .or. j==lx)then
            node(ii)%hold=2 ;node(ii)%repcel=1d1 ; !node(ii)%border=2 ;
          end if
          !if(i==1 .or. i==ly .or. j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
          !if(j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
        end do
      end do
    end do
  end if

  !m_xwall=maxval(node(1:nd)%x)+0.01
  !mi_xwall=minval(node(1:nd)%x)-0.01
  !m_ywall=maxval(node(1:nd)%y)+0.01
  !mi_ywall=minval(node(1:nd)%y)-0.01
  
!print*,"tall2"

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
  !node(:)%hold=2
  !node(:)%hold=0  

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    !if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    !end if
    !Mesenchymal
    if(layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=25
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=0d0 ; gen(1)%mu=5d-1 ;gen(1)%name="epithelial basal marker"
      gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=5d-1 ;gen(2)%name="epithelial suprabasal marker"
      gen(3)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=5d-1 ;gen(3)%name="mesenchymal marker"
      gen(4)%kindof=1 ; gen(4)%diffu=0d0 ; gen(4)%mu=5d-1 ;gen(4)%name="epi. placode marker"
      gen(5)%kindof=1 ; gen(5)%diffu=0d0 ; gen(5)%mu=5d-1 ;gen(5)%name="epi. interplacode marker"
      gen(6)%kindof=1 ; gen(6)%diffu=0d0 ; gen(6)%mu=5d-1 ;gen(6)%name="mes. placode marker"
      gen(7)%kindof=1 ; gen(7)%diffu=0d0 ; gen(7)%mu=5d-1 ;gen(7)%name="mes. interplacode marker"

      gen(8)%kindof=1 ; gen(8)%diffu=1d0 ; gen(8)%mu=2.5d-1 ;gen(8)%name="epi. basal adhesion molecule"
      gen(9)%kindof=1 ; gen(9)%diffu=1d0 ; gen(9)%mu=5d-1 ;gen(9)%name="epi. suprabasal adhesion molecule"
      gen(10)%kindof=1 ; gen(10)%diffu=0d0 ; gen(10)%mu=5d-1 ;gen(10)%name="mesench. adhesion molecule"

      gen(11)%kindof=1 ; gen(11)%diffu=0d0 ; gen(11)%mu=5d-1 ;gen(11)%name="epi. signalling centre marker"
      gen(12)%kindof=4 ; gen(12)%diffu=5d-1 ; gen(12)%mu=5d-1 ;gen(12)%name="epi. signal"

      gen(13)%kindof=2 ; gen(13)%diffu=1d0 ; gen(13)%mu=2.5d-1 ;gen(13)%name="epi. basal receptor inactive"
      gen(14)%kindof=2 ; gen(14)%diffu=1d0 ; gen(14)%mu=5d-1 ;gen(14)%name="epi. suprabasal receptor inactive"
      gen(15)%kindof=2 ; gen(15)%diffu=1d0 ; gen(15)%mu=5d-1 ;gen(15)%name="mesench. receptor inactive"

      gen(16)%kindof=8 ; gen(16)%diffu=1d0 ; gen(16)%mu=2.5d-1 ;gen(16)%name="epi. basal receptor active"
      gen(16)%npre=2 ; allocate(gen(16)%pre(gen(16)%npre)) ; gen(16)%pre(1)=12 ; gen(16)%pre(2)=13
      gen(16)%npost=2 ; allocate(gen(16)%post(gen(16)%npost)) ; gen(16)%post(1)=12 ; gen(16)%post(2)=13

      gen(17)%kindof=8 ; gen(17)%diffu=1d0 ; gen(17)%mu=5d-1 ;gen(17)%name="epi. suprabasal receptor active"
      gen(17)%npre=2 ; allocate(gen(17)%pre(gen(17)%npre)) ; gen(17)%pre(1)=12 ; gen(17)%pre(2)=14
      gen(17)%npost=2 ; allocate(gen(17)%post(gen(17)%npost)) ; gen(17)%post(1)=12 ; gen(17)%post(2)=14

      gen(18)%kindof=8 ; gen(18)%diffu=1d0 ; gen(18)%mu=5d-1 ;gen(18)%name="mesench. receptor active"
      gen(18)%npre=2 ; allocate(gen(18)%pre(gen(18)%npre)) ; gen(18)%pre(1)=12 ; gen(18)%pre(2)=15
      gen(18)%npost=2 ; allocate(gen(18)%post(gen(18)%npost)) ; gen(18)%post(1)=12 ; gen(18)%post(2)=15

      gen(19)%kindof=1 ; gen(19)%diffu=0d0 ; gen(19)%mu=5d-1 ;gen(19)%name="epi. basal effector"
      gen(20)%kindof=1 ; gen(20)%diffu=0d0 ; gen(20)%mu=5d-1 ;gen(20)%name="epi. suprabasal effector"
      gen(21)%kindof=1 ; gen(21)%diffu=0d0 ; gen(21)%mu=5d-1 ;gen(21)%name="mesenchymal effector"

      gen(22)%kindof=1 ; gen(22)%diffu=0d0 ; gen(22)%mu=5d-1 ;gen(22)%name="signalling centre adh. mol."
      gen(23)%kindof=1 ; gen(23)%diffu=0d0 ; gen(23)%mu=5d-1 ;gen(23)%name="distal (posterior) gene product"
      
      gen(24)%kindof=1 ; gen(24)%diffu=1d0 ; gen(24)%mu=2.5d-1 ;gen(24)%name="s2 IEE adhesion molecule."
      gen(25)%kindof=1 ; gen(25)%diffu=0d0 ; gen(25)%mu=5d-1 ;gen(25)%name="s2 dermal papilla adhesion molecule"
    
 
    !Gene-behavior interactions
    
   
      gen(8)%wa(1)=1
      gen(9)%wa(1)=2
      gen(10)%wa(1)=3
      gen(22)%wa(1)=4
      gen(24)%wa(1)=5
      gen(25)%wa(1)=6

      !gen(19)%wa(1)=5
      !gen(20)%wa(1)=6
      !gen(21)%wa(1)=7
      
      gen(22)%wa(10)=1d0
      
      gen(11)%wa(nparam_per_node+2)=-5d5

      gen(19)%wa(nparam_per_node+2)=0.0d0
      gen(20)%wa(nparam_per_node+2)=0.0d0
      gen(21)%wa(nparam_per_node+2)=0.0d0
      
      !growth parameters for inward bud
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=3.2d-1
      !gen(21)%wa(nparam_per_node+2)=5d-2
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !growth parameters for outward bud
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=5.0d-3!5.0d-2
      !gen(21)%wa(nparam_per_node+2)=3.2d-1
      !gen(23)%wa(nparam_per_node+2)=-5.0d0
      !!!!!!!!!!!!!!!!!!!!!!!!!!   
      !growth parameters for rectangle scale
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=0.0d0
      !gen(21)%wa(nparam_per_node+2)=1.0d-3
      !gen(23)%wa(nparam_per_node+2)=-7d1
      !gen(6)%wa(nparam_per_node+2)=1.0d-1
      !!!!!!!!!!!!!!!!!!!!!!!!!!   
      
    !Gene-gene interactions

      !structural genes
      !auto-maintenance
      gen(1)%w(1)=1d0
      gen(2)%w(2)=1d0
      gen(3)%w(3)=1d0
      gen(4)%w(4)=1d0
      gen(5)%w(5)=1d0
      gen(6)%w(6)=1d0
      gen(7)%w(7)=1d0
      gen(11)%w(11)=1d0
      !maintenance of stable gene expression
      gen(13)%w(4)=1d0
      gen(14)%w(4)=1d0
      gen(15)%w(6)=1d0
      gen(8)%w(1)=1d0
      gen(9)%w(2)=1d0
      gen(10)%w(3)=1d0
      gen(22)%w(11)=1d0
      gen(8)%w(11)=-1d5
      gen(9)%w(11)=-1d5
      !!!!!!!!!!!!!!!!
      gen(12)%w(11)=1d0
      gen(13)%w(11)=-1d5
      gen(14)%w(11)=-1d5
      gen(14)%w(1)=-1d5
      gen(13)%w(2)=-1d5
      !gen(24)%w(11)=1d0 !s2
      !!!!!!!!!!!!!!!!! effector activation thresholds
      gen(19)%w(1)=-0.05d0
      gen(20)%w(2)=-0.05d0
      gen(21)%w(3)=-0.05d0
      !gen(19)%w(1)=0d0
      !gen(20)%w(2)=0d0
      !gen(21)%w(3)=0d0
      gen(19)%w(16)=0d0   !tooth
      gen(20)%w(17)=0d0
      gen(21)%w(18)=0d0
      !gen(19)%w(16)=1d2    !scale
      !gen(20)%w(17)=1d2
      !gen(21)%w(18)=1d1
      !gen(19)%w(16)=1d1   !tooth
      !gen(20)%w(17)=1d1
      !gen(21)%w(18)=1d1
      !gen(19)%w(16)=1d2    !rectangle scale
      !gen(20)%w(17)=1d2
      !gen(21)%w(18)=1d1

      !gen(24)%w(19)=1d0 !s2
      !gen(25)%w(21)=1d0 !s2


      !effector supression of adhesion molecules
      !gen(8)%w(19)=-1d1
      !gen(9)%w(20)=-1d1
      !gen(10)%w(21)=-1d1

      !!!!!!!!!!!!!!!! a gene expressed in the posterior half of the placode supresses signalling
      gen(23)%w(23)=1d0
      !gen(13)%w(23)=-1d5
      !gen(14)%w(23)=-1d5

    


      !signal-receptor interactions
      !epithelial basal binding
      gen(16)%nww=4
      gen(16)%ww(1,1)=12
      gen(16)%ww(1,2)=16
      gen(16)%ww(1,3)=1d0
      gen(16)%ww(2,1)=16
      gen(16)%ww(2,2)=12
      gen(16)%ww(2,3)=1d0
      gen(16)%ww(3,1)=13
      gen(16)%ww(3,2)=16
      gen(16)%ww(3,3)=1d0
      gen(16)%ww(4,1)=16
      gen(16)%ww(4,2)=13
      gen(16)%ww(4,3)=1d0
      !epithelial suprabasal receptor binding
      gen(17)%nww=4
      gen(17)%ww(1,1)=12
      gen(17)%ww(1,2)=17
      gen(17)%ww(1,3)=1d0
      gen(17)%ww(2,1)=17
      gen(17)%ww(2,2)=12
      gen(17)%ww(2,3)=1d0
      gen(17)%ww(3,1)=14
      gen(17)%ww(3,2)=15
      gen(17)%ww(3,3)=1d0
      gen(17)%ww(4,1)=17
      gen(17)%ww(4,2)=14
      gen(17)%ww(4,3)=1d0
      !mesenchymal receptor binding
      gen(18)%nww=4
      gen(18)%ww(1,1)=12
      gen(18)%ww(1,2)=18
      gen(18)%ww(1,3)=1d0
      gen(18)%ww(2,1)=18
      gen(18)%ww(2,2)=12
      gen(18)%ww(2,3)=1d0
      gen(18)%ww(3,1)=15
      gen(18)%ww(3,2)=18
      gen(18)%ww(3,3)=1d0
      gen(18)%ww(4,1)=18
      gen(18)%ww(4,2)=15
      gen(18)%ww(4,3)=1d0
      


    !Adhesion molecules

	ntipusadh=6
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      !as default, 1 and 2 should act as the same (basal and suprabasal)
      kadh(1,1)=1d0

      kadh(1,2)=1d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=1d0

      kadh(1,3)=5d-1 ; kadh(3,1)=kadh(1,3)
      !rectangle scale
      !kadh(1,3)=3d0 ; kadh(3,1)=kadh(1,3)
      kadh(2,3)=0d0 ; kadh(3,2)=kadh(2,3)
      kadh(3,3)=1d0

      kadh(1,4)=5d-1 ; kadh(4,1)=kadh(1,4)
      kadh(2,4)=5d-1 ; kadh(4,2)=kadh(2,4)
      kadh(3,4)=5d-1 ; kadh(4,3)=kadh(3,4)
      kadh(4,4)=1d1 !tooth
      !kadh(4,4)=5d0  !scale

      kadh(1,5)=0d0 ; kadh(5,1)=kadh(1,5)
      kadh(2,5)=0d0 ; kadh(5,2)=kadh(2,5)
      kadh(3,5)=0d0 ; kadh(5,3)=kadh(3,5)
      kadh(4,5)=0d0 ; kadh(5,4)=kadh(4,5)
      kadh(5,5)=0d0

      kadh(1,6)=0d0 ; kadh(6,1)=kadh(1,6)
      kadh(2,6)=0d0 ; kadh(6,2)=kadh(2,6)
      kadh(3,6)=0d0 ; kadh(6,3)=kadh(3,6)
      kadh(4,6)=0d0 ; kadh(6,4)=kadh(4,6)
      kadh(5,6)=1d1 ; kadh(6,5)=kadh(5,6)
      kadh(6,6)=1d1

      !kadh(1,7)=5d-1 ; kadh(7,1)=kadh(1,7)
      !kadh(2,7)=0d0 ; kadh(7,2)=kadh(2,7)
      !kadh(3,7)=1d0 ; kadh(7,3)=kadh(3,7)
      !kadh(4,7)=5d-1 ; kadh(7,4)=kadh(4,7)
      !kadh(5,7)=5d-1 ; kadh(7,5)=kadh(5,7)
      !kadh(6,7)=0d0 ; kadh(7,6)=kadh(6,7)
      !kadh(7,7)=1d0

    end if

    !Gene expression on nodes

    
    !if(geometry==1)then
    !  gex(1,2)=1d0
    !  gex(ndepi+1,2)=1d0
    !  !gex(ndepi+1,2)=1d0
    !end if
    
    if(geometry==2)then
    !draw rectangles
    !!!!!!!!!!!!!!!!!!!!!!!!!
    jj1=1 ; jj2=ly
    ii1=1 ; ii2=lx
    !placode range
    px1=3 ; px2=lx-(px1-1)
    py1=3 ; py2=ly-(py1-1)
    !signalling centre
    sx1=7 ; sx2=lx-(sx1-1)
    sy1=5 ; sy2=ly-(sy1-1)
    !AP proliferation supressor
    apx1=7 ; apx2=lx
    apy1=1 ; apy2=ly

    do i=ii1,ii2
      do j=jj1,jj2
        !epithelium
        k=cels(cell_grid_epi(i,j))%node(1)

        gex(k,1)=1d0
        gex(k,8)=1d0
        !borders
        if(i<=1.or.i>=lx.or.j<=1.or.j>=ly)then
          node(k)%hold=3
          node(node(k)%altre)%hold=3
        end if
        !placode
        if((i>px1.and.i<px2).and.(j>py1.and.j<py2))then
          gex(k,4)=1d0
          gex(k,13)=1d0
        else
          gex(k,5)=1d0
        end if
        !signalling centre
        if((i>sx1.and.i<sx2).and.(j>sy1.and.j<sy2))then
          gex(k,11)=1d0 ; gex(node(k)%altre,11)=1d0
          gex(k,22)=1d0 ; gex(node(k)%altre,22)=1d0
          gex(k,13)=0d0
        end if
        !AP proliferation supressor
        if((i>apx1.and.i<apx2).and.(j>apy1.and.j<apy2))then
          gex(k,23)=1d0
        end if


        !mesench
        do ii=1,layer
          k=cels(cell_grid_mes(i,j,ii))%node(1)
          gex(k,3)=1d0
          gex(k,10)=1d0
          !borders
          if(i<=1.or.i>=lx.or.j<=1.or.j>=ly)then
            node(k)%hold=3
          end if
          !placode
          if((i>px1.and.i<px2).and.(j>py1.and.j<py2))then
            gex(k,6)=1d0
            gex(k,15)=1d0
          else
            gex(k,7)=1d0
          end if
        end do
      end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!
    end if

  if(geometry==1)then  
    !ring=2  !scale
    ring=2 !tooth

    ring2=7

    polarized=0 !tooth
    !polarized=1 !scale

    do i=1,nd
      if(node(i)%tipus<3)then
        if(node(i)%tipus==2)then
          !epithelial basal layer, basal side
          gex(i,1)=1d0  !receptors

          d=sqrt(node(i)%x**2+node(i)%y**2)
          !setting expression on placode and interplacode
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,4)=1d0
            if(polarized==0.or.node(i)%x<=0.15)then
              gex(i,13)=1d0 ; gex(node(i)%altre,13)=1d0
            else
              gex(i,23)=1d0
              gex(i,13)=1d0 ; gex(node(i)%altre,13)=1d0
            end if
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,5)=1d0
          end if
          !!!!!!!!!!!!!
          !signalling centre
          if(d-2*node(i)%req*(ring-1)<epsilod)then
            node(i)%marge=0 ; node(node(i)%altre)%marge=0
            gex(i,11)=1d0 ; gex(node(i)%altre,11)=1d0
            gex(i,22)=1d0 ; gex(node(i)%altre,22)=1d0
            gex(i,13)=0d0 ; gex(node(i)%altre,13)=0d0
          end if

        end if
        !epithelial basal layer, whole layer       
        gex(i,8)=1d0

      end if
      if(node(i)%tipus==3)then
        if(node(i)%z<=0)then
          !mesenchymal layer
          gex(i,3)=1d0
          gex(i,10)=1d0

          !setting expression on placode and interplacode
          d=sqrt(node(i)%x**2+node(i)%y**2)
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,6)=1d0
            gex(i,15)=1d0
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,7)=1d0
          end if
          !!!!!!!!!!!!!

        else
          !epithelial suprabasal layer
          gex(i,2)=1d0
          gex(i,9)=1d0

          !setting expression on placode and interplacode
          d=sqrt(node(i)%x**2+node(i)%y**2)
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,4)=1d0
            if(polarized==0.or.node(i)%x<=0.15)then
              gex(i,14)=1d0
            else
              gex(i,23)=1d0
              gex(i,14)=1d0
            end if
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,5)=1d0
          end if
          !!!!!!!!!!!!!
          !signalling centre
          !if(d-2*node(i)%req*(ring-1)<epsilod)then
          !  gex(i,11)=1d0
          !  gex(i,22)=1d0
          !  gex(i,14)=0d0 !; gex(node(i)%altre,13)=0d0
          !end if


        end if
      end if


    end do
    end if

    call update_npag

    node(:)%talone=0.0d0

    !ramax=maxval(node(:)%da)*3

    
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine epidermal_organ_in_tooth_bud_overall

!real*8::pericel,radio2,req1,req2,adjustreq,di
!integer::kk2

integer:: lx,ly
integer:: x1,x2,y1,y2,ii1,ii2,jj1,jj2
integer:: x3,x4,y3,y4,ii3,ii4,jj3,jj4,jj5,y5
integer, allocatable::cell_grid_epi(:,:),cell_grid_mes(:,:,:)
!integer:: periodic_borders,nborders
integer::geometry
integer::xlayer
integer::ring,ring2,polarized
integer::px1,px2,py1,py2,sx1,sx2,sy1,sy2,apx1,apx2,apy1,apy2

!******* #1 DEFINING SPATIAL DIMENSIONS *******

  geometry=1 !geometry=1 means hexagonal epithelium ; geometry=2 means rectangular epithelium
  
  if (geometry==1)then
	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=7    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=7   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=0.75d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    do i=1,layer
      mesradicel(i)=mradicel
    end do
    
    xlayer=1 !this will make a layer of ECM beneath the layers of mesenchymal cells
    
  elseif(geometry==2)then
    layer=1      !number of planar cell layers
    lx=16
    ly=32
    
    allocate(cell_grid_epi(lx,ly),cell_grid_mes(lx,ly,layer))

  end if  
    !periodic_borders=1
    


    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(geometry==1)then
      if(radi>0.and.radicel>0)then
        j=0
        do i=1,radi-1
          j=j+i
        end do
        nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
        nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13
  
        j=0
        do i=1,radicel-1
          j=j+i
        end do
        ncelsepi=(6*j+1)
        ndepi=nodecel*ncelsepi
      else
        ncelsepi=0
        ndepi=0
        nodecel=0
        nodecela=0
      end if
  
      if(mradi>0.and.mradicel>0.and.layer>0)then
        if(packed==1)then !if mesenchyme is packed, we use mradi differently
          j=0
          do i=1,mradi-1
            j=j+i
          end do
          nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
          nodecela=2*nodecel+1
        else
          nodecel=mradi
          nodecela=2*nodecel+1
        end if
  
        ncelsmes=0
        do k=1,layer
          j=0 !;print*,"mesradicel",mesradicel(k)
          do i=1,mesradicel(k)-1
            j=j+i
          end do
          ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
        end do
        ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
        if(radi==0.and.radicel==0)then
          nodecel=mradi
          nodecela=2*nodecel+1
        else if(nodecel<mradi)then
          nodecel=radi
          nodecela=nodecel*2+1
        end if
      else
        ndmes=0
        ncelsmes=0
      end if
      ndx=(6*j+1)*xlayer
      nd=ndepi+ndmes+ndx
      !ncels=ncelsepi+ncelsmes
      ncels=ncelsepi+ncelsmes+ndx
   print*,"ncels",ncels,ncelsepi,ncelsmes,ndx
      nda=nd+10
      ncals=ncels+10
    elseif(geometry==2)then

      ncelsepi=lx*ly
      ndepi=ncelsepi*2
      ncelsmes=ncelsepi*layer
      ndmes=ncelsmes
      nd=ndepi+ndmes
      ncels=ncelsepi+ncelsmes
      nodecel=2 ; nodecela=5
      nda=nd+10
      ncals=ncels+10
    end if
    

  !End initializing dimension parameters

  print*,"inicial nd",nd,"ncels",ncels,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=1
    !End Allocatations


    

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=0.85d0
    deltamin=1d-2
    khold=1d0
    angletor=0.05



    k_bu=5d0
    ramax=0.35d0
    k_press=5d-1
    m_xwall=zmes !tooth
    mi_xwall=0d0 !tooth

    !m_xwall=0d0  !SCALE
    !mi_xwall=-5d-1 !SCALE

    ndmax=9d4
    
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=1
    ffu(2)=1 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !buoyancy

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(14)=1 !physical boundaries (walls)
    ffu(17)=1 !volume conservation (epithelial)
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise
    ffu(24)=0


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    !if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(mod(i,2)/=0)then  !basal
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
          !node(i)%reqp=-0.05d0
          !node(i)%req=0.20d0
        !node(i)%kplast=1d0
        else                      !apical
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !node(i)%kplast=1d1
        !  !node(i)%reqp=0.05d0
        !  !node(i)%req=0.30d0
        end if
        node(i)%req=0.25d0
        node(i)%reqcr=0.25d0 !; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.125
        node(i)%da=node(i)%req*2.0; 
        node(i)%ke=5d0
        node(i)%tor=3d0
        node(i)%stor=1d0
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=1d1
        node(i)%kvol=5d1
      end do
    !end if

	!Mesenchyme
    if(layer>0)then
      do i=ndepi+1,ndepi+ndmes
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !scale rect
        !node(i)%rep=0d0 ; node(i)%repcel=1d0
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*2.0
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

	!Suprabasal cells
    if(xlayer>0)then
      do i=ndepi+ndmes+1,nd
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*2.0
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if
    
    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)

    node(:)%hold=0


  if(geometry==1)then
    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    !if(radi>0.and.radicel>0)               call epiteli_sphere1(radi,radicel,radius,di,zepi)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if
!print*,"tall0"

    !Mesenchymal layer
    ii=0
    do i=ndepi+ndmes+1,nd
      ii=ii+1
      j=i-ndx
      node(i)%x=node(j)%x ; node(i)%y=node(j)%y ; node(i)%z=-0.5d0
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
      node(i)%icel=node(j)%icel+ndx ; node(i)%tipus=3
      !node(i)%icel=ncelsepi+ncelsmes+ii ; node(i)%tipus=3

      node(i)%repcel=5d-1
      !allocate(cels(node(i)%icel)%node(nodela))
      cels(node(i)%icel)%nunodes=1
      cels(node(i)%icel)%node(1)=i
      node(i)%marge=0
      !print*,"i",i,"j",j,"icel i",node(i)%icel,"icel j",node(j)%icel
      !node(i)%border=1
    end do
    ncelsmes=ncelsmes+ndx ; ndx=0
    ndmes=ncelsmes
!print*,"tall1"
      
    j=0
    do i=1,radicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    do i=j+1,ndepi
      node(i)%hold=2 ;node(i)%repcel=1d0
    end do
    j=0
    do i=1,mradicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    k=0
    do i=1,mradicel-1
      k=k+i
    end do
    k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
      node(i)%hold=2;node(i)%repcel=1d0
      !node(i)%border=1
      !node(i)%da=node(i)%req*2.0
      !node(i)%orix=0 ; node(i)%oriy=0
      !node(i)%rep=1d1;node(i)%repcel=1d1
    end do
    do i=ndepi+ndmes/2+j+1,nd
      node(i)%hold=2;node(i)%repcel=1d0
      !node(i)%border=1
    end do
    !do i=ndepi+k+1,nd
    !  node(i)%hold=2;node(i)%repcel=1d0
    !  node(i)%border=1
    !end do

    
  elseif(geometry==2)then
	do i=1,ncels
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do


    !epithelial cells
    !origin of rectangle
    aa=node(1)%req*2
    bb=sqrt(aa**2-(aa*0.5)**2)
    a=-real(lx)*0.5*aa ; b=-real(ly)*0.5*bb
    !print*,"a",a,"b",b,"aa",aa,"bb",bb,"lx",lx,"ly",ly
    ii=0 ; jj=0
    do i=1,ly
      do j=1,lx
        ii=ii+1 ;jj=jj+1 ; cell_grid_epi(j,i)=jj
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii+1 ;node(ii)%tipus=2
        ii=ii+1
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii-1 ;node(ii)%tipus=1
        cels(jj)%node(1)=ii-1 ; cels(jj)%node(2)=ii ;cels(jj)%nunodes=2; cels(jj)%ctipus=1
        if(i==1.or.i==ly.or.j==1.or.j==lx)then
          node(ii)%hold=2;node(ii-1)%hold=2 ;
          node(ii)%repcel=1d1;node(ii-1)%repcel=1d1;
          !node(ii)%border=1;node(ii-1)%border=1;
        end if
      end do
    end do

    !mesenchymals
    do k=1,layer
      !print*,"k",k
      do i=1,ly
        do j=1,lx
          ii=ii+1 ; jj=jj+1 ; cell_grid_mes(j,i,k)=jj
          if(mod(i,2)==0)then
            node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          else
            node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          end if
          node(ii)%icel=jj ; node(ii)%altre=0 ;node(ii)%tipus=3
          cels(jj)%node(1)=ii ;cels(jj)%nunodes=1 ; cels(jj)%ctipus=3
          if(i==1 .or. i==ly .or. j==1 .or. j==lx)then
            node(ii)%hold=2 ;node(ii)%repcel=1d1 ; !node(ii)%border=2 ;
          end if
          !if(i==1 .or. i==ly .or. j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
          !if(j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
        end do
      end do
    end do
  end if

  !m_xwall=maxval(node(1:nd)%x)+0.01
  !mi_xwall=minval(node(1:nd)%x)-0.01
  !m_ywall=maxval(node(1:nd)%y)+0.01
  !mi_ywall=minval(node(1:nd)%y)-0.01
  
!print*,"tall2"

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
  !node(:)%hold=2
  !node(:)%hold=0  

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    !if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    !end if
    !Mesenchymal
    if(layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=25
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=0d0 ; gen(1)%mu=5d-1 ;gen(1)%name="epithelial basal marker"
      gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=5d-1 ;gen(2)%name="epithelial suprabasal marker"
      gen(3)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=5d-1 ;gen(3)%name="mesenchymal marker"
      gen(4)%kindof=1 ; gen(4)%diffu=0d0 ; gen(4)%mu=5d-1 ;gen(4)%name="epi. placode marker"
      gen(5)%kindof=1 ; gen(5)%diffu=0d0 ; gen(5)%mu=5d-1 ;gen(5)%name="epi. interplacode marker"
      gen(6)%kindof=1 ; gen(6)%diffu=0d0 ; gen(6)%mu=5d-1 ;gen(6)%name="mes. placode marker"
      gen(7)%kindof=1 ; gen(7)%diffu=0d0 ; gen(7)%mu=5d-1 ;gen(7)%name="mes. interplacode marker"

      gen(8)%kindof=1 ; gen(8)%diffu=0d0 ; gen(8)%mu=5d-1 ;gen(8)%name="epi. basal adhesion molecule"
      gen(9)%kindof=1 ; gen(9)%diffu=0d0 ; gen(9)%mu=5d-1 ;gen(9)%name="epi. suprabasal adhesion molecule"
      gen(10)%kindof=1 ; gen(10)%diffu=0d0 ; gen(10)%mu=5d-1 ;gen(10)%name="mesench. adhesion molecule"

      gen(11)%kindof=1 ; gen(11)%diffu=0d0 ; gen(11)%mu=5d-1 ;gen(11)%name="epi. signalling centre marker"
      gen(12)%kindof=4 ; gen(12)%diffu=5d-1 ; gen(12)%mu=5d-1 ;gen(12)%name="epi. signal"

      gen(13)%kindof=2 ; gen(13)%diffu=0d0 ; gen(13)%mu=5d-1 ;gen(13)%name="epi. basal receptor inactive"
      gen(14)%kindof=2 ; gen(14)%diffu=0d0 ; gen(14)%mu=5d-1 ;gen(14)%name="epi. suprabasal receptor inactive"
      gen(15)%kindof=2 ; gen(15)%diffu=0d0 ; gen(15)%mu=5d-1 ;gen(15)%name="mesench. receptor inactive"

      gen(16)%kindof=8 ; gen(16)%diffu=0d0 ; gen(16)%mu=5d-1 ;gen(16)%name="epi. basal receptor active"
      gen(16)%npre=2 ; allocate(gen(16)%pre(gen(16)%npre)) ; gen(16)%pre(1)=12 ; gen(16)%pre(2)=13
      gen(16)%npost=2 ; allocate(gen(16)%post(gen(16)%npost)) ; gen(16)%post(1)=12 ; gen(16)%post(2)=13

      gen(17)%kindof=8 ; gen(17)%diffu=0d0 ; gen(17)%mu=5d-1 ;gen(17)%name="epi. suprabasal receptor active"
      gen(17)%npre=2 ; allocate(gen(17)%pre(gen(17)%npre)) ; gen(17)%pre(1)=12 ; gen(17)%pre(2)=14
      gen(17)%npost=2 ; allocate(gen(17)%post(gen(17)%npost)) ; gen(17)%post(1)=12 ; gen(17)%post(2)=14

      gen(18)%kindof=8 ; gen(18)%diffu=0d0 ; gen(18)%mu=5d-1 ;gen(18)%name="mesench. receptor active"
      gen(18)%npre=2 ; allocate(gen(18)%pre(gen(18)%npre)) ; gen(18)%pre(1)=12 ; gen(18)%pre(2)=15
      gen(18)%npost=2 ; allocate(gen(18)%post(gen(18)%npost)) ; gen(18)%post(1)=12 ; gen(18)%post(2)=15

      gen(19)%kindof=1 ; gen(19)%diffu=0d0 ; gen(19)%mu=5d-1 ;gen(19)%name="epi. basal effector"
      gen(20)%kindof=1 ; gen(20)%diffu=0d0 ; gen(20)%mu=5d-1 ;gen(20)%name="epi. suprabasal effector"
      gen(21)%kindof=1 ; gen(21)%diffu=0d0 ; gen(21)%mu=5d-1 ;gen(21)%name="mesenchymal effector"

      gen(22)%kindof=1 ; gen(22)%diffu=0d0 ; gen(22)%mu=5d-1 ;gen(22)%name="signalling centre adh. mol."
      gen(23)%kindof=1 ; gen(23)%diffu=0d0 ; gen(23)%mu=5d-1 ;gen(23)%name="distal (posterior) gene product"
      
      gen(24)%kindof=1 ; gen(24)%diffu=0d0 ; gen(24)%mu=5d-1 ;gen(24)%name="s2 IEE adhesion molecule."
      gen(25)%kindof=1 ; gen(25)%diffu=0d0 ; gen(25)%mu=5d-1 ;gen(25)%name="s2 dermal papilla adhesion molecule"
    
 
    !Gene-behavior interactions
    
   
      gen(8)%wa(1)=1
      gen(9)%wa(1)=2
      gen(10)%wa(1)=3
      gen(22)%wa(1)=4
      gen(24)%wa(1)=5
      gen(25)%wa(1)=6

      !gen(19)%wa(1)=5
      !gen(20)%wa(1)=6
      !gen(21)%wa(1)=7
      
      gen(22)%wa(10)=1d0
      
      gen(11)%wa(nparam_per_node+2)=-5d5

      gen(19)%wa(nparam_per_node+2)=2.0d-1/2d0
      gen(20)%wa(nparam_per_node+2)=2.0d-1
      gen(21)%wa(nparam_per_node+2)=5.0d-2


      !gen(13)%wa(nparam_per_node+2)=2.0d-1/2d0
      !gen(14)%wa(nparam_per_node+2)=2.0d-1
      !gen(15)%wa(nparam_per_node+2)=5.0d-2

      
      !growth parameters for inward bud
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=3.2d-1
      !gen(21)%wa(nparam_per_node+2)=5d-2
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !growth parameters for outward bud
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=5.0d-3!5.0d-2
      !gen(21)%wa(nparam_per_node+2)=3.2d-1
      !gen(23)%wa(nparam_per_node+2)=-5.0d0
      !!!!!!!!!!!!!!!!!!!!!!!!!!   
      !growth parameters for rectangle scale
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=0.0d0
      !gen(21)%wa(nparam_per_node+2)=1.0d-3
      !gen(23)%wa(nparam_per_node+2)=-7d1
      !gen(6)%wa(nparam_per_node+2)=1.0d-1
      !!!!!!!!!!!!!!!!!!!!!!!!!!   
      
    !Gene-gene interactions

      !structural genes
      !auto-maintenance
      gen(1)%w(1)=1d0
      gen(2)%w(2)=1d0
      gen(3)%w(3)=1d0
      gen(4)%w(4)=1d0
      gen(5)%w(5)=1d0
      gen(6)%w(6)=1d0
      gen(7)%w(7)=1d0
      gen(11)%w(11)=1d0
      !maintenance of stable gene expression
      gen(13)%w(4)=1d0
      gen(14)%w(4)=1d0
      gen(15)%w(6)=1d0
      gen(8)%w(1)=1d0
      gen(9)%w(2)=1d0
      gen(10)%w(3)=1d0
      gen(22)%w(11)=1d0
      gen(8)%w(11)=-1d5
      gen(9)%w(11)=-1d5
      !!!!!!!!!!!!!!!!
      gen(12)%w(11)=1d0
      gen(13)%w(11)=-1d5
      gen(14)%w(11)=-1d5
      gen(14)%w(1)=-1d5
      gen(13)%w(2)=-1d5
      !gen(24)%w(11)=1d0 !s2
      !!!!!!!!!!!!!!!!! effector activation thresholds
      gen(19)%w(1)=-0.05d0
      gen(20)%w(2)=-0.05d0
      gen(21)%w(3)=-0.05d0
      !gen(19)%w(1)=0d0
      !gen(20)%w(2)=0d0
      !gen(21)%w(3)=0d0
      gen(19)%w(16)=5d1   !tooth
      gen(20)%w(17)=5d1
      gen(21)%w(18)=5d1
      !gen(19)%w(16)=1d2    !scale
      !gen(20)%w(17)=1d2
      !gen(21)%w(18)=1d1
      !gen(19)%w(16)=1d1   !tooth
      !gen(20)%w(17)=1d1
      !gen(21)%w(18)=1d1
      !gen(19)%w(16)=1d2    !rectangle scale
      !gen(20)%w(17)=1d2
      !gen(21)%w(18)=1d1

      !gen(24)%w(19)=1d0 !s2
      !gen(25)%w(21)=1d0 !s2


      !effector supression of adhesion molecules
      !gen(8)%w(19)=-1d1
      !gen(9)%w(20)=-1d1
      !gen(10)%w(21)=-1d1

      !!!!!!!!!!!!!!!! a gene expressed in the posterior half of the placode supresses signalling
      gen(23)%w(23)=1d0
      !gen(13)%w(23)=-1d5
      !gen(14)%w(23)=-1d5

    


      !signal-receptor interactions
      !epithelial basal binding
      gen(16)%nww=4
      gen(16)%ww(1,1)=12
      gen(16)%ww(1,2)=16
      gen(16)%ww(1,3)=1d0
      gen(16)%ww(2,1)=16
      gen(16)%ww(2,2)=12
      gen(16)%ww(2,3)=1d0
      gen(16)%ww(3,1)=13
      gen(16)%ww(3,2)=16
      gen(16)%ww(3,3)=1d0
      gen(16)%ww(4,1)=16
      gen(16)%ww(4,2)=13
      gen(16)%ww(4,3)=1d0
      !epithelial suprabasal receptor binding
      gen(17)%nww=4
      gen(17)%ww(1,1)=12
      gen(17)%ww(1,2)=17
      gen(17)%ww(1,3)=1d0
      gen(17)%ww(2,1)=17
      gen(17)%ww(2,2)=12
      gen(17)%ww(2,3)=1d0
      gen(17)%ww(3,1)=14
      gen(17)%ww(3,2)=15
      gen(17)%ww(3,3)=1d0
      gen(17)%ww(4,1)=17
      gen(17)%ww(4,2)=14
      gen(17)%ww(4,3)=1d0
      !mesenchymal receptor binding
      gen(18)%nww=4
      gen(18)%ww(1,1)=12
      gen(18)%ww(1,2)=18
      gen(18)%ww(1,3)=1d0
      gen(18)%ww(2,1)=18
      gen(18)%ww(2,2)=12
      gen(18)%ww(2,3)=1d0
      gen(18)%ww(3,1)=15
      gen(18)%ww(3,2)=18
      gen(18)%ww(3,3)=1d0
      gen(18)%ww(4,1)=18
      gen(18)%ww(4,2)=15
      gen(18)%ww(4,3)=1d0
      


    !Adhesion molecules

	ntipusadh=6
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      !as default, 1 and 2 should act as the same (basal and suprabasal)
      kadh(1,1)=1d0

      kadh(1,2)=1d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=1d0

      kadh(1,3)=5d-1 ; kadh(3,1)=kadh(1,3)
      !rectangle scale
      !kadh(1,3)=3d0 ; kadh(3,1)=kadh(1,3)
      kadh(2,3)=0d0 ; kadh(3,2)=kadh(2,3)
      kadh(3,3)=1d0

      kadh(1,4)=5d-1 ; kadh(4,1)=kadh(1,4)
      kadh(2,4)=5d-1 ; kadh(4,2)=kadh(2,4)
      kadh(3,4)=5d-1 ; kadh(4,3)=kadh(3,4)
      kadh(4,4)=1d1 !tooth
      !kadh(4,4)=5d0  !scale

      kadh(1,5)=0d0 ; kadh(5,1)=kadh(1,5)
      kadh(2,5)=0d0 ; kadh(5,2)=kadh(2,5)
      kadh(3,5)=0d0 ; kadh(5,3)=kadh(3,5)
      kadh(4,5)=0d0 ; kadh(5,4)=kadh(4,5)
      kadh(5,5)=0d0

      kadh(1,6)=0d0 ; kadh(6,1)=kadh(1,6)
      kadh(2,6)=0d0 ; kadh(6,2)=kadh(2,6)
      kadh(3,6)=0d0 ; kadh(6,3)=kadh(3,6)
      kadh(4,6)=0d0 ; kadh(6,4)=kadh(4,6)
      kadh(5,6)=1d1 ; kadh(6,5)=kadh(5,6)
      kadh(6,6)=1d1

      !kadh(1,7)=5d-1 ; kadh(7,1)=kadh(1,7)
      !kadh(2,7)=0d0 ; kadh(7,2)=kadh(2,7)
      !kadh(3,7)=1d0 ; kadh(7,3)=kadh(3,7)
      !kadh(4,7)=5d-1 ; kadh(7,4)=kadh(4,7)
      !kadh(5,7)=5d-1 ; kadh(7,5)=kadh(5,7)
      !kadh(6,7)=0d0 ; kadh(7,6)=kadh(6,7)
      !kadh(7,7)=1d0

    end if

    !Gene expression on nodes

    
    !if(geometry==1)then
    !  gex(1,2)=1d0
    !  gex(ndepi+1,2)=1d0
    !  !gex(ndepi+1,2)=1d0
    !end if
    
    if(geometry==2)then
    !draw rectangles
    !!!!!!!!!!!!!!!!!!!!!!!!!
    jj1=1 ; jj2=ly
    ii1=1 ; ii2=lx
    !placode range
    px1=3 ; px2=lx-(px1-1)
    py1=3 ; py2=ly-(py1-1)
    !signalling centre
    sx1=7 ; sx2=lx-(sx1-1)
    sy1=5 ; sy2=ly-(sy1-1)
    !AP proliferation supressor
    apx1=7 ; apx2=lx
    apy1=1 ; apy2=ly

    do i=ii1,ii2
      do j=jj1,jj2
        !epithelium
        k=cels(cell_grid_epi(i,j))%node(1)

        gex(k,1)=1d0
        gex(k,8)=1d0
        !borders
        if(i<=1.or.i>=lx.or.j<=1.or.j>=ly)then
          node(k)%hold=3
          node(node(k)%altre)%hold=3
        end if
        !placode
        if((i>px1.and.i<px2).and.(j>py1.and.j<py2))then
          gex(k,4)=1d0
          gex(k,13)=1d0
        else
          gex(k,5)=1d0
        end if
        !signalling centre
        if((i>sx1.and.i<sx2).and.(j>sy1.and.j<sy2))then
          gex(k,11)=1d0 ; gex(node(k)%altre,11)=1d0
          gex(k,22)=1d0 ; gex(node(k)%altre,22)=1d0
          gex(k,13)=0d0
        end if
        !AP proliferation supressor
        if((i>apx1.and.i<apx2).and.(j>apy1.and.j<apy2))then
          gex(k,23)=1d0
        end if


        !mesench
        do ii=1,layer
          k=cels(cell_grid_mes(i,j,ii))%node(1)
          gex(k,3)=1d0
          gex(k,10)=1d0
          !borders
          if(i<=1.or.i>=lx.or.j<=1.or.j>=ly)then
            node(k)%hold=3
          end if
          !placode
          if((i>px1.and.i<px2).and.(j>py1.and.j<py2))then
            gex(k,6)=1d0
            gex(k,15)=1d0
          else
            gex(k,7)=1d0
          end if
        end do
      end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!
    end if

  if(geometry==1)then  
    !ring=2  !scale
    ring=2 !tooth

    ring2=5

    polarized=0 !tooth
    !polarized=1 !scale

    do i=1,nd
      if(node(i)%tipus<3)then
        !if(node(i)%tipus==2)then
          !epithelial basal layer, basal side
          gex(i,1)=1d0  !receptors
          node(i)%marge=0
          d=sqrt(node(i)%x**2+node(i)%y**2)
          !setting expression on placode and interplacode
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,4)=1d0
            if(polarized==0.or.node(i)%x<=0.15)then
              gex(i,13)=1d0 ; gex(node(i)%altre,13)=1d0
            else
              gex(i,23)=1d0
              gex(i,13)=1d0 ; gex(node(i)%altre,13)=1d0
            end if
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,5)=1d0
          end if
          !!!!!!!!!!!!!
          !signalling centre
          if(d-2*node(i)%req*(ring-1)<epsilod)then
            node(i)%marge=0 ; node(node(i)%altre)%marge=0
            gex(i,11)=1d0 ; gex(node(i)%altre,11)=1d0
            gex(i,22)=1d0 ; gex(node(i)%altre,22)=1d0
            gex(i,13)=0d0 ; gex(node(i)%altre,13)=0d0
          end if

        !end if
        !epithelial basal layer, whole layer       
        gex(i,8)=1d0

      end if
      if(node(i)%tipus==3)then
        if(node(i)%z<=0)then
          !mesenchymal layer
          gex(i,3)=1d0
          gex(i,10)=1d0

          !setting expression on placode and interplacode
          d=sqrt(node(i)%x**2+node(i)%y**2)
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,6)=1d0
            gex(i,15)=1d0
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,7)=1d0
          end if
          !!!!!!!!!!!!!

        else
          !epithelial suprabasal layer
          gex(i,2)=1d0
          gex(i,9)=1d0

          !setting expression on placode and interplacode
          d=sqrt(node(i)%x**2+node(i)%y**2)
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,4)=1d0
            if(polarized==0.or.node(i)%x<=0.15)then
              gex(i,14)=1d0
            else
              gex(i,23)=1d0
              gex(i,14)=1d0
            end if
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,5)=1d0
          end if
          !!!!!!!!!!!!!
          !signalling centre
          !if(d-2*node(i)%req*(ring-1)<epsilod)then
          !  gex(i,11)=1d0
          !  gex(i,22)=1d0
          !  gex(i,14)=0d0 !; gex(node(i)%altre,13)=0d0
          !end if


        end if
      end if


    end do
    end if

    call update_npag

    node(:)%talone=0.0d0

    !ramax=maxval(node(:)%da)*3

    
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1



subroutine epidermal_organ_in_tooth_bud_s2

!real*8::pericel,radio2,req1,req2,adjustreq,di
!integer::kk2

integer:: lx,ly
integer:: x1,x2,y1,y2,ii1,ii2,jj1,jj2
integer:: x3,x4,y3,y4,ii3,ii4,jj3,jj4,jj5,y5
integer, allocatable::cell_grid_epi(:,:),cell_grid_mes(:,:,:)
!integer:: periodic_borders,nborders
integer::geometry
integer::xlayer
integer::ring,ring2,polarized
integer::px1,px2,py1,py2,sx1,sx2,sy1,sy2,apx1,apx2,apy1,apy2

!******* #1 DEFINING SPATIAL DIMENSIONS *******

  geometry=1 !geometry=1 means hexagonal epithelium ; geometry=2 means rectangular epithelium
  
  if (geometry==1)then
	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=10    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=10   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=0.75d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    do i=1,layer
      mesradicel(i)=mradicel
    end do
    
    xlayer=1 !this will make a layer of ECM beneath the layers of mesenchymal cells
    
  elseif(geometry==2)then
    layer=1      !number of planar cell layers
    lx=16
    ly=32
    
    allocate(cell_grid_epi(lx,ly),cell_grid_mes(lx,ly,layer))

  end if  
    !periodic_borders=1
    


    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(geometry==1)then
      if(radi>0.and.radicel>0)then
        j=0
        do i=1,radi-1
          j=j+i
        end do
        nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
        nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13
  
        j=0
        do i=1,radicel-1
          j=j+i
        end do
        ncelsepi=(6*j+1)
        ndepi=nodecel*ncelsepi
      else
        ncelsepi=0
        ndepi=0
        nodecel=0
        nodecela=0
      end if
  
      if(mradi>0.and.mradicel>0.and.layer>0)then
        if(packed==1)then !if mesenchyme is packed, we use mradi differently
          j=0
          do i=1,mradi-1
            j=j+i
          end do
          nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
          nodecela=2*nodecel+1
        else
          nodecel=mradi
          nodecela=2*nodecel+1
        end if
  
        ncelsmes=0
        do k=1,layer
          j=0 !;print*,"mesradicel",mesradicel(k)
          do i=1,mesradicel(k)-1
            j=j+i
          end do
          ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
        end do
        ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
        if(radi==0.and.radicel==0)then
          nodecel=mradi
          nodecela=2*nodecel+1
        else if(nodecel<mradi)then
          nodecel=radi
          nodecela=nodecel*2+1
        end if
      else
        ndmes=0
        ncelsmes=0
      end if
      ndx=(6*j+1)*xlayer
      nd=ndepi+ndmes+ndx
      !ncels=ncelsepi+ncelsmes
      ncels=ncelsepi+ncelsmes+ndx
   print*,"ncels",ncels,ncelsepi,ncelsmes,ndx
      nda=nd+10
      ncals=ncels+10
    elseif(geometry==2)then

      ncelsepi=lx*ly
      ndepi=ncelsepi*2
      ncelsmes=ncelsepi*layer
      ndmes=ncelsmes
      nd=ndepi+ndmes
      ncels=ncelsepi+ncelsmes
      nodecel=2 ; nodecela=5
      nda=nd+10
      ncals=ncels+10
    end if
    

  !End initializing dimension parameters

  print*,"inicial nd",nd,"ncels",ncels,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=1
    !End Allocatations


    

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=0.85d0
    deltamin=1d-2
    khold=1d0
    angletor=0.05



    k_bu=5d0
    ramax=0.35d0
    k_press=5d-1
    m_xwall=zmes !tooth
    mi_xwall=0d0 !tooth

    !m_xwall=0d0  !SCALE
    !mi_xwall=-5d-1 !SCALE

    ndmax=9d4
    
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=1
    ffu(2)=1 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !buoyancy

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(14)=1 !physical boundaries (walls)
    ffu(17)=1 !volume conservation (epithelial)
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise
    ffu(24)=0


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    !if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(mod(i,2)/=0)then  !basal
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
          !node(i)%reqp=-0.05d0
          !node(i)%req=0.20d0
        !node(i)%kplast=1d0
        else                      !apical
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !node(i)%kplast=1d1
        !  !node(i)%reqp=0.05d0
        !  !node(i)%req=0.30d0
        end if
        node(i)%req=0.25d0
        node(i)%reqcr=0.25d0 !; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.125
        node(i)%da=node(i)%req*2.0; 
        node(i)%ke=5d0
        node(i)%tor=3d0
        node(i)%stor=1d0
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=1d1
        node(i)%kvol=5d1
      end do
    !end if

	!Mesenchyme
    if(layer>0)then
      do i=ndepi+1,ndepi+ndmes
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !scale rect
        !node(i)%rep=0d0 ; node(i)%repcel=1d0
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*2.0
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

	!Suprabasal cells
    if(xlayer>0)then
      do i=ndepi+ndmes+1,nd
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*2.0
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if
    
    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)

    node(:)%hold=0


  if(geometry==1)then
    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    !if(radi>0.and.radicel>0)               call epiteli_sphere1(radi,radicel,radius,di,zepi)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if
!print*,"tall0"

    !Mesenchymal layer
    ii=0
    do i=ndepi+ndmes+1,nd
      ii=ii+1
      j=i-ndx
      node(i)%x=node(j)%x ; node(i)%y=node(j)%y ; node(i)%z=-0.5d0
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
      node(i)%icel=node(j)%icel+ndx ; node(i)%tipus=3
      !node(i)%icel=ncelsepi+ncelsmes+ii ; node(i)%tipus=3

      node(i)%repcel=5d-1
      !allocate(cels(node(i)%icel)%node(nodela))
      cels(node(i)%icel)%nunodes=1
      cels(node(i)%icel)%node(1)=i
      node(i)%marge=0
      !print*,"i",i,"j",j,"icel i",node(i)%icel,"icel j",node(j)%icel
      !node(i)%border=1
    end do
    ncelsmes=ncelsmes+ndx ; ndx=0
    ndmes=ncelsmes
!print*,"tall1"
      
    j=0
    do i=1,radicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    do i=j+1,ndepi
      node(i)%hold=2 ;node(i)%repcel=1d0
    end do
    j=0
    do i=1,mradicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    k=0
    do i=1,mradicel-1
      k=k+i
    end do
    k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
      node(i)%hold=2;node(i)%repcel=1d0
      !node(i)%border=1
      !node(i)%da=node(i)%req*2.0
      !node(i)%orix=0 ; node(i)%oriy=0
      !node(i)%rep=1d1;node(i)%repcel=1d1
    end do
    do i=ndepi+ndmes/2+j+1,nd
      node(i)%hold=2;node(i)%repcel=1d0
      !node(i)%border=1
    end do
    !do i=ndepi+k+1,nd
    !  node(i)%hold=2;node(i)%repcel=1d0
    !  node(i)%border=1
    !end do

    
  elseif(geometry==2)then
	do i=1,ncels
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do


    !epithelial cells
    !origin of rectangle
    aa=node(1)%req*2
    bb=sqrt(aa**2-(aa*0.5)**2)
    a=-real(lx)*0.5*aa ; b=-real(ly)*0.5*bb
    !print*,"a",a,"b",b,"aa",aa,"bb",bb,"lx",lx,"ly",ly
    ii=0 ; jj=0
    do i=1,ly
      do j=1,lx
        ii=ii+1 ;jj=jj+1 ; cell_grid_epi(j,i)=jj
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii+1 ;node(ii)%tipus=2
        ii=ii+1
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii-1 ;node(ii)%tipus=1
        cels(jj)%node(1)=ii-1 ; cels(jj)%node(2)=ii ;cels(jj)%nunodes=2; cels(jj)%ctipus=1
        if(i==1.or.i==ly.or.j==1.or.j==lx)then
          node(ii)%hold=2;node(ii-1)%hold=2 ;
          node(ii)%repcel=1d1;node(ii-1)%repcel=1d1;
          !node(ii)%border=1;node(ii-1)%border=1;
        end if
      end do
    end do

    !mesenchymals
    do k=1,layer
      !print*,"k",k
      do i=1,ly
        do j=1,lx
          ii=ii+1 ; jj=jj+1 ; cell_grid_mes(j,i,k)=jj
          if(mod(i,2)==0)then
            node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          else
            node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          end if
          node(ii)%icel=jj ; node(ii)%altre=0 ;node(ii)%tipus=3
          cels(jj)%node(1)=ii ;cels(jj)%nunodes=1 ; cels(jj)%ctipus=3
          if(i==1 .or. i==ly .or. j==1 .or. j==lx)then
            node(ii)%hold=2 ;node(ii)%repcel=1d1 ; !node(ii)%border=2 ;
          end if
          !if(i==1 .or. i==ly .or. j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
          !if(j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
        end do
      end do
    end do
  end if

  !m_xwall=maxval(node(1:nd)%x)+0.01
  !mi_xwall=minval(node(1:nd)%x)-0.01
  !m_ywall=maxval(node(1:nd)%y)+0.01
  !mi_ywall=minval(node(1:nd)%y)-0.01
  
!print*,"tall2"

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
  !node(:)%hold=2
  !node(:)%hold=0  

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    !if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    !end if
    !Mesenchymal
    if(layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=31
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=0d0 ; gen(1)%mu=5d-1 ;gen(1)%name="epithelial basal marker"
      gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=5d-1 ;gen(2)%name="epithelial suprabasal marker"
      gen(3)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=5d-1 ;gen(3)%name="mesenchymal marker"
      gen(4)%kindof=1 ; gen(4)%diffu=0d0 ; gen(4)%mu=5d-1 ;gen(4)%name="epi. placode marker"
      gen(5)%kindof=1 ; gen(5)%diffu=0d0 ; gen(5)%mu=5d-1 ;gen(5)%name="epi. interplacode marker"
      gen(6)%kindof=1 ; gen(6)%diffu=0d0 ; gen(6)%mu=5d-1 ;gen(6)%name="mes. placode marker"
      gen(7)%kindof=1 ; gen(7)%diffu=0d0 ; gen(7)%mu=5d-1 ;gen(7)%name="mes. interplacode marker"

      gen(8)%kindof=1 ; gen(8)%diffu=0d0 ; gen(8)%mu=5d-1 ;gen(8)%name="epi. basal adhesion molecule"
      gen(9)%kindof=1 ; gen(9)%diffu=0d0 ; gen(9)%mu=5d-1 ;gen(9)%name="epi. suprabasal adhesion molecule"
      gen(10)%kindof=1 ; gen(10)%diffu=0d0 ; gen(10)%mu=5d-1 ;gen(10)%name="mesench. adhesion molecule"

      gen(11)%kindof=1 ; gen(11)%diffu=0d0 ; gen(11)%mu=5d-1 ;gen(11)%name="epi. signalling centre marker"
      gen(12)%kindof=4 ; gen(12)%diffu=5d-1 ; gen(12)%mu=5d-1 ;gen(12)%name="epi. signal"

      gen(13)%kindof=2 ; gen(13)%diffu=0d0 ; gen(13)%mu=5d-1 ;gen(13)%name="epi. basal receptor inactive"
      gen(14)%kindof=2 ; gen(14)%diffu=0d0 ; gen(14)%mu=5d-1 ;gen(14)%name="epi. suprabasal receptor inactive"
      gen(15)%kindof=2 ; gen(15)%diffu=0d0 ; gen(15)%mu=5d-1 ;gen(15)%name="mesench. receptor inactive"

      gen(16)%kindof=8 ; gen(16)%diffu=0d0 ; gen(16)%mu=5d-1 ;gen(16)%name="epi. basal receptor active"
      gen(16)%npre=2 ; allocate(gen(16)%pre(gen(16)%npre)) ; gen(16)%pre(1)=12 ; gen(16)%pre(2)=13
      gen(16)%npost=2 ; allocate(gen(16)%post(gen(16)%npost)) ; gen(16)%post(1)=12 ; gen(16)%post(2)=13

      gen(17)%kindof=8 ; gen(17)%diffu=0d0 ; gen(17)%mu=5d-1 ;gen(17)%name="epi. suprabasal receptor active"
      gen(17)%npre=2 ; allocate(gen(17)%pre(gen(17)%npre)) ; gen(17)%pre(1)=12 ; gen(17)%pre(2)=14
      gen(17)%npost=2 ; allocate(gen(17)%post(gen(17)%npost)) ; gen(17)%post(1)=12 ; gen(17)%post(2)=14

      gen(18)%kindof=8 ; gen(18)%diffu=0d0 ; gen(18)%mu=5d-1 ;gen(18)%name="mesench. receptor active"
      gen(18)%npre=2 ; allocate(gen(18)%pre(gen(18)%npre)) ; gen(18)%pre(1)=12 ; gen(18)%pre(2)=15
      gen(18)%npost=2 ; allocate(gen(18)%post(gen(18)%npost)) ; gen(18)%post(1)=12 ; gen(18)%post(2)=15

      gen(19)%kindof=1 ; gen(19)%diffu=0d0 ; gen(19)%mu=5d-1 ;gen(19)%name="epi. basal effector"
      gen(20)%kindof=1 ; gen(20)%diffu=0d0 ; gen(20)%mu=5d-1 ;gen(20)%name="epi. suprabasal effector"
      gen(21)%kindof=1 ; gen(21)%diffu=0d0 ; gen(21)%mu=5d-1 ;gen(21)%name="mesenchymal effector"

      gen(22)%kindof=1 ; gen(22)%diffu=0d0 ; gen(22)%mu=5d-1 ;gen(22)%name="signalling centre adh. mol."
      
      gen(23)%kindof=1 ; gen(23)%diffu=0d0 ; gen(23)%mu=5d-1 ;gen(23)%name="s2 IEE adhesion molecule."
      gen(24)%kindof=1 ; gen(24)%diffu=0d0 ; gen(24)%mu=5d-1 ;gen(24)%name="s2 dermal papilla adhesion molecule"

      gen(25)%kindof=4 ; gen(25)%diffu=5d-1 ; gen(25)%mu=5d-1 ;gen(25)%name="s2 secondary signal"
      gen(26)%kindof=1 ; gen(26)%diffu=0d0 ; gen(26)%mu=5d-1 ;gen(26)%name="s2 enamel knot marker basal"
      gen(27)%kindof=1 ; gen(27)%diffu=0d0 ; gen(27)%mu=5d-1 ;gen(27)%name="s2 enamel knot marker supra"

      gen(28)%kindof=2 ; gen(28)%diffu=0d0 ; gen(28)%mu=5d-1 ;gen(28)%name="s2 epi. knot basal receptor inactive"
      gen(29)%kindof=2 ; gen(29)%diffu=0d0 ; gen(29)%mu=5d-1 ;gen(29)%name="s2 epi. knot supra receptor inactive"

      gen(30)%kindof=8 ; gen(30)%diffu=0d0 ; gen(30)%mu=5d-1 ;gen(30)%name="s2 epi. knot basal receptor active"
      gen(30)%npre=2 ; allocate(gen(30)%pre(gen(30)%npre)) ; gen(30)%pre(1)=25 ; gen(30)%pre(2)=28
      gen(30)%npost=2 ; allocate(gen(30)%post(gen(30)%npost)) ; gen(30)%post(1)=25 ; gen(30)%post(2)=28

      gen(31)%kindof=8 ; gen(31)%diffu=0d0 ; gen(31)%mu=5d-1 ;gen(31)%name="s2 epi. knot supra receptor active"
      gen(31)%npre=2 ; allocate(gen(31)%pre(gen(31)%npre)) ; gen(31)%pre(1)=25 ; gen(31)%pre(2)=29
      gen(31)%npost=2 ; allocate(gen(31)%post(gen(31)%npost)) ; gen(31)%post(1)=25 ; gen(31)%post(2)=29
    


 
    !Gene-behavior interactions
    
   
      gen(8)%wa(1)=1
      gen(9)%wa(1)=2
      gen(10)%wa(1)=3
      gen(22)%wa(1)=4
      gen(23)%wa(1)=5
      gen(24)%wa(1)=6

      !gen(19)%wa(1)=5
      !gen(20)%wa(1)=6
      !gen(21)%wa(1)=7
      
      gen(22)%wa(10)=1d0
      
      gen(11)%wa(nparam_per_node+2)=-5d5

      gen(19)%wa(nparam_per_node+2)=7d-2
      gen(20)%wa(nparam_per_node+2)=9d-2
      gen(21)%wa(nparam_per_node+2)=1.2d-1
      
      !growth parameters for inward bud
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=3.2d-1
      !gen(21)%wa(nparam_per_node+2)=5d-2
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !growth parameters for outward bud
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=5.0d-3!5.0d-2
      !gen(21)%wa(nparam_per_node+2)=3.2d-1
      !gen(23)%wa(nparam_per_node+2)=-5.0d0
      !!!!!!!!!!!!!!!!!!!!!!!!!!   
      !growth parameters for rectangle scale
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=0.0d0
      !gen(21)%wa(nparam_per_node+2)=1.0d-3
      !gen(23)%wa(nparam_per_node+2)=-7d1
      !gen(6)%wa(nparam_per_node+2)=1.0d-1
      !!!!!!!!!!!!!!!!!!!!!!!!!!   
      
    !Gene-gene interactions

      !structural genes
      !auto-maintenance
      gen(1)%w(1)=1d0
      gen(2)%w(2)=1d0
      gen(3)%w(3)=1d0
      gen(4)%w(4)=1d0
      gen(5)%w(5)=1d0
      gen(6)%w(6)=1d0
      gen(7)%w(7)=1d0
      gen(11)%w(11)=1d0
      !maintenance of stable gene expression
      gen(13)%w(4)=1d0
      gen(14)%w(4)=1d0
      gen(15)%w(6)=1d0
      gen(8)%w(1)=1d0
      gen(9)%w(2)=1d0
      gen(10)%w(3)=1d0
      gen(22)%w(11)=1d0
      gen(8)%w(11)=-1d5
      gen(9)%w(11)=-1d5
      !!!!!!!!!!!!!!!!
      gen(12)%w(11)=1d0
      gen(13)%w(11)=-1d5
      gen(14)%w(11)=-1d5
      gen(14)%w(1)=-1d5
      gen(13)%w(2)=-1d5

      !!!!!!!!!!!!!!!!! effector activation thresholds
      gen(19)%w(1)=-0.05d0
      gen(20)%w(2)=-0.05d0
      gen(21)%w(3)=-0.05d0
      !gen(19)%w(1)=0d0
      !gen(20)%w(2)=0d0
      !gen(21)%w(3)=0d0
      gen(19)%w(16)=1d2   !tooth
      gen(20)%w(17)=1d2
      gen(21)%w(18)=1d1
      !gen(19)%w(16)=1d2    !scale
      !gen(20)%w(17)=1d2
      !gen(21)%w(18)=1d1
      !gen(19)%w(16)=1d1   !tooth
      !gen(20)%w(17)=1d1
      !gen(21)%w(18)=1d1
      !gen(19)%w(16)=1d2    !rectangle scale
      !gen(20)%w(17)=1d2
      !gen(21)%w(18)=1d1


      !new interactions for s2    

      !gen(12)%w(26)=1d0 !s2 !the knot markers secrete the primary signal
      !gen(12)%w(27)=1d0 !s2
      !gen(12)%w(11)=0d0 !s2  !primary signal secretion is transfered to the knot markers
      !gen(25)%w(11)=1d0 !s2


      gen(23)%w(11)=1d0 !s2  !adhesion molecule activation
      gen(23)%w(26)=1d0 !s2
      gen(23)%w(19)=1d0 !s2
      gen(24)%w(21)=1d0 !s2

      gen(28)%w(4)=1d0 !s2 !activation of receptors
      gen(29)%w(4)=1d0 !s2
      gen(28)%w(2)=-1d5 !s2 !we restrict every receptor to a specific compartment
      gen(29)%w(1)=-1d5 !s2 !we restrict every receptor to a specific compartment


      gen(26)%w(1)=-0.05d0 !s2 !threshold activation for the knot markers
      gen(27)%w(2)=-0.05d0 !s2
      gen(26)%w(30)=1d2 !s2
      gen(27)%w(31)=1d2 !s2
      gen(19)%w(26)=-1d5 !s2
      gen(20)%w(27)=-1d5 !s2
 
 

      !signal-receptor interactions
      !epithelial basal binding
      gen(16)%nww=4
      gen(16)%ww(1,1)=12
      gen(16)%ww(1,2)=16
      gen(16)%ww(1,3)=1d0
      gen(16)%ww(2,1)=16
      gen(16)%ww(2,2)=12
      gen(16)%ww(2,3)=1d0
      gen(16)%ww(3,1)=13
      gen(16)%ww(3,2)=16
      gen(16)%ww(3,3)=1d0
      gen(16)%ww(4,1)=16
      gen(16)%ww(4,2)=13
      gen(16)%ww(4,3)=1d0
      !epithelial suprabasal receptor binding
      gen(17)%nww=4
      gen(17)%ww(1,1)=12
      gen(17)%ww(1,2)=17
      gen(17)%ww(1,3)=1d0
      gen(17)%ww(2,1)=17
      gen(17)%ww(2,2)=12
      gen(17)%ww(2,3)=1d0
      gen(17)%ww(3,1)=14
      gen(17)%ww(3,2)=15
      gen(17)%ww(3,3)=1d0
      gen(17)%ww(4,1)=17
      gen(17)%ww(4,2)=14
      gen(17)%ww(4,3)=1d0
      !mesenchymal receptor binding
      gen(18)%nww=4
      gen(18)%ww(1,1)=12
      gen(18)%ww(1,2)=18
      gen(18)%ww(1,3)=1d0
      gen(18)%ww(2,1)=18
      gen(18)%ww(2,2)=12
      gen(18)%ww(2,3)=1d0
      gen(18)%ww(3,1)=15
      gen(18)%ww(3,2)=18
      gen(18)%ww(3,3)=1d0
      gen(18)%ww(4,1)=18
      gen(18)%ww(4,2)=15
      gen(18)%ww(4,3)=1d0
     

      !knot receptors
      !epithelial basal binding
      gen(30)%nww=4
      gen(30)%ww(1,1)=25
      gen(30)%ww(1,2)=30
      gen(30)%ww(1,3)=1d0
      gen(30)%ww(2,1)=30
      gen(30)%ww(2,2)=25
      gen(30)%ww(2,3)=1d0
      gen(30)%ww(3,1)=28
      gen(30)%ww(3,2)=30
      gen(30)%ww(3,3)=1d0
      gen(30)%ww(4,1)=30
      gen(30)%ww(4,2)=28
      gen(30)%ww(4,3)=1d0
      !epithelial suprabasal receptor binding
      gen(31)%nww=4
      gen(31)%ww(1,1)=25
      gen(31)%ww(1,2)=31
      gen(31)%ww(1,3)=1d0
      gen(31)%ww(2,1)=31
      gen(31)%ww(2,2)=25
      gen(31)%ww(2,3)=1d0
      gen(31)%ww(3,1)=29
      gen(31)%ww(3,2)=31
      gen(31)%ww(3,3)=1d0
      gen(31)%ww(4,1)=31
      gen(31)%ww(4,2)=29
      gen(31)%ww(4,3)=1d0



    !Adhesion molecules

	ntipusadh=6
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      !as default, 1 and 2 should act as the same (basal and suprabasal)
      kadh(1,1)=1d0

      kadh(1,2)=1d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=1d0

      kadh(1,3)=5d-1 ; kadh(3,1)=kadh(1,3)
      !rectangle scale
      !kadh(1,3)=3d0 ; kadh(3,1)=kadh(1,3)
      kadh(2,3)=0d0 ; kadh(3,2)=kadh(2,3)
      kadh(3,3)=1d0

      kadh(1,4)=5d-1 ; kadh(4,1)=kadh(1,4)
      kadh(2,4)=5d-1 ; kadh(4,2)=kadh(2,4)
      kadh(3,4)=5d-1 ; kadh(4,3)=kadh(3,4)
      kadh(4,4)=1d1 !tooth
      !kadh(4,4)=5d0  !scale

      kadh(1,5)=0d0 ; kadh(5,1)=kadh(1,5)
      kadh(2,5)=0d0 ; kadh(5,2)=kadh(2,5)
      kadh(3,5)=0d0 ; kadh(5,3)=kadh(3,5)
      kadh(4,5)=0d0 ; kadh(5,4)=kadh(4,5)
      kadh(5,5)=0d0

      kadh(1,6)=0d0 ; kadh(6,1)=kadh(1,6)
      kadh(2,6)=0d0 ; kadh(6,2)=kadh(2,6)
      kadh(3,6)=0d0 ; kadh(6,3)=kadh(3,6)
      kadh(4,6)=0d0 ; kadh(6,4)=kadh(4,6)
      kadh(5,6)=1d1 ; kadh(6,5)=kadh(5,6)
      kadh(6,6)=1d1

      !kadh(1,7)=5d-1 ; kadh(7,1)=kadh(1,7)
      !kadh(2,7)=0d0 ; kadh(7,2)=kadh(2,7)
      !kadh(3,7)=1d0 ; kadh(7,3)=kadh(3,7)
      !kadh(4,7)=5d-1 ; kadh(7,4)=kadh(4,7)
      !kadh(5,7)=5d-1 ; kadh(7,5)=kadh(5,7)
      !kadh(6,7)=0d0 ; kadh(7,6)=kadh(6,7)
      !kadh(7,7)=1d0

    end if

    !Gene expression on nodes

    
    !if(geometry==1)then
    !  gex(1,2)=1d0
    !  gex(ndepi+1,2)=1d0
    !  !gex(ndepi+1,2)=1d0
    !end if
    
    if(geometry==2)then
    !draw rectangles
    !!!!!!!!!!!!!!!!!!!!!!!!!
    jj1=1 ; jj2=ly
    ii1=1 ; ii2=lx
    !placode range
    px1=3 ; px2=lx-(px1-1)
    py1=3 ; py2=ly-(py1-1)
    !signalling centre
    sx1=7 ; sx2=lx-(sx1-1)
    sy1=5 ; sy2=ly-(sy1-1)
    !AP proliferation supressor
    apx1=7 ; apx2=lx
    apy1=1 ; apy2=ly

    do i=ii1,ii2
      do j=jj1,jj2
        !epithelium
        k=cels(cell_grid_epi(i,j))%node(1)

        gex(k,1)=1d0
        gex(k,8)=1d0
        !borders
        if(i<=1.or.i>=lx.or.j<=1.or.j>=ly)then
          node(k)%hold=3
          node(node(k)%altre)%hold=3
        end if
        !placode
        if((i>px1.and.i<px2).and.(j>py1.and.j<py2))then
          gex(k,4)=1d0
          gex(k,13)=1d0
        else
          gex(k,5)=1d0
        end if
        !signalling centre
        if((i>sx1.and.i<sx2).and.(j>sy1.and.j<sy2))then
          gex(k,11)=1d0 ; gex(node(k)%altre,11)=1d0
          gex(k,22)=1d0 ; gex(node(k)%altre,22)=1d0
          gex(k,13)=0d0
        end if
        !AP proliferation supressor
        if((i>apx1.and.i<apx2).and.(j>apy1.and.j<apy2))then
          gex(k,23)=1d0
        end if


        !mesench
        do ii=1,layer
          k=cels(cell_grid_mes(i,j,ii))%node(1)
          gex(k,3)=1d0
          gex(k,10)=1d0
          !borders
          if(i<=1.or.i>=lx.or.j<=1.or.j>=ly)then
            node(k)%hold=3
          end if
          !placode
          if((i>px1.and.i<px2).and.(j>py1.and.j<py2))then
            gex(k,6)=1d0
            gex(k,15)=1d0
          else
            gex(k,7)=1d0
          end if
        end do
      end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!
    end if

  if(geometry==1)then  
    !ring=2  !scale
    ring=2 !tooth

    ring2=7

    polarized=0 !tooth
    !polarized=1 !scale

    do i=1,nd
      if(node(i)%tipus<3)then
        !if(node(i)%tipus==2)then
          !epithelial basal layer, basal side
          gex(i,1)=1d0  !receptors
          node(i)%marge=0
          d=sqrt(node(i)%x**2+node(i)%y**2)
          !setting expression on placode and interplacode
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,4)=1d0
            if(polarized==0.or.node(i)%x<=0.15)then
              gex(i,13)=1d0 !; gex(node(i)%altre,13)=1d0
            else
              !gex(i,23)=1d0
              gex(i,13)=1d0 !; gex(node(i)%altre,13)=1d0
            end if
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,5)=1d0
          end if
          !!!!!!!!!!!!!
          !signalling centre
          if(d-2*node(i)%req*(ring-1)<epsilod)then
            node(i)%marge=0 !; node(node(i)%altre)%marge=0
            gex(i,11)=1d0 !; gex(node(i)%altre,11)=1d0
            gex(i,22)=1d0 !; gex(node(i)%altre,22)=1d0
            gex(i,13)=0d0 !; gex(node(i)%altre,13)=0d0
          end if

        !end if
        !epithelial basal layer, whole layer       
        gex(i,8)=1d0

      end if
      if(node(i)%tipus==3)then
        if(node(i)%z<=0)then
          !mesenchymal layer
          gex(i,3)=1d0
          gex(i,10)=1d0

          !setting expression on placode and interplacode
          d=sqrt(node(i)%x**2+node(i)%y**2)
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,6)=1d0
            gex(i,15)=1d0
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,7)=1d0
          end if
          !!!!!!!!!!!!!

        else
          !epithelial suprabasal layer
          gex(i,2)=1d0
          gex(i,9)=1d0

          !setting expression on placode and interplacode
          d=sqrt(node(i)%x**2+node(i)%y**2)
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,4)=1d0
            if(polarized==0.or.node(i)%x<=0.15)then
              gex(i,14)=1d0
            else
              gex(i,23)=1d0
              gex(i,14)=1d0
            end if
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,5)=1d0
          end if
          !!!!!!!!!!!!!
          !signalling centre
          !if(d-2*node(i)%req*(ring-1)<epsilod)then
          !  gex(i,11)=1d0
          !  gex(i,22)=1d0
          !  gex(i,14)=0d0 !; gex(node(i)%altre,13)=0d0
          !end if


        end if
      end if


    end do
    end if

    call update_npag

    node(:)%talone=0.0d0

    !ramax=maxval(node(:)%da)*3

    
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine epidermal_organ_in_tooth_bud_s2_small

!real*8::pericel,radio2,req1,req2,adjustreq,di
!integer::kk2

integer:: lx,ly
integer:: x1,x2,y1,y2,ii1,ii2,jj1,jj2
integer:: x3,x4,y3,y4,ii3,ii4,jj3,jj4,jj5,y5
integer, allocatable::cell_grid_epi(:,:),cell_grid_mes(:,:,:)
!integer:: periodic_borders,nborders
integer::geometry
integer::xlayer
integer::ring,ring2,polarized
integer::px1,px2,py1,py2,sx1,sx2,sy1,sy2,apx1,apx2,apy1,apy2

!******* #1 DEFINING SPATIAL DIMENSIONS *******

  geometry=1 !geometry=1 means hexagonal epithelium ; geometry=2 means rectangular epithelium
  
  if (geometry==1)then
	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=5    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=5   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=0.75d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    do i=1,layer
      mesradicel(i)=mradicel
    end do
    
    xlayer=1 !this will make a layer of ECM beneath the layers of mesenchymal cells
    
  elseif(geometry==2)then
    layer=1      !number of planar cell layers
    lx=16
    ly=32
    
    allocate(cell_grid_epi(lx,ly),cell_grid_mes(lx,ly,layer))

  end if  
    !periodic_borders=1
    


    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(geometry==1)then
      if(radi>0.and.radicel>0)then
        j=0
        do i=1,radi-1
          j=j+i
        end do
        nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
        nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13
  
        j=0
        do i=1,radicel-1
          j=j+i
        end do
        ncelsepi=(6*j+1)
        ndepi=nodecel*ncelsepi
      else
        ncelsepi=0
        ndepi=0
        nodecel=0
        nodecela=0
      end if
  
      if(mradi>0.and.mradicel>0.and.layer>0)then
        if(packed==1)then !if mesenchyme is packed, we use mradi differently
          j=0
          do i=1,mradi-1
            j=j+i
          end do
          nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
          nodecela=2*nodecel+1
        else
          nodecel=mradi
          nodecela=2*nodecel+1
        end if
  
        ncelsmes=0
        do k=1,layer
          j=0 !;print*,"mesradicel",mesradicel(k)
          do i=1,mesradicel(k)-1
            j=j+i
          end do
          ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
        end do
        ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
        if(radi==0.and.radicel==0)then
          nodecel=mradi
          nodecela=2*nodecel+1
        else if(nodecel<mradi)then
          nodecel=radi
          nodecela=nodecel*2+1
        end if
      else
        ndmes=0
        ncelsmes=0
      end if
      ndx=(6*j+1)*xlayer
      nd=ndepi+ndmes+ndx
      !ncels=ncelsepi+ncelsmes
      ncels=ncelsepi+ncelsmes+ndx
   print*,"ncels",ncels,ncelsepi,ncelsmes,ndx
      nda=nd+10
      ncals=ncels+10
    elseif(geometry==2)then

      ncelsepi=lx*ly
      ndepi=ncelsepi*2
      ncelsmes=ncelsepi*layer
      ndmes=ncelsmes
      nd=ndepi+ndmes
      ncels=ncelsepi+ncelsmes
      nodecel=2 ; nodecela=5
      nda=nd+10
      ncals=ncels+10
    end if
    

  !End initializing dimension parameters

  print*,"inicial nd",nd,"ncels",ncels,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=1
    !End Allocatations


    

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=0.85d0
    deltamin=1d-2
    khold=1d0
    angletor=0.05



    k_bu=5d0
    ramax=0.35d0
    k_press=5d-1
    m_xwall=zmes !tooth
    mi_xwall=0d0 !tooth

    !m_xwall=0d0  !SCALE
    !mi_xwall=-5d-1 !SCALE

    ndmax=9d4
    
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=1
    ffu(2)=1 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !buoyancy

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(14)=1 !physical boundaries (walls)
    ffu(17)=1 !volume conservation (epithelial)
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise
    ffu(24)=0


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    !if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(mod(i,2)/=0)then  !basal
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
          !node(i)%reqp=-0.05d0
          !node(i)%req=0.20d0
        !node(i)%kplast=1d0
        else                      !apical
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !node(i)%kplast=1d1
        !  !node(i)%reqp=0.05d0
        !  !node(i)%req=0.30d0
        end if
        node(i)%req=0.25d0
        node(i)%reqcr=0.25d0 !; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.125
        node(i)%da=node(i)%req*2.0; 
        node(i)%ke=5d0
        node(i)%tor=3d0
        node(i)%stor=1d0
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=1d1
        node(i)%kvol=5d1
      end do
    !end if

	!Mesenchyme
    if(layer>0)then
      do i=ndepi+1,ndepi+ndmes
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !scale rect
        !node(i)%rep=0d0 ; node(i)%repcel=1d0
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*2.0
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

	!Suprabasal cells
    if(xlayer>0)then
      do i=ndepi+ndmes+1,nd
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*2.0
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if
    
    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)

    node(:)%hold=0


  if(geometry==1)then
    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    !if(radi>0.and.radicel>0)               call epiteli_sphere1(radi,radicel,radius,di,zepi)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if
!print*,"tall0"

    !Mesenchymal layer
    ii=0
    do i=ndepi+ndmes+1,nd
      ii=ii+1
      j=i-ndx
      node(i)%x=node(j)%x ; node(i)%y=node(j)%y ; node(i)%z=-0.5d0
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
      node(i)%icel=node(j)%icel+ndx ; node(i)%tipus=3
      !node(i)%icel=ncelsepi+ncelsmes+ii ; node(i)%tipus=3

      node(i)%repcel=5d-1
      !allocate(cels(node(i)%icel)%node(nodela))
      cels(node(i)%icel)%nunodes=1
      cels(node(i)%icel)%node(1)=i
      node(i)%marge=0
      !print*,"i",i,"j",j,"icel i",node(i)%icel,"icel j",node(j)%icel
      !node(i)%border=1
    end do
    ncelsmes=ncelsmes+ndx ; ndx=0
    ndmes=ncelsmes
!print*,"tall1"
      
    j=0
    do i=1,radicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    do i=j+1,ndepi
      node(i)%hold=2 ;node(i)%repcel=1d0
    end do
    j=0
    do i=1,mradicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    k=0
    do i=1,mradicel-1
      k=k+i
    end do
    k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
      node(i)%hold=2;node(i)%repcel=1d0
      !node(i)%border=1
      !node(i)%da=node(i)%req*2.0
      !node(i)%orix=0 ; node(i)%oriy=0
      !node(i)%rep=1d1;node(i)%repcel=1d1
    end do
    do i=ndepi+ndmes/2+j+1,nd
      node(i)%hold=2;node(i)%repcel=1d0
      !node(i)%border=1
    end do
    !do i=ndepi+k+1,nd
    !  node(i)%hold=2;node(i)%repcel=1d0
    !  node(i)%border=1
    !end do

    
  elseif(geometry==2)then
	do i=1,ncels
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do


    !epithelial cells
    !origin of rectangle
    aa=node(1)%req*2
    bb=sqrt(aa**2-(aa*0.5)**2)
    a=-real(lx)*0.5*aa ; b=-real(ly)*0.5*bb
    !print*,"a",a,"b",b,"aa",aa,"bb",bb,"lx",lx,"ly",ly
    ii=0 ; jj=0
    do i=1,ly
      do j=1,lx
        ii=ii+1 ;jj=jj+1 ; cell_grid_epi(j,i)=jj
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii+1 ;node(ii)%tipus=2
        ii=ii+1
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii-1 ;node(ii)%tipus=1
        cels(jj)%node(1)=ii-1 ; cels(jj)%node(2)=ii ;cels(jj)%nunodes=2; cels(jj)%ctipus=1
        if(i==1.or.i==ly.or.j==1.or.j==lx)then
          node(ii)%hold=2;node(ii-1)%hold=2 ;
          node(ii)%repcel=1d1;node(ii-1)%repcel=1d1;
          !node(ii)%border=1;node(ii-1)%border=1;
        end if
      end do
    end do

    !mesenchymals
    do k=1,layer
      !print*,"k",k
      do i=1,ly
        do j=1,lx
          ii=ii+1 ; jj=jj+1 ; cell_grid_mes(j,i,k)=jj
          if(mod(i,2)==0)then
            node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          else
            node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          end if
          node(ii)%icel=jj ; node(ii)%altre=0 ;node(ii)%tipus=3
          cels(jj)%node(1)=ii ;cels(jj)%nunodes=1 ; cels(jj)%ctipus=3
          if(i==1 .or. i==ly .or. j==1 .or. j==lx)then
            node(ii)%hold=2 ;node(ii)%repcel=1d1 ; !node(ii)%border=2 ;
          end if
          !if(i==1 .or. i==ly .or. j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
          !if(j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
        end do
      end do
    end do
  end if

  !m_xwall=maxval(node(1:nd)%x)+0.01
  !mi_xwall=minval(node(1:nd)%x)-0.01
  !m_ywall=maxval(node(1:nd)%y)+0.01
  !mi_ywall=minval(node(1:nd)%y)-0.01
  
!print*,"tall2"

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
  !node(:)%hold=2
  !node(:)%hold=0  

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    !if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    !end if
    !Mesenchymal
    if(layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=31
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=0d0 ; gen(1)%mu=5d-1 ;gen(1)%name="epithelial basal marker"
      gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=5d-1 ;gen(2)%name="epithelial suprabasal marker"
      gen(3)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=5d-1 ;gen(3)%name="mesenchymal marker"
      gen(4)%kindof=1 ; gen(4)%diffu=0d0 ; gen(4)%mu=5d-1 ;gen(4)%name="epi. placode marker"
      gen(5)%kindof=1 ; gen(5)%diffu=0d0 ; gen(5)%mu=5d-1 ;gen(5)%name="epi. interplacode marker"
      gen(6)%kindof=1 ; gen(6)%diffu=0d0 ; gen(6)%mu=5d-1 ;gen(6)%name="mes. placode marker"
      gen(7)%kindof=1 ; gen(7)%diffu=0d0 ; gen(7)%mu=5d-1 ;gen(7)%name="mes. interplacode marker"

      gen(8)%kindof=1 ; gen(8)%diffu=0d0 ; gen(8)%mu=5d-1 ;gen(8)%name="epi. basal adhesion molecule"
      gen(9)%kindof=1 ; gen(9)%diffu=0d0 ; gen(9)%mu=5d-1 ;gen(9)%name="epi. suprabasal adhesion molecule"
      gen(10)%kindof=1 ; gen(10)%diffu=0d0 ; gen(10)%mu=5d-1 ;gen(10)%name="mesench. adhesion molecule"

      gen(11)%kindof=1 ; gen(11)%diffu=0d0 ; gen(11)%mu=5d-1 ;gen(11)%name="epi. signalling centre marker"
      gen(12)%kindof=4 ; gen(12)%diffu=5d-1 ; gen(12)%mu=5d-1 ;gen(12)%name="epi. signal"

      gen(13)%kindof=2 ; gen(13)%diffu=0d0 ; gen(13)%mu=5d-1 ;gen(13)%name="epi. basal receptor inactive"
      gen(14)%kindof=2 ; gen(14)%diffu=0d0 ; gen(14)%mu=5d-1 ;gen(14)%name="epi. suprabasal receptor inactive"
      gen(15)%kindof=2 ; gen(15)%diffu=0d0 ; gen(15)%mu=5d-1 ;gen(15)%name="mesench. receptor inactive"

      gen(16)%kindof=8 ; gen(16)%diffu=0d0 ; gen(16)%mu=5d-1 ;gen(16)%name="epi. basal receptor active"
      gen(16)%npre=2 ; allocate(gen(16)%pre(gen(16)%npre)) ; gen(16)%pre(1)=12 ; gen(16)%pre(2)=13
      gen(16)%npost=2 ; allocate(gen(16)%post(gen(16)%npost)) ; gen(16)%post(1)=12 ; gen(16)%post(2)=13

      gen(17)%kindof=8 ; gen(17)%diffu=0d0 ; gen(17)%mu=5d-1 ;gen(17)%name="epi. suprabasal receptor active"
      gen(17)%npre=2 ; allocate(gen(17)%pre(gen(17)%npre)) ; gen(17)%pre(1)=12 ; gen(17)%pre(2)=14
      gen(17)%npost=2 ; allocate(gen(17)%post(gen(17)%npost)) ; gen(17)%post(1)=12 ; gen(17)%post(2)=14

      gen(18)%kindof=8 ; gen(18)%diffu=0d0 ; gen(18)%mu=5d-1 ;gen(18)%name="mesench. receptor active"
      gen(18)%npre=2 ; allocate(gen(18)%pre(gen(18)%npre)) ; gen(18)%pre(1)=12 ; gen(18)%pre(2)=15
      gen(18)%npost=2 ; allocate(gen(18)%post(gen(18)%npost)) ; gen(18)%post(1)=12 ; gen(18)%post(2)=15

      gen(19)%kindof=1 ; gen(19)%diffu=0d0 ; gen(19)%mu=5d-1 ;gen(19)%name="epi. basal effector"
      gen(20)%kindof=1 ; gen(20)%diffu=0d0 ; gen(20)%mu=5d-1 ;gen(20)%name="epi. suprabasal effector"
      gen(21)%kindof=1 ; gen(21)%diffu=0d0 ; gen(21)%mu=5d-1 ;gen(21)%name="mesenchymal effector"

      gen(22)%kindof=1 ; gen(22)%diffu=0d0 ; gen(22)%mu=5d-1 ;gen(22)%name="signalling centre adh. mol."
      
      gen(23)%kindof=1 ; gen(23)%diffu=0d0 ; gen(23)%mu=5d-1 ;gen(23)%name="s2 IEE adhesion molecule."
      gen(24)%kindof=1 ; gen(24)%diffu=0d0 ; gen(24)%mu=5d-1 ;gen(24)%name="s2 dermal papilla adhesion molecule"

      gen(25)%kindof=4 ; gen(25)%diffu=5d-1 ; gen(25)%mu=5d-1 ;gen(25)%name="s2 secondary signal"
      gen(26)%kindof=1 ; gen(26)%diffu=0d0 ; gen(26)%mu=5d-1 ;gen(26)%name="s2 enamel knot marker basal"
      gen(27)%kindof=1 ; gen(27)%diffu=0d0 ; gen(27)%mu=5d-1 ;gen(27)%name="s2 enamel knot marker supra"

      gen(28)%kindof=2 ; gen(28)%diffu=0d0 ; gen(28)%mu=5d-1 ;gen(28)%name="s2 epi. knot basal receptor inactive"
      gen(29)%kindof=2 ; gen(29)%diffu=0d0 ; gen(29)%mu=5d-1 ;gen(29)%name="s2 epi. knot supra receptor inactive"

      gen(30)%kindof=8 ; gen(30)%diffu=0d0 ; gen(30)%mu=5d-1 ;gen(30)%name="s2 epi. knot basal receptor active"
      gen(30)%npre=2 ; allocate(gen(30)%pre(gen(30)%npre)) ; gen(30)%pre(1)=25 ; gen(30)%pre(2)=28
      gen(30)%npost=2 ; allocate(gen(30)%post(gen(30)%npost)) ; gen(30)%post(1)=25 ; gen(30)%post(2)=28

      gen(31)%kindof=8 ; gen(31)%diffu=0d0 ; gen(31)%mu=5d-1 ;gen(31)%name="s2 epi. knot supra receptor active"
      gen(31)%npre=2 ; allocate(gen(31)%pre(gen(31)%npre)) ; gen(31)%pre(1)=25 ; gen(31)%pre(2)=29
      gen(31)%npost=2 ; allocate(gen(31)%post(gen(31)%npost)) ; gen(31)%post(1)=25 ; gen(31)%post(2)=29
    


 
    !Gene-behavior interactions
    
   
      gen(8)%wa(1)=1
      gen(9)%wa(1)=2
      gen(10)%wa(1)=3
      gen(22)%wa(1)=4
      gen(23)%wa(1)=5
      gen(24)%wa(1)=6

      !gen(19)%wa(1)=5
      !gen(20)%wa(1)=6
      !gen(21)%wa(1)=7
      
      gen(22)%wa(10)=1d0
      
      gen(11)%wa(nparam_per_node+2)=-5d5

      gen(19)%wa(nparam_per_node+2)=1.2d-1
      gen(20)%wa(nparam_per_node+2)=3.2d-1
      gen(21)%wa(nparam_per_node+2)=1.2d-1
      
      !growth parameters for inward bud
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=3.2d-1
      !gen(21)%wa(nparam_per_node+2)=5d-2
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !growth parameters for outward bud
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=5.0d-3!5.0d-2
      !gen(21)%wa(nparam_per_node+2)=3.2d-1
      !gen(23)%wa(nparam_per_node+2)=-5.0d0
      !!!!!!!!!!!!!!!!!!!!!!!!!!   
      !growth parameters for rectangle scale
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=0.0d0
      !gen(21)%wa(nparam_per_node+2)=1.0d-3
      !gen(23)%wa(nparam_per_node+2)=-7d1
      !gen(6)%wa(nparam_per_node+2)=1.0d-1
      !!!!!!!!!!!!!!!!!!!!!!!!!!   
      
    !Gene-gene interactions

      !structural genes
      !auto-maintenance
      gen(1)%w(1)=1d0
      gen(2)%w(2)=1d0
      gen(3)%w(3)=1d0
      gen(4)%w(4)=1d0
      gen(5)%w(5)=1d0
      gen(6)%w(6)=1d0
      gen(7)%w(7)=1d0
      gen(11)%w(11)=1d0
      !maintenance of stable gene expression
      gen(13)%w(4)=1d0
      gen(14)%w(4)=1d0
      gen(15)%w(6)=1d0
      gen(8)%w(1)=1d0
      gen(9)%w(2)=1d0
      gen(10)%w(3)=1d0
      gen(22)%w(11)=1d0
      gen(8)%w(11)=-1d5
      gen(9)%w(11)=-1d5
      !!!!!!!!!!!!!!!!
      gen(12)%w(11)=1d0
      gen(13)%w(11)=-1d5
      gen(14)%w(11)=-1d5
      gen(14)%w(1)=-1d5
      gen(13)%w(2)=-1d5

      !!!!!!!!!!!!!!!!! effector activation thresholds
      gen(19)%w(1)=-0.05d0
      gen(20)%w(2)=-0.05d0
      gen(21)%w(3)=-0.05d0
      !gen(19)%w(1)=0d0
      !gen(20)%w(2)=0d0
      !gen(21)%w(3)=0d0
      gen(19)%w(16)=5d1   !tooth
      gen(20)%w(17)=5d1
      gen(21)%w(18)=1d1
      !gen(19)%w(16)=1d2    !scale
      !gen(20)%w(17)=1d2
      !gen(21)%w(18)=1d1
      !gen(19)%w(16)=1d1   !tooth
      !gen(20)%w(17)=1d1
      !gen(21)%w(18)=1d1
      !gen(19)%w(16)=1d2    !rectangle scale
      !gen(20)%w(17)=1d2
      !gen(21)%w(18)=1d1


      !new interactions for s2    

      !gen(12)%w(26)=1d0 !s2 !the knot markers secrete the primary signal
      !gen(12)%w(27)=1d0 !s2
      !gen(12)%w(11)=0d0 !s2  !primary signal secretion is transfered to the knot markers
      !gen(25)%w(11)=1d0 !s2


      gen(23)%w(11)=1d0 !s2  !adhesion molecule activation
      gen(23)%w(26)=1d0 !s2
      gen(23)%w(19)=1d0 !s2
      gen(24)%w(21)=1d0 !s2

      gen(28)%w(4)=1d0 !s2 !activation of receptors
      gen(29)%w(4)=1d0 !s2
      gen(28)%w(2)=-1d5 !s2 !we restrict every receptor to a specific compartment
      gen(29)%w(1)=-1d5 !s2 !we restrict every receptor to a specific compartment


      gen(26)%w(1)=-0.05d0 !s2 !threshold activation for the knot markers
      gen(27)%w(2)=-0.05d0 !s2
      gen(26)%w(30)=2d1 !s2
      gen(27)%w(31)=2d1 !s2
      gen(19)%w(26)=-1d5 !s2
      gen(20)%w(27)=-1d5 !s2
 
 

      !signal-receptor interactions
      !epithelial basal binding
      gen(16)%nww=4
      gen(16)%ww(1,1)=12
      gen(16)%ww(1,2)=16
      gen(16)%ww(1,3)=1d0
      gen(16)%ww(2,1)=16
      gen(16)%ww(2,2)=12
      gen(16)%ww(2,3)=1d0
      gen(16)%ww(3,1)=13
      gen(16)%ww(3,2)=16
      gen(16)%ww(3,3)=1d0
      gen(16)%ww(4,1)=16
      gen(16)%ww(4,2)=13
      gen(16)%ww(4,3)=1d0
      !epithelial suprabasal receptor binding
      gen(17)%nww=4
      gen(17)%ww(1,1)=12
      gen(17)%ww(1,2)=17
      gen(17)%ww(1,3)=1d0
      gen(17)%ww(2,1)=17
      gen(17)%ww(2,2)=12
      gen(17)%ww(2,3)=1d0
      gen(17)%ww(3,1)=14
      gen(17)%ww(3,2)=15
      gen(17)%ww(3,3)=1d0
      gen(17)%ww(4,1)=17
      gen(17)%ww(4,2)=14
      gen(17)%ww(4,3)=1d0
      !mesenchymal receptor binding
      gen(18)%nww=4
      gen(18)%ww(1,1)=12
      gen(18)%ww(1,2)=18
      gen(18)%ww(1,3)=1d0
      gen(18)%ww(2,1)=18
      gen(18)%ww(2,2)=12
      gen(18)%ww(2,3)=1d0
      gen(18)%ww(3,1)=15
      gen(18)%ww(3,2)=18
      gen(18)%ww(3,3)=1d0
      gen(18)%ww(4,1)=18
      gen(18)%ww(4,2)=15
      gen(18)%ww(4,3)=1d0
     

      !knot receptors
      !epithelial basal binding
      gen(30)%nww=4
      gen(30)%ww(1,1)=25
      gen(30)%ww(1,2)=30
      gen(30)%ww(1,3)=1d0
      gen(30)%ww(2,1)=30
      gen(30)%ww(2,2)=25
      gen(30)%ww(2,3)=1d0
      gen(30)%ww(3,1)=28
      gen(30)%ww(3,2)=30
      gen(30)%ww(3,3)=1d0
      gen(30)%ww(4,1)=30
      gen(30)%ww(4,2)=28
      gen(30)%ww(4,3)=1d0
      !epithelial suprabasal receptor binding
      gen(31)%nww=4
      gen(31)%ww(1,1)=25
      gen(31)%ww(1,2)=31
      gen(31)%ww(1,3)=1d0
      gen(31)%ww(2,1)=31
      gen(31)%ww(2,2)=25
      gen(31)%ww(2,3)=1d0
      gen(31)%ww(3,1)=29
      gen(31)%ww(3,2)=31
      gen(31)%ww(3,3)=1d0
      gen(31)%ww(4,1)=31
      gen(31)%ww(4,2)=29
      gen(31)%ww(4,3)=1d0



    !Adhesion molecules

	ntipusadh=6
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      !as default, 1 and 2 should act as the same (basal and suprabasal)
      kadh(1,1)=1d0

      kadh(1,2)=1d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=1d0

      kadh(1,3)=5d-1 ; kadh(3,1)=kadh(1,3)
      !rectangle scale
      !kadh(1,3)=3d0 ; kadh(3,1)=kadh(1,3)
      kadh(2,3)=0d0 ; kadh(3,2)=kadh(2,3)
      kadh(3,3)=1d0

      kadh(1,4)=5d-1 ; kadh(4,1)=kadh(1,4)
      kadh(2,4)=5d-1 ; kadh(4,2)=kadh(2,4)
      kadh(3,4)=5d-1 ; kadh(4,3)=kadh(3,4)
      kadh(4,4)=1d1 !tooth
      !kadh(4,4)=5d0  !scale

      kadh(1,5)=0d0 ; kadh(5,1)=kadh(1,5)
      kadh(2,5)=0d0 ; kadh(5,2)=kadh(2,5)
      kadh(3,5)=0d0 ; kadh(5,3)=kadh(3,5)
      kadh(4,5)=0d0 ; kadh(5,4)=kadh(4,5)
      kadh(5,5)=0d0

      kadh(1,6)=0d0 ; kadh(6,1)=kadh(1,6)
      kadh(2,6)=0d0 ; kadh(6,2)=kadh(2,6)
      kadh(3,6)=0d0 ; kadh(6,3)=kadh(3,6)
      kadh(4,6)=0d0 ; kadh(6,4)=kadh(4,6)
      kadh(5,6)=1d1 ; kadh(6,5)=kadh(5,6)
      kadh(6,6)=1d1

      !kadh(1,7)=5d-1 ; kadh(7,1)=kadh(1,7)
      !kadh(2,7)=0d0 ; kadh(7,2)=kadh(2,7)
      !kadh(3,7)=1d0 ; kadh(7,3)=kadh(3,7)
      !kadh(4,7)=5d-1 ; kadh(7,4)=kadh(4,7)
      !kadh(5,7)=5d-1 ; kadh(7,5)=kadh(5,7)
      !kadh(6,7)=0d0 ; kadh(7,6)=kadh(6,7)
      !kadh(7,7)=1d0

    end if

    !Gene expression on nodes

    
    !if(geometry==1)then
    !  gex(1,2)=1d0
    !  gex(ndepi+1,2)=1d0
    !  !gex(ndepi+1,2)=1d0
    !end if
    
    if(geometry==2)then
    !draw rectangles
    !!!!!!!!!!!!!!!!!!!!!!!!!
    jj1=1 ; jj2=ly
    ii1=1 ; ii2=lx
    !placode range
    px1=3 ; px2=lx-(px1-1)
    py1=3 ; py2=ly-(py1-1)
    !signalling centre
    sx1=7 ; sx2=lx-(sx1-1)
    sy1=5 ; sy2=ly-(sy1-1)
    !AP proliferation supressor
    apx1=7 ; apx2=lx
    apy1=1 ; apy2=ly

    do i=ii1,ii2
      do j=jj1,jj2
        !epithelium
        k=cels(cell_grid_epi(i,j))%node(1)

        gex(k,1)=1d0
        gex(k,8)=1d0
        !borders
        if(i<=1.or.i>=lx.or.j<=1.or.j>=ly)then
          node(k)%hold=3
          node(node(k)%altre)%hold=3
        end if
        !placode
        if((i>px1.and.i<px2).and.(j>py1.and.j<py2))then
          gex(k,4)=1d0
          gex(k,13)=1d0
        else
          gex(k,5)=1d0
        end if
        !signalling centre
        if((i>sx1.and.i<sx2).and.(j>sy1.and.j<sy2))then
          gex(k,11)=1d0 ; gex(node(k)%altre,11)=1d0
          gex(k,22)=1d0 ; gex(node(k)%altre,22)=1d0
          gex(k,13)=0d0
        end if
        !AP proliferation supressor
        if((i>apx1.and.i<apx2).and.(j>apy1.and.j<apy2))then
          gex(k,23)=1d0
        end if


        !mesench
        do ii=1,layer
          k=cels(cell_grid_mes(i,j,ii))%node(1)
          gex(k,3)=1d0
          gex(k,10)=1d0
          !borders
          if(i<=1.or.i>=lx.or.j<=1.or.j>=ly)then
            node(k)%hold=3
          end if
          !placode
          if((i>px1.and.i<px2).and.(j>py1.and.j<py2))then
            gex(k,6)=1d0
            gex(k,15)=1d0
          else
            gex(k,7)=1d0
          end if
        end do
      end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!
    end if

  if(geometry==1)then  
    !ring=2  !scale
    ring=1 !tooth

    ring2=3

    polarized=0 !tooth
    !polarized=1 !scale

    do i=1,nd
      if(node(i)%tipus<3)then
        !if(node(i)%tipus==2)then
          !epithelial basal layer, basal side
          gex(i,1)=1d0  !receptors
          node(i)%marge=0
          d=sqrt(node(i)%x**2+node(i)%y**2)
          !setting expression on placode and interplacode
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,4)=1d0
            if(polarized==0.or.node(i)%x<=0.15)then
              gex(i,13)=1d0 !; gex(node(i)%altre,13)=1d0
            else
              !gex(i,23)=1d0
              gex(i,13)=1d0 !; gex(node(i)%altre,13)=1d0
            end if
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,5)=1d0
          end if
          !!!!!!!!!!!!!
          !signalling centre
          if(d-2*node(i)%req*(ring-1)<epsilod)then
            node(i)%marge=0 !; node(node(i)%altre)%marge=0
            gex(i,11)=1d0 !; gex(node(i)%altre,11)=1d0
            gex(i,22)=1d0 !; gex(node(i)%altre,22)=1d0
            gex(i,13)=0d0 !; gex(node(i)%altre,13)=0d0
          end if

        !end if
        !epithelial basal layer, whole layer       
        gex(i,8)=1d0

      end if
      if(node(i)%tipus==3)then
        if(node(i)%z<=0)then
          !mesenchymal layer
          gex(i,3)=1d0
          gex(i,10)=1d0

          !setting expression on placode and interplacode
          d=sqrt(node(i)%x**2+node(i)%y**2)
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,6)=1d0
            gex(i,15)=1d0
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,7)=1d0
          end if
          !!!!!!!!!!!!!

        else
          !epithelial suprabasal layer
          gex(i,2)=1d0
          gex(i,9)=1d0

          !setting expression on placode and interplacode
          d=sqrt(node(i)%x**2+node(i)%y**2)
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,4)=1d0
            if(polarized==0.or.node(i)%x<=0.15)then
              gex(i,14)=1d0
            else
              gex(i,23)=1d0
              gex(i,14)=1d0
            end if
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,5)=1d0
          end if
          !!!!!!!!!!!!!
          !signalling centre
          !if(d-2*node(i)%req*(ring-1)<epsilod)then
          !  gex(i,11)=1d0
          !  gex(i,22)=1d0
          !  gex(i,14)=0d0 !; gex(node(i)%altre,13)=0d0
          !end if


        end if
      end if


    end do
    end if

    call update_npag

    node(:)%talone=0.0d0

    !ramax=maxval(node(:)%da)*3

    
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


subroutine epidermal_organ_out_bud

!real*8::pericel,radio2,req1,req2,adjustreq,di
!integer::kk2

integer:: lx,ly
integer:: x1,x2,y1,y2,ii1,ii2,jj1,jj2
integer:: x3,x4,y3,y4,ii3,ii4,jj3,jj4,jj5,y5
integer, allocatable::cell_grid_epi(:,:),cell_grid_mes(:,:,:)
!integer:: periodic_borders,nborders
integer::geometry
integer::xlayer
integer::ring,ring2,polarized
integer::px1,px2,py1,py2,sx1,sx2,sy1,sy2,apx1,apx2,apy1,apy2

!******* #1 DEFINING SPATIAL DIMENSIONS *******

  geometry=2 !geometry=1 means hexagonal epithelium ; geometry=2 means rectangular epithelium
  
  if (geometry==1)then
	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=10    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=10   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=-0.5d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    do i=1,layer
      mesradicel(i)=mradicel
    end do
    
    xlayer=0 !this will make a layer of ECM beneath the layers of mesenchymal cells
    
  elseif(geometry==2)then
    layer=1      !number of planar cell layers
    lx=16
    ly=32
    
    allocate(cell_grid_epi(lx,ly),cell_grid_mes(lx,ly,layer))

  end if  
    !periodic_borders=1
    


    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(geometry==1)then
      if(radi>0.and.radicel>0)then
        j=0
        do i=1,radi-1
          j=j+i
        end do
        nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
        nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13
  
        j=0
        do i=1,radicel-1
          j=j+i
        end do
        ncelsepi=(6*j+1)
        ndepi=nodecel*ncelsepi
      else
        ncelsepi=0
        ndepi=0
        nodecel=0
        nodecela=0
      end if
  
      if(mradi>0.and.mradicel>0.and.layer>0)then
        if(packed==1)then !if mesenchyme is packed, we use mradi differently
          j=0
          do i=1,mradi-1
            j=j+i
          end do
          nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
          nodecela=2*nodecel+1
        else
          nodecel=mradi
          nodecela=2*nodecel+1
        end if
  
        ncelsmes=0
        do k=1,layer
          j=0 !;print*,"mesradicel",mesradicel(k)
          do i=1,mesradicel(k)-1
            j=j+i
          end do
          ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
        end do
        ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
        if(radi==0.and.radicel==0)then
          nodecel=mradi
          nodecela=2*nodecel+1
        else if(nodecel<mradi)then
          nodecel=radi
          nodecela=nodecel*2+1
        end if
      else
        ndmes=0
        ncelsmes=0
      end if
      ndx=(6*j+1)*xlayer
      nd=ndepi+ndmes+ndx
      !ncels=ncelsepi+ncelsmes
      ncels=ncelsepi+ncelsmes+ndx
   print*,"ncels",ncels,ncelsepi,ncelsmes,ndx
      nda=nd+10
      ncals=ncels+10
    elseif(geometry==2)then

      ncelsepi=lx*ly
      ndepi=ncelsepi*2
      ncelsmes=ncelsepi*layer
      ndmes=ncelsmes
      nd=ndepi+ndmes
      ncels=ncelsepi+ncelsmes
      nodecel=2 ; nodecela=5
      nda=nd+10
      ncals=ncels+10
    end if
    

  !End initializing dimension parameters

  print*,"inicial nd",nd,"ncels",ncels,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=1
    !End Allocatations


    

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=0.85d0
    deltamin=1d-2
    khold=1d0
    angletor=0.05



    k_bu=5d0
    ramax=0.35d0
    k_press=5d-1
    m_xwall=zmes !tooth
    mi_xwall=0d0 !tooth

    !m_xwall=0d0  !SCALE
    !mi_xwall=-5d-1 !SCALE

    ndmax=9d4
    
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=1
    ffu(2)=1 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !buoyancy

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(14)=0 !physical boundaries (walls)
    ffu(17)=1 !volume conservation (epithelial)
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise
    ffu(24)=0


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    !if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(mod(i,2)/=0)then  !basal
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
          !node(i)%reqp=-0.05d0
          !node(i)%req=0.20d0
        !node(i)%kplast=1d0
        else                      !apical
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !node(i)%kplast=1d1
        !  !node(i)%reqp=0.05d0
        !  !node(i)%req=0.30d0
        end if
        node(i)%req=0.25d0
        node(i)%reqcr=0.25d0 !; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.125
        node(i)%da=node(i)%req*2.0; 
        node(i)%ke=5d0
        node(i)%tor=3d0
        node(i)%stor=1d0
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=1d1
        node(i)%kvol=5d1
      end do
    !end if

	!Mesenchyme
    if(layer>0)then
      do i=ndepi+1,ndepi+ndmes
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !scale rect
        !node(i)%rep=0d0 ; node(i)%repcel=1d0
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*2.0
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

	!Suprabasal cells
    if(xlayer>0)then
      do i=ndepi+ndmes+1,nd
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*2.0
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if
    
    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)

    node(:)%hold=0


  if(geometry==1)then
    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    !if(radi>0.and.radicel>0)               call epiteli_sphere1(radi,radicel,radius,di,zepi)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if
!print*,"tall0"

    !Mesenchymal layer
    ii=0
    do i=ndepi+ndmes+1,nd
      ii=ii+1
      j=i-ndx
      node(i)%x=node(j)%x ; node(i)%y=node(j)%y ; node(i)%z=-0.5d0
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
      node(i)%icel=node(j)%icel+ndx ; node(i)%tipus=3
      !node(i)%icel=ncelsepi+ncelsmes+ii ; node(i)%tipus=3

      node(i)%repcel=5d-1
      !allocate(cels(node(i)%icel)%node(nodela))
      cels(node(i)%icel)%nunodes=1
      cels(node(i)%icel)%node(1)=i
      node(i)%marge=0
      !print*,"i",i,"j",j,"icel i",node(i)%icel,"icel j",node(j)%icel
      !node(i)%border=1
    end do
    ncelsmes=ncelsmes+ndx ; ndx=0
    ndmes=ncelsmes
!print*,"tall1"
      
    j=0
    do i=1,radicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    do i=j+1,ndepi
      node(i)%hold=2 ;node(i)%repcel=1d0
    end do
    j=0
    do i=1,mradicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    k=0
    do i=1,mradicel-1
      k=k+i
    end do
    k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
      node(i)%hold=2;node(i)%repcel=1d0
      !node(i)%border=1
      !node(i)%da=node(i)%req*2.0
      !node(i)%orix=0 ; node(i)%oriy=0
      !node(i)%rep=1d1;node(i)%repcel=1d1
    end do
    do i=ndepi+ndmes/2+j+1,nd
      node(i)%hold=2;node(i)%repcel=1d0
      !node(i)%border=1
    end do
    !do i=ndepi+k+1,nd
    !  node(i)%hold=2;node(i)%repcel=1d0
    !  node(i)%border=1
    !end do

    
  elseif(geometry==2)then
	do i=1,ncels
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do


    !epithelial cells
    !origin of rectangle
    aa=node(1)%req*2
    bb=sqrt(aa**2-(aa*0.5)**2)
    a=-real(lx)*0.5*aa ; b=-real(ly)*0.5*bb
    !print*,"a",a,"b",b,"aa",aa,"bb",bb,"lx",lx,"ly",ly
    ii=0 ; jj=0
    do i=1,ly
      do j=1,lx
        ii=ii+1 ;jj=jj+1 ; cell_grid_epi(j,i)=jj
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii+1 ;node(ii)%tipus=2
        ii=ii+1
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii-1 ;node(ii)%tipus=1
        cels(jj)%node(1)=ii-1 ; cels(jj)%node(2)=ii ;cels(jj)%nunodes=2; cels(jj)%ctipus=1
        if(i==1.or.i==ly.or.j==1.or.j==lx)then
          node(ii)%hold=2;node(ii-1)%hold=2 ;
          node(ii)%repcel=1d1;node(ii-1)%repcel=1d1;
          !node(ii)%border=1;node(ii-1)%border=1;
        end if
      end do
    end do

    !mesenchymals
    do k=1,layer
      !print*,"k",k
      do i=1,ly
        do j=1,lx
          ii=ii+1 ; jj=jj+1 ; cell_grid_mes(j,i,k)=jj
          if(mod(i,2)==0)then
            node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          else
            node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          end if
          node(ii)%icel=jj ; node(ii)%altre=0 ;node(ii)%tipus=3
          cels(jj)%node(1)=ii ;cels(jj)%nunodes=1 ; cels(jj)%ctipus=3
          if(i==1 .or. i==ly .or. j==1 .or. j==lx)then
            node(ii)%hold=2 ;node(ii)%repcel=1d1 ; !node(ii)%border=2 ;
          end if
          !if(i==1 .or. i==ly .or. j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
          !if(j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
        end do
      end do
    end do
  end if

  !m_xwall=maxval(node(1:nd)%x)+0.01
  !mi_xwall=minval(node(1:nd)%x)-0.01
  !m_ywall=maxval(node(1:nd)%y)+0.01
  !mi_ywall=minval(node(1:nd)%y)-0.01
  
!print*,"tall2"

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
  !node(:)%hold=2
  !node(:)%hold=0  

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    !if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    !end if
    !Mesenchymal
    if(layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=23
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=0d0 ; gen(1)%mu=5d-1 ;gen(1)%name="epithelial basal marker"
      gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=5d-1 ;gen(2)%name="epithelial suprabasal marker"
      gen(3)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=5d-1 ;gen(3)%name="mesenchymal marker"
      gen(4)%kindof=1 ; gen(4)%diffu=0d0 ; gen(4)%mu=5d-1 ;gen(4)%name="epi. placode marker"
      gen(5)%kindof=1 ; gen(5)%diffu=0d0 ; gen(5)%mu=5d-1 ;gen(5)%name="epi. interplacode marker"
      gen(6)%kindof=1 ; gen(6)%diffu=0d0 ; gen(6)%mu=5d-1 ;gen(6)%name="mes. placode marker"
      gen(7)%kindof=1 ; gen(7)%diffu=0d0 ; gen(7)%mu=5d-1 ;gen(7)%name="mes. interplacode marker"

      gen(8)%kindof=1 ; gen(8)%diffu=0d0 ; gen(8)%mu=5d-1 ;gen(8)%name="epi. basal adhesion molecule"
      gen(9)%kindof=1 ; gen(9)%diffu=0d0 ; gen(9)%mu=5d-1 ;gen(9)%name="epi. suprabasal adhesion molecule"
      gen(10)%kindof=1 ; gen(10)%diffu=0d0 ; gen(10)%mu=5d-1 ;gen(10)%name="mesench. adhesion molecule"

      gen(11)%kindof=1 ; gen(11)%diffu=0d0 ; gen(11)%mu=5d-1 ;gen(11)%name="epi. signalling centre marker"
      gen(12)%kindof=4 ; gen(12)%diffu=5d-1 ; gen(12)%mu=5d-1 ;gen(12)%name="epi. signal"

      gen(13)%kindof=2 ; gen(13)%diffu=0d0 ; gen(13)%mu=5d-1 ;gen(13)%name="epi. basal receptor inactive"
      gen(14)%kindof=2 ; gen(14)%diffu=0d0 ; gen(14)%mu=5d-1 ;gen(14)%name="epi. suprabasal receptor inactive"
      gen(15)%kindof=2 ; gen(15)%diffu=0d0 ; gen(15)%mu=5d-1 ;gen(15)%name="mesench. receptor inactive"

      gen(16)%kindof=8 ; gen(16)%diffu=0d0 ; gen(16)%mu=5d-1 ;gen(16)%name="epi. basal receptor active"
      gen(16)%npre=2 ; allocate(gen(16)%pre(gen(16)%npre)) ; gen(16)%pre(1)=12 ; gen(16)%pre(2)=13
      gen(16)%npost=2 ; allocate(gen(16)%post(gen(16)%npost)) ; gen(16)%post(1)=12 ; gen(16)%post(2)=13

      gen(17)%kindof=8 ; gen(17)%diffu=0d0 ; gen(17)%mu=5d-1 ;gen(17)%name="epi. suprabasal receptor active"
      gen(17)%npre=2 ; allocate(gen(17)%pre(gen(17)%npre)) ; gen(17)%pre(1)=12 ; gen(17)%pre(2)=14
      gen(17)%npost=2 ; allocate(gen(17)%post(gen(17)%npost)) ; gen(17)%post(1)=12 ; gen(17)%post(2)=14

      gen(18)%kindof=8 ; gen(18)%diffu=0d0 ; gen(18)%mu=5d-1 ;gen(18)%name="mesench. receptor active"
      gen(18)%npre=2 ; allocate(gen(18)%pre(gen(18)%npre)) ; gen(18)%pre(1)=12 ; gen(18)%pre(2)=15
      gen(18)%npost=2 ; allocate(gen(18)%post(gen(18)%npost)) ; gen(18)%post(1)=12 ; gen(18)%post(2)=15

      gen(19)%kindof=1 ; gen(19)%diffu=0d0 ; gen(19)%mu=5d-1 ;gen(19)%name="epi. basal effector"
      gen(20)%kindof=1 ; gen(20)%diffu=0d0 ; gen(20)%mu=5d-1 ;gen(20)%name="epi. suprabasal effector"
      gen(21)%kindof=1 ; gen(21)%diffu=0d0 ; gen(21)%mu=5d-1 ;gen(21)%name="mesenchymal effector"

      gen(22)%kindof=1 ; gen(22)%diffu=0d0 ; gen(22)%mu=5d-1 ;gen(22)%name="signalling centre adh. mol."
      gen(23)%kindof=1 ; gen(23)%diffu=0d0 ; gen(23)%mu=5d-1 ;gen(23)%name="distal (posterior) gene product"
      
      !gen(24)%kindof=1 ; gen(24)%diffu=0d0 ; gen(24)%mu=5d-1 ;gen(24)%name="s2 IEE adhesion molecule."
      !gen(25)%kindof=1 ; gen(25)%diffu=0d0 ; gen(25)%mu=5d-1 ;gen(25)%name="s2 dermal papilla adhesion molecule"
    
 
    !Gene-behavior interactions
    
   
      gen(8)%wa(1)=1
      gen(9)%wa(1)=2
      gen(10)%wa(1)=3
      gen(22)%wa(1)=4
      !gen(24)%wa(1)=5
      !gen(25)%wa(1)=6

      !gen(19)%wa(1)=5
      !gen(20)%wa(1)=6
      !gen(21)%wa(1)=7
      
      gen(22)%wa(10)=1d0
      
      gen(11)%wa(nparam_per_node+2)=-5d5

      gen(19)%wa(nparam_per_node+2)=1.4d-1/2d0
      !gen(20)%wa(nparam_per_node+2)=2.0d-1
      gen(21)%wa(nparam_per_node+2)=0.7d-1

      gen(23)%wa(nparam_per_node+2)=-5.0d-2


      !gen(13)%wa(nparam_per_node+2)=2.0d-1/2d0
      !gen(14)%wa(nparam_per_node+2)=2.0d-1
      !gen(15)%wa(nparam_per_node+2)=5.0d-2

      
      !growth parameters for inward bud
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=3.2d-1
      !gen(21)%wa(nparam_per_node+2)=5d-2
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !growth parameters for outward bud
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=5.0d-3!5.0d-2
      !gen(21)%wa(nparam_per_node+2)=3.2d-1
      !gen(23)%wa(nparam_per_node+2)=-5.0d0
      !!!!!!!!!!!!!!!!!!!!!!!!!!   
      !growth parameters for rectangle scale
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=0.0d0
      !gen(21)%wa(nparam_per_node+2)=1.0d-3
      !gen(23)%wa(nparam_per_node+2)=-7d1
      !gen(6)%wa(nparam_per_node+2)=1.0d-1
      !!!!!!!!!!!!!!!!!!!!!!!!!!   
      
    !Gene-gene interactions

      !structural genes
      !auto-maintenance
      gen(1)%w(1)=1d0
      gen(2)%w(2)=1d0
      gen(3)%w(3)=1d0
      gen(4)%w(4)=1d0
      gen(5)%w(5)=1d0
      gen(6)%w(6)=1d0
      gen(7)%w(7)=1d0
      gen(11)%w(11)=1d0
      !maintenance of stable gene expression
      gen(13)%w(4)=1d0
      gen(14)%w(4)=1d0
      gen(15)%w(6)=1d0
      gen(8)%w(1)=1d0
      gen(9)%w(2)=1d0
      gen(10)%w(3)=1d0
      gen(22)%w(11)=1d0
      gen(8)%w(11)=-1d5
      gen(9)%w(11)=-1d5
      !!!!!!!!!!!!!!!!
      gen(12)%w(11)=1d0
      gen(13)%w(11)=-1d5
      gen(14)%w(11)=-1d5
      gen(14)%w(1)=-1d5
      gen(13)%w(2)=-1d5
      !gen(24)%w(11)=1d0 !s2
      !!!!!!!!!!!!!!!!! effector activation thresholds
      gen(19)%w(1)=-0.05d0
      gen(20)%w(2)=-0.05d0
      gen(21)%w(3)=-0.05d0
      !gen(19)%w(1)=0d0
      !gen(20)%w(2)=0d0
      !gen(21)%w(3)=0d0
      gen(19)%w(16)=5d1   !tooth
      gen(20)%w(17)=5d1
      gen(21)%w(18)=5d1
      !gen(19)%w(16)=1d2    !scale
      !gen(20)%w(17)=1d2
      !gen(21)%w(18)=1d1
      !gen(19)%w(16)=1d1   !tooth
      !gen(20)%w(17)=1d1
      !gen(21)%w(18)=1d1
      !gen(19)%w(16)=1d2    !rectangle scale
      !gen(20)%w(17)=1d2
      !gen(21)%w(18)=1d1

      !gen(24)%w(19)=1d0 !s2
      !gen(25)%w(21)=1d0 !s2


      !effector supression of adhesion molecules
      !gen(8)%w(19)=-1d1
      !gen(9)%w(20)=-1d1
      !gen(10)%w(21)=-1d1

      !!!!!!!!!!!!!!!! a gene expressed in the posterior half of the placode supresses signalling
      gen(23)%w(23)=1d0
      !gen(13)%w(23)=-1d5
      !gen(14)%w(23)=-1d5

    


      !signal-receptor interactions
      !epithelial basal binding
      gen(16)%nww=4
      gen(16)%ww(1,1)=12
      gen(16)%ww(1,2)=16
      gen(16)%ww(1,3)=1d0
      gen(16)%ww(2,1)=16
      gen(16)%ww(2,2)=12
      gen(16)%ww(2,3)=1d0
      gen(16)%ww(3,1)=13
      gen(16)%ww(3,2)=16
      gen(16)%ww(3,3)=1d0
      gen(16)%ww(4,1)=16
      gen(16)%ww(4,2)=13
      gen(16)%ww(4,3)=1d0
      !epithelial suprabasal receptor binding
      gen(17)%nww=4
      gen(17)%ww(1,1)=12
      gen(17)%ww(1,2)=17
      gen(17)%ww(1,3)=1d0
      gen(17)%ww(2,1)=17
      gen(17)%ww(2,2)=12
      gen(17)%ww(2,3)=1d0
      gen(17)%ww(3,1)=14
      gen(17)%ww(3,2)=15
      gen(17)%ww(3,3)=1d0
      gen(17)%ww(4,1)=17
      gen(17)%ww(4,2)=14
      gen(17)%ww(4,3)=1d0
      !mesenchymal receptor binding
      gen(18)%nww=4
      gen(18)%ww(1,1)=12
      gen(18)%ww(1,2)=18
      gen(18)%ww(1,3)=1d0
      gen(18)%ww(2,1)=18
      gen(18)%ww(2,2)=12
      gen(18)%ww(2,3)=1d0
      gen(18)%ww(3,1)=15
      gen(18)%ww(3,2)=18
      gen(18)%ww(3,3)=1d0
      gen(18)%ww(4,1)=18
      gen(18)%ww(4,2)=15
      gen(18)%ww(4,3)=1d0
      


    !Adhesion molecules

	ntipusadh=6
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      !as default, 1 and 2 should act as the same (basal and suprabasal)
      kadh(1,1)=1d0

      kadh(1,2)=0d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=0d0

      kadh(1,3)=5d-1 ; kadh(3,1)=kadh(1,3)
      !rectangle scale
      !kadh(1,3)=3d0 ; kadh(3,1)=kadh(1,3)
      kadh(2,3)=0d0 ; kadh(3,2)=kadh(2,3)
      kadh(3,3)=1d1

      kadh(1,4)=5d-1 ; kadh(4,1)=kadh(1,4)
      kadh(2,4)=0d0 ; kadh(4,2)=kadh(2,4)
      kadh(3,4)=5d-1 ; kadh(4,3)=kadh(3,4)
      !kadh(4,4)=1d1 !tooth
      kadh(4,4)=1d1  !scale

      !kadh(1,5)=0d0 ; kadh(5,1)=kadh(1,5)
      !kadh(2,5)=0d0 ; kadh(5,2)=kadh(2,5)
      !kadh(3,5)=0d0 ; kadh(5,3)=kadh(3,5)
      !kadh(4,5)=0d0 ; kadh(5,4)=kadh(4,5)
      !kadh(5,5)=0d0

      !kadh(1,6)=0d0 ; kadh(6,1)=kadh(1,6)
      !kadh(2,6)=0d0 ; kadh(6,2)=kadh(2,6)
      !kadh(3,6)=0d0 ; kadh(6,3)=kadh(3,6)
      !kadh(4,6)=0d0 ; kadh(6,4)=kadh(4,6)
      !kadh(5,6)=1d1 ; kadh(6,5)=kadh(5,6)
      !kadh(6,6)=1d1

      !kadh(1,7)=5d-1 ; kadh(7,1)=kadh(1,7)
      !kadh(2,7)=0d0 ; kadh(7,2)=kadh(2,7)
      !kadh(3,7)=1d0 ; kadh(7,3)=kadh(3,7)
      !kadh(4,7)=5d-1 ; kadh(7,4)=kadh(4,7)
      !kadh(5,7)=5d-1 ; kadh(7,5)=kadh(5,7)
      !kadh(6,7)=0d0 ; kadh(7,6)=kadh(6,7)
      !kadh(7,7)=1d0

    end if

    !Gene expression on nodes

    
    !if(geometry==1)then
    !  gex(1,2)=1d0
    !  gex(ndepi+1,2)=1d0
    !  !gex(ndepi+1,2)=1d0
    !end if
    
    if(geometry==2)then
    !draw rectangles
    !!!!!!!!!!!!!!!!!!!!!!!!!
    jj1=1 ; jj2=ly
    ii1=1 ; ii2=lx
    !placode range
    px1=3 ; px2=lx-(px1-1)
    py1=3 ; py2=ly-(py1-1)
    !signalling centre
    sx1=6 ; sx2=lx-(sx1-1)
    sy1=5 ; sy2=ly-(sy1-1)
    !AP proliferation supressor
    apx1=6 ; apx2=lx
    apy1=1 ; apy2=ly

    do i=ii1,ii2
      do j=jj1,jj2
        !epithelium
        k=cels(cell_grid_epi(i,j))%node(1)

        gex(k,1)=1d0 ; gex(node(k)%altre,1)=1d0
        gex(k,8)=1d0 ; gex(node(k)%altre,8)=1d0
        !borders
        if(i<=1.or.i>=lx.or.j<=1.or.j>=ly)then
          node(k)%hold=2
          node(node(k)%altre)%hold=2
        end if
        !placode
        if((i>px1.and.i<px2).and.(j>py1.and.j<py2))then
          gex(k,4)=1d0 ; gex(node(k)%altre,4)=1d0
          gex(k,13)=1d0 ; gex(node(k)%altre,13)=1d0
        else
          gex(k,5)=1d0 ; gex(node(k)%altre,5)=1d0
        end if
        !signalling centre
        if((i>sx1.and.i<sx2).and.(j>sy1.and.j<sy2))then
          gex(k,11)=1d0 ; gex(node(k)%altre,11)=1d0
          gex(k,22)=1d0 ; gex(node(k)%altre,22)=1d0
          gex(k,13)=0d0 ; gex(node(k)%altre,13)=0d0
        end if
        !AP proliferation supressor
        if((i>apx1.and.i<apx2).and.(j>apy1.and.j<apy2))then
          gex(k,23)=1d0 ; gex(node(k)%altre,23)=1d0
        end if


        !mesench
        do ii=1,layer
          k=cels(cell_grid_mes(i,j,ii))%node(1)
          gex(k,3)=1d0
          gex(k,10)=1d0
          !borders
          if(i<=1.or.i>=lx.or.j<=1.or.j>=ly)then
            node(k)%hold=2
          end if
          !placode
          if((i>px1.and.i<px2).and.(j>py1.and.j<py2))then
            gex(k,6)=1d0
            gex(k,15)=1d0
          else
            gex(k,7)=1d0
          end if
        end do
      end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!
    end if

  if(geometry==1)then  

    ring=3 

    ring2=7

    !polarized=0 !tooth
    polarized=1 !scale

    do i=1,nd
      if(node(i)%tipus<3)then
        !if(node(i)%tipus==2)then
          !epithelial basal layer, basal side
          gex(i,1)=1d0  !receptors
          node(i)%marge=0
          d=sqrt(node(i)%x**2+node(i)%y**2)
          !setting expression on placode and interplacode
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,4)=1d0
            if(polarized==0.or.node(i)%x<=-0.75)then
              gex(i,13)=1d0 ; gex(node(i)%altre,13)=1d0
            else
              gex(i,23)=1d0
              gex(i,13)=1d0 ; gex(node(i)%altre,13)=1d0
            end if
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,5)=1d0
          end if
          !!!!!!!!!!!!!
          !signalling centre
          if(d-2*node(i)%req*(ring-1)<epsilod)then
            node(i)%marge=0 ; node(node(i)%altre)%marge=0
            gex(i,11)=1d0 ; gex(node(i)%altre,11)=1d0
            gex(i,22)=1d0 ; gex(node(i)%altre,22)=1d0
            gex(i,13)=0d0 ; gex(node(i)%altre,13)=0d0
          end if

        !end if
        !epithelial basal layer, whole layer       
        gex(i,8)=1d0

      end if
      if(node(i)%tipus==3)then
        if(node(i)%z<=0)then
          !mesenchymal layer
          gex(i,3)=1d0
          gex(i,10)=1d0

          !setting expression on placode and interplacode
          d=sqrt(node(i)%x**2+node(i)%y**2)
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,6)=1d0
            gex(i,15)=1d0
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,7)=1d0
          end if
          !!!!!!!!!!!!!

        else
          !epithelial suprabasal layer
          gex(i,2)=1d0
          gex(i,9)=1d0

          !setting expression on placode and interplacode
          d=sqrt(node(i)%x**2+node(i)%y**2)
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,4)=1d0
            if(polarized==0.or.node(i)%x<=0.15)then
              gex(i,14)=1d0
            else
              gex(i,23)=1d0
              gex(i,14)=1d0
            end if
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,5)=1d0
          end if
          !!!!!!!!!!!!!
          !signalling centre
          !if(d-2*node(i)%req*(ring-1)<epsilod)then
          !  gex(i,11)=1d0
          !  gex(i,22)=1d0
          !  gex(i,14)=0d0 !; gex(node(i)%altre,13)=0d0
          !end if


        end if
      end if


    end do
    end if

    call update_npag

    node(:)%talone=0.0d0

    !ramax=maxval(node(:)%da)*3

    
end subroutine

subroutine epidermal_scale

!real*8::pericel,radio2,req1,req2,adjustreq,di
!integer::kk2

integer:: lx,ly
integer:: x1,x2,y1,y2,ii1,ii2,jj1,jj2
integer:: x3,x4,y3,y4,ii3,ii4,jj3,jj4,jj5,y5
integer, allocatable::cell_grid_epi(:,:),cell_grid_mes(:,:,:)
!integer:: periodic_borders,nborders
integer::geometry
integer::xlayer
integer::ring,ring2,polarized
integer::px1,px2,py1,py2,sx1,sx2,sy1,sy2,apx1,apx2,apy1,apy2

!******* #1 DEFINING SPATIAL DIMENSIONS *******

  geometry=1 !geometry=1 means hexagonal epithelium ; geometry=2 means rectangular epithelium
  
  if (geometry==1)then
	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=10    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=10   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=-0.5d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    do i=1,layer
      mesradicel(i)=mradicel
    end do
    
    xlayer=0 !this will make a layer of ECM beneath the layers of mesenchymal cells
    
  elseif(geometry==2)then
    layer=1      !number of planar cell layers
    lx=16
    ly=32
    
    allocate(cell_grid_epi(lx,ly),cell_grid_mes(lx,ly,layer))

  end if  
    !periodic_borders=1
    


    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(geometry==1)then
      if(radi>0.and.radicel>0)then
        j=0
        do i=1,radi-1
          j=j+i
        end do
        nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
        nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13
  
        j=0
        do i=1,radicel-1
          j=j+i
        end do
        ncelsepi=(6*j+1)
        ndepi=nodecel*ncelsepi
      else
        ncelsepi=0
        ndepi=0
        nodecel=0
        nodecela=0
      end if
  
      if(mradi>0.and.mradicel>0.and.layer>0)then
        if(packed==1)then !if mesenchyme is packed, we use mradi differently
          j=0
          do i=1,mradi-1
            j=j+i
          end do
          nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
          nodecela=2*nodecel+1
        else
          nodecel=mradi
          nodecela=2*nodecel+1
        end if
  
        ncelsmes=0
        do k=1,layer
          j=0 !;print*,"mesradicel",mesradicel(k)
          do i=1,mesradicel(k)-1
            j=j+i
          end do
          ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
        end do
        ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
        if(radi==0.and.radicel==0)then
          nodecel=mradi
          nodecela=2*nodecel+1
        else if(nodecel<mradi)then
          nodecel=radi
          nodecela=nodecel*2+1
        end if
      else
        ndmes=0
        ncelsmes=0
      end if
      ndx=(6*j+1)*xlayer
      nd=ndepi+ndmes+ndx
      !ncels=ncelsepi+ncelsmes
      ncels=ncelsepi+ncelsmes+ndx
   print*,"ncels",ncels,ncelsepi,ncelsmes,ndx
      nda=nd+10
      ncals=ncels+10
    elseif(geometry==2)then

      ncelsepi=lx*ly
      ndepi=ncelsepi*2
      ncelsmes=ncelsepi*layer
      ndmes=ncelsmes
      nd=ndepi+ndmes
      ncels=ncelsepi+ncelsmes
      nodecel=2 ; nodecela=5
      nda=nd+10
      ncals=ncels+10
    end if
    

  !End initializing dimension parameters

  print*,"inicial nd",nd,"ncels",ncels,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=1
    !End Allocatations


    

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=0.85d0
    deltamin=1d-2
    khold=1d0
    angletor=0.05



    k_bu=5d0
    ramax=0.35d0
    k_press=5d-1
    m_xwall=zmes !tooth
    mi_xwall=0d0 !tooth

    !m_xwall=0d0  !SCALE
    !mi_xwall=-5d-1 !SCALE

    ndmax=9d4
    
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=1
    ffu(2)=1 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !buoyancy

    ffu(11)=1 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(14)=0 !physical boundaries (walls)
    ffu(17)=1 !volume conservation (epithelial)
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise
    ffu(24)=0


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    !if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(mod(i,2)/=0)then  !basal
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
          !node(i)%reqp=-0.05d0
          !node(i)%req=0.20d0
        !node(i)%kplast=1d0
        else                      !apical
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !node(i)%kplast=1d1
        !  !node(i)%reqp=0.05d0
        !  !node(i)%req=0.30d0
        end if
        node(i)%req=0.25d0
        node(i)%reqcr=0.25d0 !; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.125
        node(i)%da=node(i)%req*2.0; 
        node(i)%ke=5d0
        node(i)%tor=3d0
        node(i)%stor=1d0
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=1d1
        node(i)%kvol=5d1
      end do
    !end if

	!Mesenchyme
    if(layer>0)then
      do i=ndepi+1,ndepi+ndmes
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !scale rect
        !node(i)%rep=0d0 ; node(i)%repcel=1d0
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*2.0
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

	!Suprabasal cells
    if(xlayer>0)then
      do i=ndepi+ndmes+1,nd
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*2.0
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if
    
    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)

    node(:)%hold=0


  if(geometry==1)then
    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    !if(radi>0.and.radicel>0)               call epiteli_sphere1(radi,radicel,radius,di,zepi)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if
!print*,"tall0"

    !Mesenchymal layer
    ii=0
    do i=ndepi+ndmes+1,nd
      ii=ii+1
      j=i-ndx
      node(i)%x=node(j)%x ; node(i)%y=node(j)%y ; node(i)%z=-0.5d0
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
      node(i)%icel=node(j)%icel+ndx ; node(i)%tipus=3
      !node(i)%icel=ncelsepi+ncelsmes+ii ; node(i)%tipus=3

      node(i)%repcel=5d-1
      !allocate(cels(node(i)%icel)%node(nodela))
      cels(node(i)%icel)%nunodes=1
      cels(node(i)%icel)%node(1)=i
      node(i)%marge=0
      !print*,"i",i,"j",j,"icel i",node(i)%icel,"icel j",node(j)%icel
      !node(i)%border=1
    end do
    ncelsmes=ncelsmes+ndx ; ndx=0
    ndmes=ncelsmes
!print*,"tall1"
      
    j=0
    do i=1,radicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    do i=j+1,ndepi
      node(i)%hold=3 ;node(i)%repcel=5d-1
    end do
    j=0
    do i=1,mradicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    k=0
    do i=1,mradicel-1
      k=k+i
    end do
    k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
      node(i)%hold=3;node(i)%repcel=5d-1
      !node(i)%border=1
      !node(i)%da=node(i)%req*2.0
      !node(i)%orix=0 ; node(i)%oriy=0
      !node(i)%rep=1d1;node(i)%repcel=1d1
    end do
    do i=ndepi+ndmes/2+j+1,nd
      node(i)%hold=3;node(i)%repcel=5d-1
      !node(i)%border=1
    end do
    !do i=ndepi+k+1,nd
    !  node(i)%hold=2;node(i)%repcel=1d0
    !  node(i)%border=1
    !end do

    
  elseif(geometry==2)then
	do i=1,ncels
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do


    !epithelial cells
    !origin of rectangle
    aa=node(1)%req*2
    bb=sqrt(aa**2-(aa*0.5)**2)
    a=-real(lx)*0.5*aa ; b=-real(ly)*0.5*bb
    !print*,"a",a,"b",b,"aa",aa,"bb",bb,"lx",lx,"ly",ly
    ii=0 ; jj=0
    do i=1,ly
      do j=1,lx
        ii=ii+1 ;jj=jj+1 ; cell_grid_epi(j,i)=jj
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii+1 ;node(ii)%tipus=2
        ii=ii+1
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii-1 ;node(ii)%tipus=1
        cels(jj)%node(1)=ii-1 ; cels(jj)%node(2)=ii ;cels(jj)%nunodes=2; cels(jj)%ctipus=1
        if(i==1.or.i==ly.or.j==1.or.j==lx)then
          node(ii)%hold=2;node(ii-1)%hold=2 ;
          node(ii)%repcel=1d1;node(ii-1)%repcel=1d1;
          !node(ii)%border=1;node(ii-1)%border=1;
        end if
      end do
    end do

    !mesenchymals
    do k=1,layer
      !print*,"k",k
      do i=1,ly
        do j=1,lx
          ii=ii+1 ; jj=jj+1 ; cell_grid_mes(j,i,k)=jj
          if(mod(i,2)==0)then
            node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          else
            node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          end if
          node(ii)%icel=jj ; node(ii)%altre=0 ;node(ii)%tipus=3
          cels(jj)%node(1)=ii ;cels(jj)%nunodes=1 ; cels(jj)%ctipus=3
          if(i==1 .or. i==ly .or. j==1 .or. j==lx)then
            node(ii)%hold=2 ;node(ii)%repcel=1d1 ; !node(ii)%border=2 ;
          end if
          !if(i==1 .or. i==ly .or. j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
          !if(j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
        end do
      end do
    end do
  end if

  !m_xwall=maxval(node(1:nd)%x)+0.01
  !mi_xwall=minval(node(1:nd)%x)-0.01
  !m_ywall=maxval(node(1:nd)%y)+0.01
  !mi_ywall=minval(node(1:nd)%y)-0.01
  
!print*,"tall2"

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
  !node(:)%hold=2
  !node(:)%hold=0  

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    !if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    !end if
    !Mesenchymal
    if(layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=30
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=0d0 ; gen(1)%mu=5d-1 ;gen(1)%name="epithelial basal marker"
      gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=5d-1 ;gen(2)%name="epithelial suprabasal marker"
      gen(3)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=5d-1 ;gen(3)%name="mesenchymal marker"
      gen(4)%kindof=1 ; gen(4)%diffu=0d0 ; gen(4)%mu=5d-1 ;gen(4)%name="epi. placode marker"
      gen(5)%kindof=1 ; gen(5)%diffu=0d0 ; gen(5)%mu=5d-1 ;gen(5)%name="epi. interplacode marker"
      gen(6)%kindof=1 ; gen(6)%diffu=0d0 ; gen(6)%mu=5d-1 ;gen(6)%name="mes. placode marker"
      gen(7)%kindof=1 ; gen(7)%diffu=0d0 ; gen(7)%mu=5d-1 ;gen(7)%name="mes. interplacode marker"

      gen(8)%kindof=1 ; gen(8)%diffu=0d0 ; gen(8)%mu=5d-1 ;gen(8)%name="epi. basal adhesion molecule"
      gen(9)%kindof=1 ; gen(9)%diffu=0d0 ; gen(9)%mu=5d-1 ;gen(9)%name="epi. suprabasal adhesion molecule"
      gen(10)%kindof=1 ; gen(10)%diffu=0d0 ; gen(10)%mu=5d-1 ;gen(10)%name="mesench. adhesion molecule"

      gen(11)%kindof=1 ; gen(11)%diffu=0d0 ; gen(11)%mu=5d-1 ;gen(11)%name="epi. signalling centre marker"
      gen(12)%kindof=4 ; gen(12)%diffu=5d-1 ; gen(12)%mu=5d-1 ;gen(12)%name="epi. signal"

      gen(13)%kindof=2 ; gen(13)%diffu=0d0 ; gen(13)%mu=5d-1 ;gen(13)%name="epi. basal receptor inactive"
      gen(14)%kindof=2 ; gen(14)%diffu=0d0 ; gen(14)%mu=5d-1 ;gen(14)%name="epi. suprabasal receptor inactive"
      gen(15)%kindof=2 ; gen(15)%diffu=0d0 ; gen(15)%mu=5d-1 ;gen(15)%name="mesench. receptor inactive"

      gen(16)%kindof=8 ; gen(16)%diffu=0d0 ; gen(16)%mu=5d-1 ;gen(16)%name="epi. basal receptor active"
      gen(16)%npre=2 ; allocate(gen(16)%pre(gen(16)%npre)) ; gen(16)%pre(1)=12 ; gen(16)%pre(2)=13
      gen(16)%npost=2 ; allocate(gen(16)%post(gen(16)%npost)) ; gen(16)%post(1)=12 ; gen(16)%post(2)=13

      gen(17)%kindof=8 ; gen(17)%diffu=0d0 ; gen(17)%mu=5d-1 ;gen(17)%name="epi. suprabasal receptor active"
      gen(17)%npre=2 ; allocate(gen(17)%pre(gen(17)%npre)) ; gen(17)%pre(1)=12 ; gen(17)%pre(2)=14
      gen(17)%npost=2 ; allocate(gen(17)%post(gen(17)%npost)) ; gen(17)%post(1)=12 ; gen(17)%post(2)=14

      gen(18)%kindof=8 ; gen(18)%diffu=0d0 ; gen(18)%mu=5d-1 ;gen(18)%name="mesench. receptor active"
      gen(18)%npre=2 ; allocate(gen(18)%pre(gen(18)%npre)) ; gen(18)%pre(1)=12 ; gen(18)%pre(2)=15
      gen(18)%npost=2 ; allocate(gen(18)%post(gen(18)%npost)) ; gen(18)%post(1)=12 ; gen(18)%post(2)=15

      gen(19)%kindof=1 ; gen(19)%diffu=0d0 ; gen(19)%mu=5d-1 ;gen(19)%name="epi. basal effector"
      gen(20)%kindof=1 ; gen(20)%diffu=0d0 ; gen(20)%mu=5d-1 ;gen(20)%name="epi. suprabasal effector"
      gen(21)%kindof=1 ; gen(21)%diffu=0d0 ; gen(21)%mu=5d-1 ;gen(21)%name="mesenchymal effector"

      gen(22)%kindof=1 ; gen(22)%diffu=0d0 ; gen(22)%mu=5d-1 ;gen(22)%name="signalling centre adh. mol."
      gen(23)%kindof=1 ; gen(23)%diffu=0d0 ; gen(23)%mu=5d-1 ;gen(23)%name="distal (posterior) gene product"

      gen(24)%kindof=4 ; gen(24)%diffu=5d-1 ; gen(24)%mu=5d-1 ;gen(24)%name="Wnt posterior signal"
      gen(25)%kindof=1 ; gen(25)%diffu=0d0 ; gen(25)%mu=5d-1 ;gen(25)%name="Shh territory epi"
      gen(26)%kindof=1 ; gen(26)%diffu=0d0 ; gen(26)%mu=5d-1 ;gen(26)%name="beta-cat positive mesench territory"


      gen(27)%kindof=2 ; gen(27)%diffu=0d0 ; gen(27)%mu=5d-1 ;gen(27)%name="Wnt epi receptor inactive"
      gen(28)%kindof=2 ; gen(28)%diffu=0d0 ; gen(28)%mu=5d-1 ;gen(28)%name="Wnt mes receptor inactive"

      gen(29)%kindof=8 ; gen(29)%diffu=0d0 ; gen(29)%mu=5d-1 ;gen(29)%name="Wnt epi receptor active"
      gen(29)%npre=2 ; allocate(gen(29)%pre(gen(29)%npre)) ; gen(29)%pre(1)=24 ; gen(29)%pre(2)=27
      gen(29)%npost=2 ; allocate(gen(29)%post(gen(29)%npost)) ; gen(29)%post(1)=24 ; gen(29)%post(2)=27

      gen(30)%kindof=8 ; gen(30)%diffu=0d0 ; gen(30)%mu=5d-1 ;gen(30)%name="Wnt mes receptor active"
      gen(30)%npre=2 ; allocate(gen(30)%pre(gen(30)%npre)) ; gen(30)%pre(1)=24 ; gen(30)%pre(2)=28
      gen(30)%npost=2 ; allocate(gen(30)%post(gen(30)%npost)) ; gen(30)%post(1)=24 ; gen(30)%post(2)=28
    
 


    !Gene-behavior interactions
    
   
      gen(8)%wa(1)=1
      gen(9)%wa(1)=2
      gen(10)%wa(1)=3
      gen(22)%wa(1)=4
      !gen(24)%wa(1)=5
      !gen(25)%wa(1)=6

      !gen(19)%wa(1)=5
      !gen(20)%wa(1)=6
      !gen(21)%wa(1)=7
      
      gen(22)%wa(10)=1d0
      
      gen(11)%wa(nparam_per_node+2)=-5d5

      gen(19)%wa(nparam_per_node+2)=7d-2
      !gen(20)%wa(nparam_per_node+2)=2.0d-1

      !gen(25)%wa(nparam_per_node+2)=0.4d-1     
      gen(21)%wa(nparam_per_node+2)=1.0d-1

      !gen(3)%wa(nparam_per_node+2)=5.0d-2

      gen(23)%wa(nparam_per_node+2)=-5.0d-2
      gen(23)%wa(8)=2.0d0
      gen(23)%wa(10)=2.0d0

      !gen(13)%wa(nparam_per_node+2)=2.0d-1/2d0
      !gen(14)%wa(nparam_per_node+2)=2.0d-1
      !gen(15)%wa(nparam_per_node+2)=5.0d-2

      
      !growth parameters for inward bud
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=3.2d-1
      !gen(21)%wa(nparam_per_node+2)=5d-2
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !growth parameters for outward bud
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=5.0d-3!5.0d-2
      !gen(21)%wa(nparam_per_node+2)=3.2d-1
      !gen(23)%wa(nparam_per_node+2)=-5.0d0
      !!!!!!!!!!!!!!!!!!!!!!!!!!   
      !growth parameters for rectangle scale
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=0.0d0
      !gen(21)%wa(nparam_per_node+2)=1.0d-3
      !gen(23)%wa(nparam_per_node+2)=-7d1
      !gen(6)%wa(nparam_per_node+2)=1.0d-1
      !!!!!!!!!!!!!!!!!!!!!!!!!!   
      
    !Gene-gene interactions

      !structural genes
      !auto-maintenance
      gen(1)%w(1)=1d0
      gen(2)%w(2)=1d0
      gen(3)%w(3)=1d0
      gen(4)%w(4)=1d0
      gen(5)%w(5)=1d0
      gen(6)%w(6)=1d0
      gen(7)%w(7)=1d0
      gen(11)%w(11)=1d0
      !maintenance of stable gene expression
      gen(13)%w(4)=1d0
      gen(14)%w(4)=1d0
      gen(15)%w(6)=1d0
      gen(8)%w(1)=1d0
      gen(9)%w(2)=1d0
      gen(10)%w(3)=1d0
      gen(22)%w(11)=1d0
      gen(8)%w(11)=-1d5
      gen(9)%w(11)=-1d5
      !!!!!!!!!!!!!!!!
      gen(12)%w(11)=1d0
      gen(13)%w(11)=-1d5
      gen(14)%w(11)=-1d5
      gen(14)%w(1)=-1d5
      gen(13)%w(2)=-1d5
      !gen(24)%w(11)=1d0 !s2
      !!!!!!!!!!!!!!!!! effector activation thresholds
      gen(19)%w(1)=-0.05d0
      gen(20)%w(2)=-0.05d0
      gen(21)%w(3)=-0.05d0
      !gen(19)%w(1)=0d0
      !gen(20)%w(2)=0d0
      !gen(21)%w(3)=0d0
      gen(19)%w(16)=5d1   !tooth
      !gen(20)%w(17)=5d1
      gen(21)%w(18)=5d1
      !gen(19)%w(16)=1d2    !scale
      !gen(20)%w(17)=1d2
      !gen(21)%w(18)=1d1
      !gen(19)%w(16)=1d1   !tooth
      !gen(20)%w(17)=1d1
      !gen(21)%w(18)=1d1
      !gen(19)%w(16)=1d2    !rectangle scale
      !gen(20)%w(17)=1d2
      !gen(21)%w(18)=1d1

      !gen(24)%w(19)=1d0 !s2
      !gen(25)%w(21)=1d0 !s2


      !effector supression of adhesion molecules
      !gen(8)%w(19)=-1d1
      !gen(9)%w(20)=-1d1
      !gen(10)%w(21)=-1d1

      !!!!!!!!!!!!!!!! a gene expressed in the posterior half of the placode supresses signalling
      gen(23)%w(23)=1d0
      !gen(13)%w(23)=-1d5
      !gen(14)%w(23)=-1d5

      !!!scale model 2
      !gen(24)%w(23)=1d0

      !gen(27)%w(1)=1d0
      !gen(28)%w(3)=1d0
      !gen(27)%w(23)=-1d5    
      !gen(28)%w(23)=-1d5


      !gen(25)%w(1)=-0.05d0
      !gen(26)%w(3)=-0.05d0
      !gen(25)%w(29)=5d1 
      !gen(26)%w(30)=5d1

      !gen(12)%w(25)=1d0


      !gen(27)%w(5)=-1d5
      !gen(28)%w(7)=-1d5

      !signal-receptor interactions
      !epithelial basal binding
      gen(16)%nww=4
      gen(16)%ww(1,1)=12
      gen(16)%ww(1,2)=16
      gen(16)%ww(1,3)=1d0
      gen(16)%ww(2,1)=16
      gen(16)%ww(2,2)=12
      gen(16)%ww(2,3)=1d0
      gen(16)%ww(3,1)=13
      gen(16)%ww(3,2)=16
      gen(16)%ww(3,3)=1d0
      gen(16)%ww(4,1)=16
      gen(16)%ww(4,2)=13
      gen(16)%ww(4,3)=1d0
      !epithelial suprabasal receptor binding
      gen(17)%nww=4
      gen(17)%ww(1,1)=12
      gen(17)%ww(1,2)=17
      gen(17)%ww(1,3)=1d0
      gen(17)%ww(2,1)=17
      gen(17)%ww(2,2)=12
      gen(17)%ww(2,3)=1d0
      gen(17)%ww(3,1)=14
      gen(17)%ww(3,2)=15
      gen(17)%ww(3,3)=1d0
      gen(17)%ww(4,1)=17
      gen(17)%ww(4,2)=14
      gen(17)%ww(4,3)=1d0
      !mesenchymal receptor binding
      gen(18)%nww=4
      gen(18)%ww(1,1)=12
      gen(18)%ww(1,2)=18
      gen(18)%ww(1,3)=1d0
      gen(18)%ww(2,1)=18
      gen(18)%ww(2,2)=12
      gen(18)%ww(2,3)=1d0
      gen(18)%ww(3,1)=15
      gen(18)%ww(3,2)=18
      gen(18)%ww(3,3)=1d0
      gen(18)%ww(4,1)=18
      gen(18)%ww(4,2)=15
      gen(18)%ww(4,3)=1d0
      


      !Wnt receptors
      gen(29)%nww=4
      gen(29)%ww(1,1)=24
      gen(29)%ww(1,2)=29
      gen(29)%ww(1,3)=1d0
      gen(29)%ww(2,1)=29
      gen(29)%ww(2,2)=24
      gen(29)%ww(2,3)=1d0
      gen(29)%ww(3,1)=27
      gen(29)%ww(3,2)=29
      gen(29)%ww(3,3)=1d0
      gen(29)%ww(4,1)=29
      gen(29)%ww(4,2)=27
      gen(29)%ww(4,3)=1d0
      !mesench receptor binding
      gen(30)%nww=4
      gen(30)%ww(1,1)=24
      gen(30)%ww(1,2)=30
      gen(30)%ww(1,3)=1d0
      gen(30)%ww(2,1)=30
      gen(30)%ww(2,2)=24
      gen(30)%ww(2,3)=1d0
      gen(30)%ww(3,1)=28
      gen(30)%ww(3,2)=30
      gen(30)%ww(3,3)=1d0
      gen(30)%ww(4,1)=30
      gen(30)%ww(4,2)=28
      gen(30)%ww(4,3)=1d0



    !Adhesion molecules

	ntipusadh=6
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      !as default, 1 and 2 should act as the same (basal and suprabasal)
      kadh(1,1)=1d0

      kadh(1,2)=0d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=0d0

      kadh(1,3)=5d0 ; kadh(3,1)=kadh(1,3)
      !rectangle scale
      !kadh(1,3)=3d0 ; kadh(3,1)=kadh(1,3)
      kadh(2,3)=0d0 ; kadh(3,2)=kadh(2,3)
      kadh(3,3)=1d0

      kadh(1,4)=5d-1 ; kadh(4,1)=kadh(1,4)
      kadh(2,4)=0d0 ; kadh(4,2)=kadh(2,4)
      kadh(3,4)=5d-1 ; kadh(4,3)=kadh(3,4)
      !kadh(4,4)=1d1 !tooth
      kadh(4,4)=1d1  !scale

      !kadh(1,5)=0d0 ; kadh(5,1)=kadh(1,5)
      !kadh(2,5)=0d0 ; kadh(5,2)=kadh(2,5)
      !kadh(3,5)=0d0 ; kadh(5,3)=kadh(3,5)
      !kadh(4,5)=0d0 ; kadh(5,4)=kadh(4,5)
      !kadh(5,5)=0d0

      !kadh(1,6)=0d0 ; kadh(6,1)=kadh(1,6)
      !kadh(2,6)=0d0 ; kadh(6,2)=kadh(2,6)
      !kadh(3,6)=0d0 ; kadh(6,3)=kadh(3,6)
      !kadh(4,6)=0d0 ; kadh(6,4)=kadh(4,6)
      !kadh(5,6)=1d1 ; kadh(6,5)=kadh(5,6)
      !kadh(6,6)=1d1

      !kadh(1,7)=5d-1 ; kadh(7,1)=kadh(1,7)
      !kadh(2,7)=0d0 ; kadh(7,2)=kadh(2,7)
      !kadh(3,7)=1d0 ; kadh(7,3)=kadh(3,7)
      !kadh(4,7)=5d-1 ; kadh(7,4)=kadh(4,7)
      !kadh(5,7)=5d-1 ; kadh(7,5)=kadh(5,7)
      !kadh(6,7)=0d0 ; kadh(7,6)=kadh(6,7)
      !kadh(7,7)=1d0

    end if

    !Gene expression on nodes

    
    !if(geometry==1)then
    !  gex(1,2)=1d0
    !  gex(ndepi+1,2)=1d0
    !  !gex(ndepi+1,2)=1d0
    !end if
    
    if(geometry==2)then
    !draw rectangles
    !!!!!!!!!!!!!!!!!!!!!!!!!
    jj1=1 ; jj2=ly
    ii1=1 ; ii2=lx
    !placode range
    px1=3 ; px2=lx-(px1-1)
    py1=3 ; py2=ly-(py1-1)
    !signalling centre
    sx1=6 ; sx2=lx-(sx1-1)
    sy1=5 ; sy2=ly-(sy1-1)
    !AP proliferation supressor
    apx1=6 ; apx2=lx
    apy1=1 ; apy2=ly

    do i=ii1,ii2
      do j=jj1,jj2
        !epithelium
        k=cels(cell_grid_epi(i,j))%node(1)

        gex(k,1)=1d0 ; gex(node(k)%altre,1)=1d0
        gex(k,8)=1d0 ; gex(node(k)%altre,8)=1d0
        !borders
        if(i<=1.or.i>=lx.or.j<=1.or.j>=ly)then
          node(k)%hold=2
          node(node(k)%altre)%hold=2
        end if
        !placode
        if((i>px1.and.i<px2).and.(j>py1.and.j<py2))then
          gex(k,4)=1d0 ; gex(node(k)%altre,4)=1d0
          gex(k,13)=1d0 ; gex(node(k)%altre,13)=1d0
        else
          gex(k,5)=1d0 ; gex(node(k)%altre,5)=1d0
        end if
        !signalling centre
        if((i>sx1.and.i<sx2).and.(j>sy1.and.j<sy2))then
          gex(k,11)=1d0 ; gex(node(k)%altre,11)=1d0
          gex(k,22)=1d0 ; gex(node(k)%altre,22)=1d0
          gex(k,13)=0d0 ; gex(node(k)%altre,13)=0d0
        end if
        !AP proliferation supressor
        if((i>apx1.and.i<apx2).and.(j>apy1.and.j<apy2))then
          gex(k,23)=1d0 ; gex(node(k)%altre,23)=1d0
        end if


        !mesench
        do ii=1,layer
          k=cels(cell_grid_mes(i,j,ii))%node(1)
          gex(k,3)=1d0
          gex(k,10)=1d0
          !borders
          if(i<=1.or.i>=lx.or.j<=1.or.j>=ly)then
            node(k)%hold=2
          end if
          !placode
          if((i>px1.and.i<px2).and.(j>py1.and.j<py2))then
            gex(k,6)=1d0
            gex(k,15)=1d0
          else
            gex(k,7)=1d0
          end if
        end do
      end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!
    end if

  if(geometry==1)then  

    ring=3 

    ring2=7

    !polarized=0 !tooth
    polarized=1 !scale

    do i=1,nd
      if(node(i)%tipus<3)then
        !if(node(i)%tipus==2)then
          !epithelial basal layer, basal side
          gex(i,1)=1d0  !receptors
          node(i)%marge=0
          d=sqrt(node(i)%x**2+node(i)%y**2)
          !setting expression on placode and interplacode
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,4)=1d0
            if(polarized==0.or.node(i)%x<=-0.1)then
              gex(i,13)=1d0 ; gex(node(i)%altre,13)=1d0
              gex(i,27)=1d0 ; gex(node(i)%altre,27)=1d0
            else
              gex(i,23)=1d0
              gex(i,13)=1d0 ; gex(node(i)%altre,13)=1d0
              gex(i,27)=0d0 ; gex(node(i)%altre,27)=0d0
            end if
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,5)=1d0
            node(i)%hold=3
          end if
          !!!!!!!!!!!!!
          !signalling centre
          if(d-2*node(i)%req*(ring-1)<epsilod)then
            node(i)%marge=0 ; node(node(i)%altre)%marge=0
            gex(i,11)=1d0 ; gex(node(i)%altre,11)=1d0
            gex(i,22)=1d0 ; gex(node(i)%altre,22)=1d0
            gex(i,13)=0d0 ; gex(node(i)%altre,13)=0d0
          end if

        !end if
        !epithelial basal layer, whole layer       
        gex(i,8)=1d0

      end if
      if(node(i)%tipus==3)then
        if(node(i)%z<=0)then
          !mesenchymal layer
          gex(i,3)=1d0
          gex(i,10)=1d0

          !setting expression on placode and interplacode
          d=sqrt(node(i)%x**2+node(i)%y**2)
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,6)=1d0
            gex(i,15)=1d0
            gex(i,28)=1d0
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,7)=1d0
            node(i)%hold=3
          end if
          !!!!!!!!!!!!!

        else
          !epithelial suprabasal layer
          gex(i,2)=1d0
          gex(i,9)=1d0

          !setting expression on placode and interplacode
          d=sqrt(node(i)%x**2+node(i)%y**2)
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,4)=1d0
            if(polarized==0.or.node(i)%x<=0.15)then
              gex(i,14)=1d0
            else
              gex(i,23)=1d0
              gex(i,14)=1d0
            end if
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,5)=1d0
          end if
          !!!!!!!!!!!!!
          !signalling centre
          !if(d-2*node(i)%req*(ring-1)<epsilod)then
          !  gex(i,11)=1d0
          !  gex(i,22)=1d0
          !  gex(i,14)=0d0 !; gex(node(i)%altre,13)=0d0
          !end if


        end if
      end if


    end do
    end if

    call update_npag

    node(:)%talone=0.0d0

    !ramax=maxval(node(:)%da)*3

    
end subroutine


subroutine dros_wing_vein_2d

!real*8::pericel,radio2,req1,req2,adjustreq,di
!integer::kk2

integer:: lx,ly
integer:: x1,x2,y1,y2,ii1,ii2,jj1,jj2
integer:: x3,x4,y3,y4,ii3,ii4,jj3,jj4,jj5,y5
integer, allocatable::cell_grid_epi(:,:),cell_grid_mes(:,:)
!integer:: periodic_borders,nborders
integer::geometry
integer::xlayer
integer::ring,ring2,polarized
integer::px1,px2,py1,py2,sx1,sx2,sy1,sy2,apx1,apx2,apy1,apy2
real*8::zepi2

!******* #1 DEFINING SPATIAL DIMENSIONS *******

  geometry=2 !geometry=1 means hexagonal epithelium ; geometry=2 means rectangular epithelium
  
  if (geometry==1)then
   
  elseif(geometry==2)then
    !layer=2      !number of planar cell layers
    lx=15
    ly=3

    zepi=0.75d0
    zepi2=-0.75d0
    
    allocate(cell_grid_epi(lx,ly),cell_grid_mes(lx,ly))

  end if  
    !periodic_borders=1
    


    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(geometry==1)then
 
    elseif(geometry==2)then

      ncelsepi=lx*ly*2
      ndepi=ncelsepi*2
      ncelsmes=0!ncelsepi*layer
      ndmes=0!ncelsmes
      nd=ndepi+ndmes
      ncels=ncelsepi+ncelsmes
      nodecel=2 ; nodecela=2
      nda=nd+10
      ncals=ncels+10
    end if
    

  !End initializing dimension parameters

  print*,"inicial nd",nd,"ncels",ncels,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=1
    !End Allocatations


    

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=1d0
    deltamin=1d-2
    khold=1d0
    angletor=0.05



    k_bu=5d0
    ramax=0.35d0
    k_press=5d-1
    m_xwall=0d0!zmes !tooth
    mi_xwall=0d0 !tooth

    mi_ywall=-0.65d0
    m_ywall=0.2165d0

    !m_xwall=0d0  !SCALE
    !mi_xwall=-5d-1 !SCALE

    ndmax=9d4
    
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=1
    ffu(2)=1 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !buoyancy

    ffu(11)=0 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(14)=1 !physical boundaries (walls)
    ffu(16)=0
    ffu(17)=0 !volume conservation (epithelial)
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise
    ffu(24)=0


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    !if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(mod(i,2)/=0)then  !basal
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
          !node(i)%reqp=-0.05d0
          !node(i)%req=0.20d0
        !node(i)%kplast=1d0
        else                      !apical
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !node(i)%kplast=1d1
        !  !node(i)%reqp=0.05d0
        !  !node(i)%req=0.30d0
        end if
        node(i)%req=0.25d0
        node(i)%reqcr=0.25d0 !; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.125
        node(i)%da=node(i)%req*2.0; 
        node(i)%ke=5d-2
        node(i)%tor=1d0
        node(i)%stor=5d0
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=1d1
        node(i)%kvol=5d1
      end do
    !end if


    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)

    node(:)%hold=0


  if(geometry==1)then


    
  elseif(geometry==2)then
	do i=1,ncels
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do


    !epithelial cells
    !origin of rectangle
    aa=node(1)%req*2
    bb=sqrt(aa**2-(aa*0.5)**2)
    a=-real(lx)*0.5*aa ; b=-real(ly)*0.5*bb
    !print*,"a",a,"b",b,"aa",aa,"bb",bb,"lx",lx,"ly",ly
    ii=0 ; jj=0
    !dorsal epithelial layer
    do i=1,ly
      do j=1,lx
        ii=ii+1 ;jj=jj+1 ; cell_grid_epi(j,i)=jj
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        end if
        !print*,"x",node(ii)%x,"y",node(ii)%y
        node(ii)%icel=jj ; node(ii)%altre=ii+1 ;node(ii)%tipus=2
        if(j==1.or.j==lx) node(ii)%hold=2
        ii=ii+1
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii-1 ;node(ii)%tipus=1
        node(ii)%hold=3
        cels(jj)%node(1)=ii-1 ; cels(jj)%node(2)=ii ;cels(jj)%nunodes=2; cels(jj)%ctipus=1
        !if(i==1.or.i==ly.or.j==1.or.j==lx)then
        !  node(ii)%hold=2;node(ii-1)%hold=2 ;
        !  node(ii)%repcel=1d1;node(ii-1)%repcel=1d1;
        !  !node(ii)%border=1;node(ii-1)%border=1;
        !end if
      end do
    end do

    !ventral epithelial layer
    do i=1,ly
      do j=1,lx
        ii=ii+1 ;jj=jj+1 ; cell_grid_mes(j,i)=jj
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi2
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi2
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii+1 ;node(ii)%tipus=2
        if(j==1.or.j==lx) node(ii)%hold=2
        ii=ii+1
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi2-2*node(ii)%reqs
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi2-2*node(ii)%reqs
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii-1 ;node(ii)%tipus=1
        node(ii)%hold=3
        cels(jj)%node(1)=ii-1 ; cels(jj)%node(2)=ii ;cels(jj)%nunodes=2; cels(jj)%ctipus=1
        !if(i==1.or.i==ly.or.j==1.or.j==lx)then
        !  node(ii)%hold=2;node(ii-1)%hold=2 ;
        !  node(ii)%repcel=1d1;node(ii-1)%repcel=1d1;
        !  !node(ii)%border=1;node(ii-1)%border=1;
        !end if
      end do
    end do

 
  end if

  !m_xwall=maxval(node(1:nd)%x)+0.01
  !mi_xwall=minval(node(1:nd)%x)-0.01
  !m_ywall=maxval(node(1:nd)%y)+0.01
  !mi_ywall=minval(node(1:nd)%y)-0.01
  
!print*,"tall2"

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
  !node(:)%hold=2
  !node(:)%hold=0  

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    !if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    !end if
    !Mesenchymal
    if(layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=4
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=0d0 ; gen(1)%mu=0d0 ;gen(1)%name="provein basal"
      gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=0d0 ;gen(2)%name="provein apical"
      gen(3)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=0d0 ;gen(3)%name="intervein basal"
      gen(4)%kindof=1 ; gen(4)%diffu=0d0 ; gen(4)%mu=0d0 ;gen(4)%name="intervein apical"

   

 
    !Gene-behavior interactions
    gen(1)%wa(1)=1
    gen(2)%wa(1)=2
    gen(3)%wa(1)=3
    gen(4)%wa(1)=4

    !gen(3)%wa(21)=-0.03d0
    !gen(4)%wa(21)=-0.03d0

    !Gene-gene interactions



    !Adhesion molecules

	ntipusadh=4
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      !as default, 1 and 2 should act as the same (basal and suprabasal)
      kadh(1,1)=0d0

      kadh(1,2)=0d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=1d0

      kadh(1,3)=1d0 ; kadh(3,1)=kadh(1,3)
      kadh(2,3)=0d0 ; kadh(3,2)=kadh(2,3)
      kadh(3,3)=5d0

      kadh(1,4)=0d0 ; kadh(4,1)=kadh(1,4)
      kadh(2,4)=0d0 ; kadh(4,2)=kadh(2,4)
      kadh(3,4)=0d0 ; kadh(4,3)=kadh(3,4)
      kadh(4,4)=1d0

    end if

    !Gene expression on nodes

    
    !if(geometry==1)then
    !  gex(1,2)=1d0
    !  gex(ndepi+1,2)=1d0
    !  !gex(ndepi+1,2)=1d0
    !end if
    
    if(geometry==2)then
    !draw rectangles
    !!!!!!!!!!!!!!!!!!!!!!!!!
    jj1=1 ; jj2=ly
    ii1=1 ; ii2=lx
    !placode range
    px1=5 ; px2=lx-(px1-1)
    py1=0 ; py2=ly+1!-(py1-1)
    !signalling centre
    sx1=7 ; sx2=lx-(sx1-1)
    sy1=5 ; sy2=ly-(sy1-1)
    !AP proliferation supressor
    apx1=7 ; apx2=lx
    apy1=1 ; apy2=ly

    do i=ii1,ii2
      do j=jj1,jj2
        !epithelium
        k=cels(cell_grid_epi(i,j))%node(1)

        !gex(k,1)=1d0
        !gex(k,8)=1d0
        !borders
        if(i==1.or.i==lx)then
          node(k)%z=0.25d0
        end if
        if(i==2.or.i==lx-1)then
          node(k)%z=0.35d0
        end if

        !placode
        if((i>px1.and.i<px2).and.(j>py1.and.j<py2))then
          gex(k,1)=1d0
          gex(node(k)%altre,2)=1d0
          !gex(k,4)=1d0
          !gex(k,13)=1d0
        else
          gex(k,3)=1d0
          gex(node(k)%altre,4)=1d0

          !gex(k,5)=1d0
        end if

        !ventral layer
        
        k=cels(cell_grid_mes(i,j))%node(1)
        !gex(k,3)=1d0
        !gex(k,10)=1d0
        !borders
        if(i==1.or.i==lx)then
          node(k)%z=-0.25d0
        end if
        if(i==2.or.i==lx-1)then
          node(k)%z=-0.35d0
        end if
        !placode
        if((i>px1.and.i<px2).and.(j>py1.and.j<py2))then
          gex(k,1)=1d0
          gex(node(k)%altre,2)=1d0
          !gex(k,6)=1d0
          !gex(k,15)=1d0
        else
          gex(k,3)=1d0
          gex(node(k)%altre,4)=1d0
          !gex(k,7)=1d0
        end if
      end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!
    end if



    call update_npag

    node(:)%talone=0.0d0

    !ramax=maxval(node(:)%da)*3

    
end subroutine




subroutine epi_folding_growth_test

!real*8::pericel,radio2,req1,req2,adjustreq,di
!integer::kk2

integer:: lx,ly
integer:: x1,x2,y1,y2,ii1,ii2,jj1,jj2
integer:: x3,x4,y3,y4,ii3,ii4,jj3,jj4,jj5,y5
integer, allocatable::cell_grid_epi(:,:),cell_grid_mes(:,:,:)
!integer:: periodic_borders,nborders
integer::geometry
integer::xlayer
integer::ring,ring2,polarized
integer::px1,px2,py1,py2,sx1,sx2,sy1,sy2,apx1,apx2,apy1,apy2

!******* #1 DEFINING SPATIAL DIMENSIONS *******

  geometry=2 !geometry=1 means hexagonal epithelium ; geometry=2 means rectangular epithelium
  
  if (geometry==1)then
	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=10    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=10   !number of radial cell layers
    layer=0      !number of planar cell layers
    zmes=0.75d0  !z-position of uppermost layer
    
    !allocate(mesradicel(layer))
    !do i=1,layer
    !  mesradicel(i)=mradicel
    !end do
    
    xlayer=0 !this will make a layer of ECM beneath the layers of mesenchymal cells
    
  elseif(geometry==2)then
    layer=0      !number of planar cell layers
    lx=16
    ly=32
    
    !allocate(cell_grid_epi(lx,ly),cell_grid_mes(lx,ly,layer))
    allocate(cell_grid_epi(lx,ly))

  end if  
    !periodic_borders=1
    


    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(geometry==1)then
      if(radi>0.and.radicel>0)then
        j=0
        do i=1,radi-1
          j=j+i
        end do
        nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
        nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13
  
        j=0
        do i=1,radicel-1
          j=j+i
        end do
        ncelsepi=(6*j+1)
        ndepi=nodecel*ncelsepi
      else
        ncelsepi=0
        ndepi=0
        nodecel=0
        nodecela=0
      end if
  
      if(mradi>0.and.mradicel>0.and.layer>0)then
        if(packed==1)then !if mesenchyme is packed, we use mradi differently
          j=0
          do i=1,mradi-1
            j=j+i
          end do
          nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
          nodecela=2*nodecel+1
        else
          nodecel=mradi
          nodecela=2*nodecel+1
        end if
  
        ncelsmes=0
        do k=1,layer
          j=0 !;print*,"mesradicel",mesradicel(k)
          do i=1,mesradicel(k)-1
            j=j+i
          end do
          ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
        end do
        ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
        if(radi==0.and.radicel==0)then
          nodecel=mradi
          nodecela=2*nodecel+1
        else if(nodecel<mradi)then
          nodecel=radi
          nodecela=nodecel*2+1
        end if
      else
        ndmes=0
        ncelsmes=0
      end if
      ndx=(6*j+1)*xlayer
      nd=ndepi+ndmes+ndx
      !ncels=ncelsepi+ncelsmes
      ncels=ncelsepi+ncelsmes+ndx
   print*,"ncels",ncels,ncelsepi,ncelsmes,ndx
      nda=nd+10
      ncals=ncels+10
    elseif(geometry==2)then

      ncelsepi=lx*ly
      ndepi=ncelsepi*2
      ncelsmes=ncelsepi*layer
      ndmes=ncelsmes
      nd=ndepi+ndmes
      ncels=ncelsepi+ncelsmes
      nodecel=2 ; nodecela=5
      nda=nd+10
      ncals=ncels+10
    end if
    

  !End initializing dimension parameters

  print*,"inicial nd",nd,"ncels",ncels,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=1
    !End Allocatations


    

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=0.85d0
    deltamin=1d-2
    khold=1d0
    angletor=0.05



    k_bu=5d0
    ramax=0.35d0
    k_press=5d-1
    m_xwall=zmes !tooth
    mi_xwall=0d0 !tooth

    !m_xwall=0d0  !SCALE
    !mi_xwall=-5d-1 !SCALE

    ndmax=9d4
    
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=1
    ffu(2)=1 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !buoyancy

    ffu(11)=0 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(14)=1 !physical boundaries (walls)
    ffu(17)=0 !volume conservation (epithelial)
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise
    ffu(24)=0


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    !if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(mod(i,2)/=0)then  !basal
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
          !node(i)%reqp=-0.05d0
          !node(i)%req=0.20d0
        !node(i)%kplast=1d0
        else                      !apical
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !node(i)%kplast=1d1
        !  !node(i)%reqp=0.05d0
        !  !node(i)%req=0.30d0
        end if
        node(i)%req=0.25d0
        node(i)%reqcr=0.25d0 !; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.125
        node(i)%da=node(i)%req*2.0; 
        node(i)%ke=5d0
        node(i)%tor=3d0
        node(i)%stor=1d0
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=1d1
        node(i)%kvol=5d1
      end do
    !end if

	!Mesenchyme
    if(layer>0)then
      do i=ndepi+1,ndepi+ndmes
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !scale rect
        !node(i)%rep=0d0 ; node(i)%repcel=1d0
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*2.0
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

	!Suprabasal cells
    if(xlayer>0)then
      do i=ndepi+ndmes+1,nd
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*2.0
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if
    
    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)

    node(:)%hold=0




  if(geometry==1)then
    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    !if(radi>0.and.radicel>0)               call epiteli_sphere1(radi,radicel,radius,di,zepi)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if
!print*,"tall0"

    !Mesenchymal layer
    !ii=0
    !do i=ndepi+ndmes+1,nd
    !  ii=ii+1
    !  j=i-ndx
    !  node(i)%x=node(j)%x ; node(i)%y=node(j)%y ; node(i)%z=-0.5d0
    !  node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    !  node(i)%icel=node(j)%icel+ndx ; node(i)%tipus=3
    !  !node(i)%icel=ncelsepi+ncelsmes+ii ; node(i)%tipus=3

    !  node(i)%repcel=5d-1
    !  !allocate(cels(node(i)%icel)%node(nodela))
    !  cels(node(i)%icel)%nunodes=1
    !  cels(node(i)%icel)%node(1)=i
    !  node(i)%marge=0
    !  !print*,"i",i,"j",j,"icel i",node(i)%icel,"icel j",node(j)%icel
    !  !node(i)%border=1
    !end do
    !ncelsmes=ncelsmes+ndx ; ndx=0
    !ndmes=ncelsmes
!print*,"tall1"


    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
      
    j=0
    do i=1,radicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    do i=j+1,ndepi
    !print*,"hola"
      node(i)%hold=1; node(i)%oriz=node(i)%z-5.0d0 !;node(i)%repcel=1d0
    end do
    !j=0
    !do i=1,mradicel-2
    !  j=j+i
    !end do
    !j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    !k=0
    !do i=1,mradicel-1
    !  k=k+i
    !end do
    !k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    !do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
    !  node(i)%hold=2;node(i)%repcel=1d0
    !  !node(i)%border=1
    !  !node(i)%da=node(i)%req*2.0
    !  !node(i)%orix=0 ; node(i)%oriy=0
    !  !node(i)%rep=1d1;node(i)%repcel=1d1
    !end do
    !do i=ndepi+ndmes/2+j+1,nd
    !  node(i)%hold=2;node(i)%repcel=1d0
    !  !node(i)%border=1
    !end do
    !do i=ndepi+k+1,nd
    !  node(i)%hold=2;node(i)%repcel=1d0
    !  node(i)%border=1
    !end do

    
  elseif(geometry==2)then
	do i=1,ncels
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do


    !epithelial cells
    !origin of rectangle
    aa=node(1)%req*2
    bb=sqrt(aa**2-(aa*0.5)**2)
    a=-real(lx)*0.5*aa ; b=-real(ly)*0.5*bb
    !print*,"a",a,"b",b,"aa",aa,"bb",bb,"lx",lx,"ly",ly
    ii=0 ; jj=0
    do i=1,ly
      do j=1,lx
        ii=ii+1 ;jj=jj+1 ; cell_grid_epi(j,i)=jj
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
          node(ii)%orix=node(ii)%x ; node(ii)%oriy=node(ii)%y ; node(ii)%oriz=node(ii)%z
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
          node(ii)%orix=node(ii)%x ; node(ii)%oriy=node(ii)%y ; node(ii)%oriz=node(ii)%z
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii+1 ;node(ii)%tipus=2
        ii=ii+1
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
          node(ii)%orix=node(ii)%x ; node(ii)%oriy=node(ii)%y ; node(ii)%oriz=node(ii)%z
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
          node(ii)%orix=node(ii)%x ; node(ii)%oriy=node(ii)%y ; node(ii)%oriz=node(ii)%z
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii-1 ;node(ii)%tipus=1
        cels(jj)%node(1)=ii-1 ; cels(jj)%node(2)=ii ;cels(jj)%nunodes=2; cels(jj)%ctipus=1
        if(i==1.or.i==ly.or.j==1.or.j==lx)then
          node(ii)%hold=1;node(ii-1)%hold=1 ;
          node(ii)%repcel=1d0;node(ii-1)%repcel=1d0;
          node(ii)%oriz=node(ii)%z-5.0d0  ;  node(ii-1)%oriz=node(ii-1)%z-5.0d0

          !node(ii)%border=1;node(ii-1)%border=1;
        end if
      end do
    end do

    !mesenchymals
    !do k=1,layer
    !  !print*,"k",k
    !  do i=1,ly
    !    do j=1,lx
    !      ii=ii+1 ; jj=jj+1 ; cell_grid_mes(j,i,k)=jj
    !      if(mod(i,2)==0)then
    !        node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
    !      else
    !        node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
    !      end if
    !      node(ii)%icel=jj ; node(ii)%altre=0 ;node(ii)%tipus=3
    !      cels(jj)%node(1)=ii ;cels(jj)%nunodes=1 ; cels(jj)%ctipus=3
    !      if(i==1 .or. i==ly .or. j==1 .or. j==lx)then
    !        node(ii)%hold=2 ;node(ii)%repcel=1d1 ; !node(ii)%border=2 ;
    !      end if
    !      !if(i==1 .or. i==ly .or. j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
    !      !if(j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
    !    end do
    !  end do
    !end do

    !do i=1,nd
    !  node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    !end do


  end if

  !m_xwall=maxval(node(1:nd)%x)+0.01
  !mi_xwall=minval(node(1:nd)%x)-0.01
  !m_ywall=maxval(node(1:nd)%y)+0.01
  !mi_ywall=minval(node(1:nd)%y)-0.01
  
!print*,"tall2"

    !do i=1,nd
    !  node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    !end do
    
  !node(:)%hold=2
  !node(:)%hold=0  

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    !if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    !end if
    !Mesenchymal
    if(layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=2
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=0d0 ; gen(1)%mu=5d-1 ;gen(1)%name="Stalk cells marker"
      gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=5d-1 ;gen(2)%name="Growing cells marker"
      
    


 
    !Gene-behavior interactions

      gen(1)%wa(1)=1       
      gen(2)%wa(1)=2   


      gen(1)%wa(nparam_per_node+2)=-1d5       
      !gen(2)%wa(nparam_per_node+2)=5d-2   
      
    !Gene-gene interactions

      !structural genes
      !auto-maintenance
      gen(1)%w(1)=1d0
      gen(2)%w(2)=1d0


    !Adhesion molecules

	ntipusadh=2
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      !as default, 1 and 2 should act as the same (basal and suprabasal)
      kadh(1,1)=5d0

      kadh(1,2)=5d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=5d0

    end if

    !Gene expression on nodes

    
    !if(geometry==1)then
    !  gex(1,2)=1d0
    !  gex(ndepi+1,2)=1d0
    !  !gex(ndepi+1,2)=1d0
    !end if
    
    if(geometry==2)then
    !draw rectangles
    !!!!!!!!!!!!!!!!!!!!!!!!!
    jj1=1 ; jj2=ly
    ii1=1 ; ii2=lx
    !placode range
    px1=3 ; px2=lx-(px1-1)
    py1=3 ; py2=ly-(py1-1)
    !signalling centre
    sx1=7 ; sx2=lx-(sx1-1)
    sy1=6 ; sy2=ly-(sy1-1)
    !AP proliferation supressor
    apx1=7 ; apx2=lx
    apy1=1 ; apy2=ly

    do i=ii1,ii2
      do j=jj1,jj2
        !epithelium
        k=cels(cell_grid_epi(i,j))%node(1)

        !gex(k,1)=1d0
        !gex(k,8)=1d0
        !borders
        !if(i<=1.or.i>=lx.or.j<=1.or.j>=ly)then
        !  node(k)%hold=3
        !  node(node(k)%altre)%hold=3
        !end if
        !placode
        if((i>sx1.and.i<sx2).and.(j>sy1.and.j<sy2))then
          gex(k,1)=1d0 ; gex(node(k)%altre,1)=1d0 ; node(k)%hold=2 ; node(node(k)%altre)%hold=2
        else
          gex(k,2)=1d0 ; gex(node(k)%altre,2)=1d0
        end if
        !signalling centre
        !if((i>sx1.and.i<sx2).and.(j>sy1.and.j<sy2))then
        !  gex(k,11)=1d0 ; gex(node(k)%altre,11)=1d0
        !  gex(k,22)=1d0 ; gex(node(k)%altre,22)=1d0
        !  gex(k,13)=0d0
        !end if
        !!AP proliferation supressor
        !if((i>apx1.and.i<apx2).and.(j>apy1.and.j<apy2))then
        !  gex(k,23)=1d0
        !end if


        !mesench
        !do ii=1,layer
        !  k=cels(cell_grid_mes(i,j,ii))%node(1)
        !  gex(k,3)=1d0
        !  gex(k,10)=1d0
        !  !borders
        !  if(i<=1.or.i>=lx.or.j<=1.or.j>=ly)then
        !    node(k)%hold=3
        !  end if
        !  !placode
        !  if((i>px1.and.i<px2).and.(j>py1.and.j<py2))then
        !    gex(k,6)=1d0
        !    gex(k,15)=1d0
        !  else
        !    gex(k,7)=1d0
        !  end if
        !end do
      end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!
    end if

  if(geometry==1)then  
    !ring=2  !scale
    ring=2 !tooth

    ring2=3

    polarized=0 !tooth
    !polarized=1 !scale

    do i=1,nd
      if(node(i)%tipus<3)then
        !if(node(i)%tipus==2)then
          !epithelial basal layer, basal side
          !gex(i,1)=1d0  !receptors
          node(i)%marge=0
          d=sqrt(node(i)%x**2+node(i)%y**2)
          !setting expression on placode and interplacode
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
          print*,"hola"
            gex(i,1)=1d0 ; node(i)%hold=2 
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
          print*,"hello"
            gex(i,2)=1d0 
          end if
          !!!!!!!!!!!!!
          !signalling centre
          !if(d-2*node(i)%req*(ring-1)<epsilod)then
          !  node(i)%marge=0 !; node(node(i)%altre)%marge=0
          !  gex(i,11)=1d0 !; gex(node(i)%altre,11)=1d0
          !  gex(i,22)=1d0 !; gex(node(i)%altre,22)=1d0
          !  gex(i,13)=0d0 !; gex(node(i)%altre,13)=0d0
          !end if

        !end if
        !epithelial basal layer, whole layer       
        !gex(i,8)=1d0

      end if

    end do
  end if

    call update_npag

    node(:)%talone=0.0d0

    !ramax=maxval(node(:)%da)*3

    
end subroutine





subroutine epidermal_organ_in_tooth_bud_s2_small_2D

!real*8::pericel,radio2,req1,req2,adjustreq,di
!integer::kk2

integer:: lx,ly
integer:: x1,x2,y1,y2,ii1,ii2,jj1,jj2
integer:: x3,x4,y3,y4,ii3,ii4,jj3,jj4,jj5,y5
integer, allocatable::cell_grid_epi(:,:),cell_grid_mes(:,:,:),cell_grid_supra(:,:,:)
!integer:: periodic_borders,nborders
integer::geometry
integer::xlayer
integer::ring,ring2,polarized
integer::px1,px2,py1,py2,sx1,sx2,sy1,sy2,apx1,apx2,apy1,apy2

!******* #1 DEFINING SPATIAL DIMENSIONS *******

  geometry=2 !geometry=1 means hexagonal epithelium ; geometry=2 means rectangular epithelium
  
  if (geometry==1)then
	!epithelium's dimension parameters
    radi=1       !number of radial layers of nodes per cell
    radicel=5    !number of radial layers of cells
    zepi=0d0     !z-position of basal layer

    !mesenchyme's dimension parameters
    packed=1     !this tells if we want the mesenchymal cells packed (hexagonally) or just loose
    mradi=1      !number of nodes per cell    !if packed this is number of radial nodes per cell
    mradicel=5   !number of radial cell layers
    layer=1      !number of planar cell layers
    zmes=0.75d0  !z-position of uppermost layer
    
    allocate(mesradicel(layer))
    do i=1,layer
      mesradicel(i)=mradicel
    end do
    
    xlayer=1 !this will make a layer of ECM beneath the layers of mesenchymal cells
    
  elseif(geometry==2)then
    layer=1      !number of planar cell layers
    lx=10
    ly=1

    xlayer=2 !this will make a layer of ECM beneath the layers of mesenchymal cells

    allocate(cell_grid_epi(lx,ly),cell_grid_mes(lx,ly,layer))
    if(xlayer>0) allocate(cell_grid_supra(lx,ly,xlayer))



  end if  
    !periodic_borders=1
    


    !ECM dimension parameters?

!************************************************


    !Initializing dimension parameters
    if(geometry==1)then
      if(radi>0.and.radicel>0)then
        j=0
        do i=1,radi-1
          j=j+i
        end do
        nodecel=(6*j+1)*2	;if(radi==1) nodecel=2!number of nodes per cell
        nodecela=2*nodecel+1	!cels()%node array's initial size the same for mesenchyme!>>>>>>>Miquel 21-3-13
  
        j=0
        do i=1,radicel-1
          j=j+i
        end do
        ncelsepi=(6*j+1)
        ndepi=nodecel*ncelsepi
      else
        ncelsepi=0
        ndepi=0
        nodecel=0
        nodecela=0
      end if
  
      if(mradi>0.and.mradicel>0.and.layer>0)then
        if(packed==1)then !if mesenchyme is packed, we use mradi differently
          j=0
          do i=1,mradi-1
            j=j+i
          end do
          nodecel=(6*j+1)*mradi  ;if(mradi==1) nodecel=1  !number of nodes per cell
          nodecela=2*nodecel+1
        else
          nodecel=mradi
          nodecela=2*nodecel+1
        end if
  
        ncelsmes=0
        do k=1,layer
          j=0 !;print*,"mesradicel",mesradicel(k)
          do i=1,mesradicel(k)-1
            j=j+i
          end do
          ncelsmes=ncelsmes+(6*j+1) ;print*,"j",j !number of mesenchymal cells
        end do
        ndmes=nodecel*ncelsmes    !number of mesenchymal nodes 
        if(radi==0.and.radicel==0)then
          nodecel=mradi
          nodecela=2*nodecel+1
        else if(nodecel<mradi)then
          nodecel=radi
          nodecela=nodecel*2+1
        end if
      else
        ndmes=0
        ncelsmes=0
      end if
      ndx=(6*j+1)*xlayer
      nd=ndepi+ndmes+ndx
      !ncels=ncelsepi+ncelsmes
      ncels=ncelsepi+ncelsmes+ndx
   print*,"ncels",ncels,ncelsepi,ncelsmes,ndx
      nda=nd+10
      ncals=ncels+10
    elseif(geometry==2)then

      ncelsepi=lx*ly
      ndepi=ncelsepi*2
      ncelsmes=ncelsepi*layer
      ndmes=ncelsmes
      nd=ndepi+ndmes+ndmes*xlayer
      ncels=ncelsepi+ncelsmes+ncelsmes*xlayer
      nodecel=2 ; nodecela=5
      nda=nd+10
      ncals=ncels+10
    end if
    

  !End initializing dimension parameters

  print*,"inicial nd",nd,"ncels",ncels,"nodecel",nodecel,"ncelsepi",ncelsepi,"ncelsmes",ncelsmes,"ndepi",ndepi,"ndmes",ndmes


    !Allocatations
    if (allocated(node)) deallocate(node)
    if (allocated(cels)) deallocate(cels)
    allocate(node(nda),cels(ncals)) 
    call iniarrays
    if(radi==1.or.mradi==1) ffu(1)=1
    !End Allocatations


    

   !******* #2 DEFINING MODEL PARAMETERS *******

    !implementation and initializations
    getot=0
    itacc=0
    nparti=1000
    idum=-11111
    idumoriginal=idum
    realtime=0

    !nparam=32       !number of params
    nvarglobal_out=5

    !physical
    temp=1d0 !low value is low temperature	
    desmax=0.01
    resmax=1d-3
    prop_noise=0.0d0
    deltamax=1d-2 ! miguel 14-10-13   
    dmax=1
    screen_radius=0.85d0
    deltamin=1d-2
    khold=1d0
    angletor=0.05



    k_bu=5d0
    ramax=0.35d0
    k_press=0.0d0!5d-1
    m_xwall=0.0d0!0.75 !tooth
    mi_xwall=0d0 !tooth

    !m_xwall=0d0  !SCALE
    !mi_xwall=-5d-1 !SCALE

    ndmax=9d4
    
    !biological
    !functions used
    
    
    ffu=0
     !spring of the ellipse
    ffu(1)=1
    ffu(2)=1 !to quit if there are too many cells
    ffu(3)=1 !screening
    ffu(4)=0 !torsion
    ffu(5)=0 !external signal source
    ffu(6)=0 !buoyancy

    ffu(11)=0 !epithelial node plastic deformation
    ffu(12)=1 !0 dynamic delta, 1 fixed delta
    ffu(13)=0 !neighboring by triangulation
    ffu(14)=1 !physical boundaries (walls)
    ffu(17)=0 !volume conservation (epithelial)
    ffu(22)=1 !0 = noise biased by energies , 1 = unbiased noise
    ffu(24)=0


   !******* #2 DEFINING NODE MECHANIC PARAMETERS *******
    !Epithelium
    !if(radi>0.and.radicel>0)then
      do i=1,ndepi
        if(mod(i,2)/=0)then  !basal
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
          !node(i)%reqp=-0.05d0
          !node(i)%req=0.20d0
        !node(i)%kplast=1d0
        else                      !apical
          node(i)%you=0d0 ; node(i)%adh=0d0
          node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !node(i)%kplast=1d1
        !  !node(i)%reqp=0.05d0
        !  !node(i)%req=0.30d0
        end if
        node(i)%req=0.25d0
        node(i)%reqcr=0.25d0 !; node(i)%reqcr=node(i)%req
        node(i)%reqs=0.125
        node(i)%da=node(i)%req*2.0; 
        node(i)%ke=5d0
        node(i)%tor=3d0
        node(i)%stor=1d0
        node(i)%mo=temp
        node(i)%dmo=0
        node(i)%diffe=0d0
        node(i)%khold=khold
        node(i)%kplast=1d1
        node(i)%kvol=5d1
      end do
    !end if

	!Mesenchyme
    if(layer>0)then
      do i=ndepi+1,ndepi+ndmes
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        !scale rect
        !node(i)%rep=0d0 ; node(i)%repcel=1d0
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*2.0
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if

	!Suprabasal cells
    if(xlayer>0)then
      do i=ndepi+ndmes+1,nd
        node(i)%you=0d0 ; node(i)%adh=0d0
        node(i)%rep=0d0 ; node(i)%repcel=5d-1
        node(i)%req=0.25d0 ; node(i)%reqcr=node(i)%req
        node(i)%reqs=0d0
        node(i)%da=node(i)%req*2.0
        node(i)%ke=0d0   !only for epithelium
        node(i)%tor=0d0  !only for epithelium
        node(i)%stor=0d0 !only for epithelium
        node(i)%altre=0  !only for epithelium
        node(i)%mo=temp
        node(i)%dmo=0.0
        node(i)%diffe=0d0
        node(i)%khold=khold
      end do
    end if
    
    !General parameters that depend on node parameters
    rv=2*maxval(node(:)%da)

    node(:)%hold=0


  if(geometry==1)then
    !Distribution of nodes in space
    if(radi>0.and.radicel>0)               call epiteli(radi,radicel,zepi,node(1)%req,node(1)%reqs)
    !if(radi>0.and.radicel>0)               call epiteli_sphere1(radi,radicel,radius,di,zepi)
    if(mradi>0.and.mradicel>0.and.layer>0)then
                       if(packed==0)then ; call mesenq(mradi,mradicel,layer,zmes,node(ndepi+1)%req)
                       else  ; call mesenq_packed(mradi,mradicel,layer,zmes,node(ndepi+1)%req) ; end if ; end if
!print*,"tall0"

    !Mesenchymal layer
    ii=0
    do i=ndepi+ndmes+1,nd
      ii=ii+1
      j=i-ndx
      node(i)%x=node(j)%x ; node(i)%y=node(j)%y ; node(i)%z=-0.5d0
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
      node(i)%icel=node(j)%icel+ndx ; node(i)%tipus=3
      !node(i)%icel=ncelsepi+ncelsmes+ii ; node(i)%tipus=3

      node(i)%repcel=5d-1
      !allocate(cels(node(i)%icel)%node(nodela))
      cels(node(i)%icel)%nunodes=1
      cels(node(i)%icel)%node(1)=i
      node(i)%marge=0
      !print*,"i",i,"j",j,"icel i",node(i)%icel,"icel j",node(j)%icel
      !node(i)%border=1
    end do
    ncelsmes=ncelsmes+ndx ; ndx=0
    ndmes=ncelsmes
!print*,"tall1"
      
    j=0
    do i=1,radicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(1)%nunodes!nodecel !this is the number of epithelial nodes wich are not external
    do i=j+1,ndepi
      node(i)%hold=2 ;node(i)%repcel=1d0
    end do
    j=0
    do i=1,mradicel-2
      j=j+i
    end do
    j=(6*j+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    k=0
    do i=1,mradicel-1
      k=k+i
    end do
    k=(6*k+1)*cels(ncelsepi+1)%nunodes!mradi !this is the number of mesenchymal nodes wich are not external
    do i=ndepi+j+1,ndepi+k  !upper mesenchymal layer
      node(i)%hold=2;node(i)%repcel=1d0
      !node(i)%border=1
      !node(i)%da=node(i)%req*2.0
      !node(i)%orix=0 ; node(i)%oriy=0
      !node(i)%rep=1d1;node(i)%repcel=1d1
    end do
    do i=ndepi+ndmes/2+j+1,nd
      node(i)%hold=2;node(i)%repcel=1d0
      !node(i)%border=1
    end do
    !do i=ndepi+k+1,nd
    !  node(i)%hold=2;node(i)%repcel=1d0
    !  node(i)%border=1
    !end do

    
  elseif(geometry==2)then
	do i=1,ncels
		cels(i)%nunodes=0
		allocate(cels(i)%node(nodecela))		!>>>>>>>>>>>>Miquel 21-3-13
		cels(i)%nodela=nodecela
	end do


    !epithelial cells
    !origin of rectangle
    aa=node(1)%req*2
    bb=sqrt(aa**2-(aa*0.5)**2)
    a=-real(lx)*0.5*aa ; b=0.0d0!b=-real(ly)*0.5*bb
    !print*,"a",a,"b",b,"aa",aa,"bb",bb,"lx",lx,"ly",ly
    ii=0 ; jj=0
    !do i=1,ly
      i=1
      do j=1,lx
        ii=ii+1 ;jj=jj+1 ; cell_grid_epi(j,i)=jj
        print*,"jj1",jj
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii+1 ;node(ii)%tipus=2
        ii=ii+1
        if(mod(i,2)==0)then
          node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        else
          node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+2*node(ii)%reqs
        end if
        node(ii)%icel=jj ; node(ii)%altre=ii-1 ;node(ii)%tipus=1
        cels(jj)%node(1)=ii-1 ; cels(jj)%node(2)=ii ;cels(jj)%nunodes=2; cels(jj)%ctipus=1
        if(j==1.or.j==lx)then
          node(ii)%hold=2;node(ii-1)%hold=2 ;
          node(ii)%repcel=1d1;node(ii-1)%repcel=1d1;
          !node(ii)%border=1;node(ii-1)%border=1;
        end if
      end do
    !end do

    !mesenchymals
    do k=1,layer
      !print*,"k",k
      !do i=1,ly
        i=1
        do j=1,lx
          ii=ii+1 ; jj=jj+1 ; cell_grid_mes(j,i,k)=jj
        print*,"jj2",jj
          if(mod(i,2)==0)then
            node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          else
            node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi-aa*k
          end if
          node(ii)%icel=jj ; node(ii)%altre=0 ;node(ii)%tipus=3
          cels(jj)%node(1)=ii ;cels(jj)%nunodes=1 ; cels(jj)%ctipus=3
          if(j==1 .or. j==lx)then
            node(ii)%hold=2 ;node(ii)%repcel=1d1 ; !node(ii)%border=2 ;
          end if
          !if(i==1 .or. i==ly .or. j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
          !if(j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
        end do
      !end do
    end do

    !suprabasal layer
    if(xlayer>0)then
    do k=1,xlayer
      !print*,"k",k
      !do i=1,ly
        i=1
        do j=1,lx
          ii=ii+1 ; jj=jj+1 ; cell_grid_supra(j,i,k)=jj
        print*,"jj3",jj
          if(mod(i,2)==0)then
            node(ii)%x=a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+0.25+aa*k
          else
            node(ii)%x=0.5*aa+a+(j-1)*aa ; node(ii)%y=b+(i-1)*bb ; node(ii)%z=zepi+0.25+aa*k
          end if
          node(ii)%icel=jj ; node(ii)%altre=0 ;node(ii)%tipus=3
          cels(jj)%node(1)=ii ;cels(jj)%nunodes=1 ; cels(jj)%ctipus=3
          if(j==1 .or. j==lx .or. k==2)then
            node(ii)%hold=2 ;node(ii)%repcel=1d1 ; !node(ii)%border=2 ;
          end if
          !if(i==1 .or. i==ly .or. j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
          !if(j==1 .or. j==lx)then ;node(ii)%hold=1 ;node(ii)%repcel=1d2 ;end if
        end do
      !end do
    end do
    end if

  end if

  !m_xwall=maxval(node(1:nd)%x)+0.01
  !mi_xwall=minval(node(1:nd)%x)-0.01
  !m_ywall=maxval(node(1:nd)%y)+0.01
  !mi_ywall=minval(node(1:nd)%y)-0.01
  
!print*,"tall2"

    do i=1,nd
      node(i)%orix=node(i)%x ; node(i)%oriy=node(i)%y ; node(i)%oriz=node(i)%z
    end do
    
  !node(:)%hold=2
  !node(:)%hold=0  

    !!end of setting boundary nodes


   !******* #3 DEFINING CELL PARAMETERS *******
    !Epithelial
    !if(radi>0.and.radicel>0)then
      do i=1,ncelsepi
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    !end if
    !Mesenchymal
    if(layer>0)then
      do i=ncelsepi+1,ncels
        cels(i)%fase=0d0
        mmae=node(cels(i)%node(1))%req
        cels(i)%minsize_for_div=cels(i)%nunodes*2
        cels(i)%maxsize_for_div=10000 !>>> Is 5-2-14
      end do
    end if

    do i=1,ncels !cells start with at random points of the cell cycle, so the divisions are not synchronous
      !if(i<=7.or.(i>ncelsepi.and.i<=ncelsepi+7))then
        call random_number(a)
        cels(i)%fase=a
      !end if
    end do
    
   !******* #4 DEFINING GENETIC PARAMETERS *******

    !Number of genes
    ng=31
    call initiate_gene
    !Gene parameters
      !gen(:)%kindof= ; gen(:)%diffu= ; gen(:)%mu=

      gen(1)%kindof=1 ; gen(1)%diffu=0d0 ; gen(1)%mu=5d-1 ;gen(1)%name="epithelial basal marker"
      gen(2)%kindof=1 ; gen(2)%diffu=0d0 ; gen(2)%mu=5d-1 ;gen(2)%name="epithelial suprabasal marker"
      gen(3)%kindof=1 ; gen(3)%diffu=0d0 ; gen(3)%mu=5d-1 ;gen(3)%name="mesenchymal marker"
      gen(4)%kindof=1 ; gen(4)%diffu=0d0 ; gen(4)%mu=5d-1 ;gen(4)%name="epi. placode marker"
      gen(5)%kindof=1 ; gen(5)%diffu=0d0 ; gen(5)%mu=5d-1 ;gen(5)%name="epi. interplacode marker"
      gen(6)%kindof=1 ; gen(6)%diffu=0d0 ; gen(6)%mu=5d-1 ;gen(6)%name="supra. placode marker"
      gen(7)%kindof=1 ; gen(7)%diffu=0d0 ; gen(7)%mu=5d-1 ;gen(7)%name="supra. interplacode marker"

      gen(8)%kindof=1 ; gen(8)%diffu=0d0 ; gen(8)%mu=5d-1 ;gen(8)%name="epi. basal adhesion molecule"
      gen(9)%kindof=1 ; gen(9)%diffu=0d0 ; gen(9)%mu=5d-1 ;gen(9)%name="epi. suprabasal adhesion molecule"
      gen(10)%kindof=1 ; gen(10)%diffu=0d0 ; gen(10)%mu=5d-1 ;gen(10)%name="mesench. adhesion molecule"

      gen(11)%kindof=1 ; gen(11)%diffu=0d0 ; gen(11)%mu=5d-1 ;gen(11)%name="epi. signalling centre marker"
      gen(12)%kindof=4 ; gen(12)%diffu=5d-1 ; gen(12)%mu=5d-1 ;gen(12)%name="epi. signal"

      gen(13)%kindof=2 ; gen(13)%diffu=0d0 ; gen(13)%mu=5d-1 ;gen(13)%name="epi. basal receptor inactive"
      gen(14)%kindof=2 ; gen(14)%diffu=0d0 ; gen(14)%mu=5d-1 ;gen(14)%name="epi. suprabasal receptor inactive"
      gen(15)%kindof=2 ; gen(15)%diffu=0d0 ; gen(15)%mu=5d-1 ;gen(15)%name="mesench. receptor inactive"

      gen(16)%kindof=8 ; gen(16)%diffu=0d0 ; gen(16)%mu=5d-1 ;gen(16)%name="epi. basal receptor active"
      gen(16)%npre=2 ; allocate(gen(16)%pre(gen(16)%npre)) ; gen(16)%pre(1)=12 ; gen(16)%pre(2)=13
      gen(16)%npost=2 ; allocate(gen(16)%post(gen(16)%npost)) ; gen(16)%post(1)=12 ; gen(16)%post(2)=13

      gen(17)%kindof=8 ; gen(17)%diffu=0d0 ; gen(17)%mu=5d-1 ;gen(17)%name="epi. suprabasal receptor active"
      gen(17)%npre=2 ; allocate(gen(17)%pre(gen(17)%npre)) ; gen(17)%pre(1)=12 ; gen(17)%pre(2)=14
      gen(17)%npost=2 ; allocate(gen(17)%post(gen(17)%npost)) ; gen(17)%post(1)=12 ; gen(17)%post(2)=14

      gen(18)%kindof=8 ; gen(18)%diffu=0d0 ; gen(18)%mu=5d-1 ;gen(18)%name="mesench. receptor active"
      gen(18)%npre=2 ; allocate(gen(18)%pre(gen(18)%npre)) ; gen(18)%pre(1)=12 ; gen(18)%pre(2)=15
      gen(18)%npost=2 ; allocate(gen(18)%post(gen(18)%npost)) ; gen(18)%post(1)=12 ; gen(18)%post(2)=15

      gen(19)%kindof=1 ; gen(19)%diffu=0d0 ; gen(19)%mu=5d-1 ;gen(19)%name="epi. basal effector"
      gen(20)%kindof=1 ; gen(20)%diffu=0d0 ; gen(20)%mu=5d-1 ;gen(20)%name="epi. suprabasal effector"
      gen(21)%kindof=1 ; gen(21)%diffu=0d0 ; gen(21)%mu=5d-1 ;gen(21)%name="mesenchymal effector"

      gen(22)%kindof=1 ; gen(22)%diffu=0d0 ; gen(22)%mu=5d-1 ;gen(22)%name="signalling centre adh. mol."
      
      gen(23)%kindof=1 ; gen(23)%diffu=0d0 ; gen(23)%mu=5d-1 ;gen(23)%name="s2 IEE adhesion molecule."
      gen(24)%kindof=1 ; gen(24)%diffu=0d0 ; gen(24)%mu=5d-1 ;gen(24)%name="s2 dermal papilla adhesion molecule"

      gen(25)%kindof=4 ; gen(25)%diffu=5d-1 ; gen(25)%mu=5d-1 ;gen(25)%name="s2 secondary signal"
      gen(26)%kindof=1 ; gen(26)%diffu=0d0 ; gen(26)%mu=5d-1 ;gen(26)%name="s2 enamel knot marker basal"
      gen(27)%kindof=1 ; gen(27)%diffu=0d0 ; gen(27)%mu=5d-1 ;gen(27)%name="s2 enamel knot marker supra"

      gen(28)%kindof=2 ; gen(28)%diffu=0d0 ; gen(28)%mu=5d-1 ;gen(28)%name="s2 epi. knot basal receptor inactive"
      gen(29)%kindof=2 ; gen(29)%diffu=0d0 ; gen(29)%mu=5d-1 ;gen(29)%name="s2 epi. knot supra receptor inactive"

      gen(30)%kindof=8 ; gen(30)%diffu=0d0 ; gen(30)%mu=5d-1 ;gen(30)%name="s2 epi. knot basal receptor active"
      gen(30)%npre=2 ; allocate(gen(30)%pre(gen(30)%npre)) ; gen(30)%pre(1)=25 ; gen(30)%pre(2)=28
      gen(30)%npost=2 ; allocate(gen(30)%post(gen(30)%npost)) ; gen(30)%post(1)=25 ; gen(30)%post(2)=28

      gen(31)%kindof=8 ; gen(31)%diffu=0d0 ; gen(31)%mu=5d-1 ;gen(31)%name="s2 epi. knot supra receptor active"
      gen(31)%npre=2 ; allocate(gen(31)%pre(gen(31)%npre)) ; gen(31)%pre(1)=25 ; gen(31)%pre(2)=29
      gen(31)%npost=2 ; allocate(gen(31)%post(gen(31)%npost)) ; gen(31)%post(1)=25 ; gen(31)%post(2)=29
    


 
    !Gene-behavior interactions
    
   
      gen(8)%wa(1)=1
      gen(9)%wa(1)=2
      gen(10)%wa(1)=3
      gen(22)%wa(1)=4
      gen(23)%wa(1)=5
      gen(24)%wa(1)=6

      !gen(19)%wa(1)=5
      !gen(20)%wa(1)=6
      !gen(21)%wa(1)=7
      
      gen(22)%wa(10)=1d0
      
      gen(11)%wa(nparam_per_node+2)=-5d5

      gen(19)%wa(nparam_per_node+2)=1.2d-2
      !gen(20)%wa(nparam_per_node+2)=3.2d-1
      gen(20)%wa(nparam_per_node+2)=2.8d-2
      gen(21)%wa(nparam_per_node+2)=1.2d-2


      !gen(4)%wa(nparam_per_node+2)=0.5d-2
      !gen(6)%wa(nparam_per_node+2)=1.1d-2
      !gen(3)%wa(nparam_per_node+2)=0.5d-2
      
      
      !growth parameters for inward bud
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=3.2d-1
      !gen(21)%wa(nparam_per_node+2)=5d-2
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !growth parameters for outward bud
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=5.0d-3!5.0d-2
      !gen(21)%wa(nparam_per_node+2)=3.2d-1
      !gen(23)%wa(nparam_per_node+2)=-5.0d0
      !!!!!!!!!!!!!!!!!!!!!!!!!!   
      !growth parameters for rectangle scale
      !gen(19)%wa(nparam_per_node+2)=1.2d-1
      !gen(20)%wa(nparam_per_node+2)=0.0d0
      !gen(21)%wa(nparam_per_node+2)=1.0d-3
      !gen(23)%wa(nparam_per_node+2)=-7d1
      !gen(6)%wa(nparam_per_node+2)=1.0d-1
      !!!!!!!!!!!!!!!!!!!!!!!!!!   
      
    !Gene-gene interactions

      !structural genes
      !auto-maintenance
      gen(1)%w(1)=1d0
      gen(2)%w(2)=1d0
      gen(3)%w(3)=1d0
      gen(4)%w(4)=1d0
      gen(5)%w(5)=1d0
      gen(6)%w(6)=1d0
      gen(7)%w(7)=1d0
      gen(11)%w(11)=1d0
      !maintenance of stable gene expression
      gen(13)%w(4)=1d0
      gen(14)%w(6)=1d0
      gen(15)%w(3)=1d0
      gen(8)%w(1)=1d0
      gen(9)%w(2)=1d0
      gen(10)%w(3)=1d0
      gen(22)%w(11)=1d0
      gen(8)%w(11)=-1d5
      gen(9)%w(11)=-1d5
      !!!!!!!!!!!!!!!!
      gen(12)%w(11)=1d0
      gen(13)%w(11)=-1d5
      gen(14)%w(11)=-1d5
      gen(14)%w(1)=-1d5
      gen(13)%w(2)=-1d5

      !!!!!!!!!!!!!!!!! effector activation thresholds
      gen(19)%w(1)=-0.05d0
      gen(20)%w(2)=-0.05d0
      gen(21)%w(3)=-0.05d0
      !gen(19)%w(1)=0d0
      !gen(20)%w(2)=0d0
      !gen(21)%w(3)=0d0
      gen(19)%w(16)=5d1   !tooth
      gen(20)%w(17)=5d1
      gen(21)%w(18)=1d1
      !gen(19)%w(16)=1d2    !scale
      !gen(20)%w(17)=1d2
      !gen(21)%w(18)=1d1
      !gen(19)%w(16)=1d1   !tooth
      !gen(20)%w(17)=1d1
      !gen(21)%w(18)=1d1
      !gen(19)%w(16)=1d2    !rectangle scale
      !gen(20)%w(17)=1d2
      !gen(21)%w(18)=1d1


      !new interactions for s2    

      !gen(12)%w(26)=1d0 !s2 !the knot markers secrete the primary signal
      !gen(12)%w(27)=1d0 !s2
      !gen(12)%w(11)=0d0 !s2  !primary signal secretion is transfered to the knot markers
      !gen(25)%w(11)=1d0 !s2


      gen(23)%w(11)=1d0 !s2  !adhesion molecule activation
      gen(23)%w(26)=1d0 !s2
      gen(23)%w(19)=1d0 !s2
      gen(24)%w(21)=1d0 !s2

      gen(28)%w(4)=1d0 !s2 !activation of receptors
      gen(29)%w(4)=1d0 !s2
      gen(28)%w(2)=-1d5 !s2 !we restrict every receptor to a specific compartment
      gen(29)%w(1)=-1d5 !s2 !we restrict every receptor to a specific compartment


      gen(26)%w(1)=-0.05d0 !s2 !threshold activation for the knot markers
      gen(27)%w(2)=-0.05d0 !s2
      gen(26)%w(30)=2d1 !s2
      gen(27)%w(31)=2d1 !s2
      gen(19)%w(26)=-1d5 !s2
      gen(20)%w(27)=-1d5 !s2
 
 

      !signal-receptor interactions
      !epithelial basal binding
      gen(16)%nww=4
      gen(16)%ww(1,1)=12
      gen(16)%ww(1,2)=16
      gen(16)%ww(1,3)=1d0
      gen(16)%ww(2,1)=16
      gen(16)%ww(2,2)=12
      gen(16)%ww(2,3)=1d0
      gen(16)%ww(3,1)=13
      gen(16)%ww(3,2)=16
      gen(16)%ww(3,3)=1d0
      gen(16)%ww(4,1)=16
      gen(16)%ww(4,2)=13
      gen(16)%ww(4,3)=1d0
      !epithelial suprabasal receptor binding
      gen(17)%nww=4
      gen(17)%ww(1,1)=12
      gen(17)%ww(1,2)=17
      gen(17)%ww(1,3)=1d0
      gen(17)%ww(2,1)=17
      gen(17)%ww(2,2)=12
      gen(17)%ww(2,3)=1d0
      gen(17)%ww(3,1)=14
      gen(17)%ww(3,2)=15
      gen(17)%ww(3,3)=1d0
      gen(17)%ww(4,1)=17
      gen(17)%ww(4,2)=14
      gen(17)%ww(4,3)=1d0
      !mesenchymal receptor binding
      gen(18)%nww=4
      gen(18)%ww(1,1)=12
      gen(18)%ww(1,2)=18
      gen(18)%ww(1,3)=1d0
      gen(18)%ww(2,1)=18
      gen(18)%ww(2,2)=12
      gen(18)%ww(2,3)=1d0
      gen(18)%ww(3,1)=15
      gen(18)%ww(3,2)=18
      gen(18)%ww(3,3)=1d0
      gen(18)%ww(4,1)=18
      gen(18)%ww(4,2)=15
      gen(18)%ww(4,3)=1d0
     

      !knot receptors
      !epithelial basal binding
      gen(30)%nww=4
      gen(30)%ww(1,1)=25
      gen(30)%ww(1,2)=30
      gen(30)%ww(1,3)=1d0
      gen(30)%ww(2,1)=30
      gen(30)%ww(2,2)=25
      gen(30)%ww(2,3)=1d0
      gen(30)%ww(3,1)=28
      gen(30)%ww(3,2)=30
      gen(30)%ww(3,3)=1d0
      gen(30)%ww(4,1)=30
      gen(30)%ww(4,2)=28
      gen(30)%ww(4,3)=1d0
      !epithelial suprabasal receptor binding
      gen(31)%nww=4
      gen(31)%ww(1,1)=25
      gen(31)%ww(1,2)=31
      gen(31)%ww(1,3)=1d0
      gen(31)%ww(2,1)=31
      gen(31)%ww(2,2)=25
      gen(31)%ww(2,3)=1d0
      gen(31)%ww(3,1)=29
      gen(31)%ww(3,2)=31
      gen(31)%ww(3,3)=1d0
      gen(31)%ww(4,1)=31
      gen(31)%ww(4,2)=29
      gen(31)%ww(4,3)=1d0



    !Adhesion molecules

	ntipusadh=6
    if(ntipusadh>0)then
      if (allocated(kadh)) deallocate(kadh)
      allocate(kadh(ntipusadh,ntipusadh))
      kadh=0d0
      !Adhesion molecules interactions

      !as default, 1 and 2 should act as the same (basal and suprabasal)
      kadh(1,1)=1d0

      kadh(1,2)=1d0 ; kadh(2,1)=kadh(1,2)
      kadh(2,2)=1d0

      kadh(1,3)=5d-1 ; kadh(3,1)=kadh(1,3)
      !rectangle scale
      !kadh(1,3)=3d0 ; kadh(3,1)=kadh(1,3)
      kadh(2,3)=0d0 ; kadh(3,2)=kadh(2,3)
      kadh(3,3)=1d0

      kadh(1,4)=5d-1 ; kadh(4,1)=kadh(1,4)
      kadh(2,4)=5d-1 ; kadh(4,2)=kadh(2,4)
      kadh(3,4)=5d-1 ; kadh(4,3)=kadh(3,4)
      kadh(4,4)=1d1 !tooth
      !kadh(4,4)=5d0  !scale

      kadh(1,5)=0d0 ; kadh(5,1)=kadh(1,5)
      kadh(2,5)=0d0 ; kadh(5,2)=kadh(2,5)
      kadh(3,5)=0d0 ; kadh(5,3)=kadh(3,5)
      kadh(4,5)=0d0 ; kadh(5,4)=kadh(4,5)
      kadh(5,5)=0d0

      kadh(1,6)=0d0 ; kadh(6,1)=kadh(1,6)
      kadh(2,6)=0d0 ; kadh(6,2)=kadh(2,6)
      kadh(3,6)=0d0 ; kadh(6,3)=kadh(3,6)
      kadh(4,6)=0d0 ; kadh(6,4)=kadh(4,6)
      kadh(5,6)=1d1 ; kadh(6,5)=kadh(5,6)
      kadh(6,6)=1d1

      !kadh(1,7)=5d-1 ; kadh(7,1)=kadh(1,7)
      !kadh(2,7)=0d0 ; kadh(7,2)=kadh(2,7)
      !kadh(3,7)=1d0 ; kadh(7,3)=kadh(3,7)
      !kadh(4,7)=5d-1 ; kadh(7,4)=kadh(4,7)
      !kadh(5,7)=5d-1 ; kadh(7,5)=kadh(5,7)
      !kadh(6,7)=0d0 ; kadh(7,6)=kadh(6,7)
      !kadh(7,7)=1d0

    end if

    !Gene expression on nodes

    
    !if(geometry==1)then
    !  gex(1,2)=1d0
    !  gex(ndepi+1,2)=1d0
    !  !gex(ndepi+1,2)=1d0
    !end if
    
    if(geometry==2)then
    !draw rectangles
    !!!!!!!!!!!!!!!!!!!!!!!!!
    jj1=1 ; jj2=ly
    ii1=1 ; ii2=lx
    !placode range
    px1=3 ; px2=lx-(px1-1)
    py1=3 ; py2=ly-(py1-1)
    !signalling centre
    sx1=4 ; sx2=lx-(sx1-1)
    sy1=5 ; sy2=ly-(sy1-1)
    !AP proliferation supressor
    apx1=7 ; apx2=lx
    apy1=1 ; apy2=ly

    do i=ii1,ii2
      do j=jj1,jj2
        !epithelium
        k=cels(cell_grid_epi(i,j))%node(1)

        gex(k,1)=1d0 ; gex(node(k)%altre,1)=1d0
        gex(k,8)=1d0 ; gex(node(k)%altre,8)=1d0
        !borders
        !if(i<=1.or.i>=lx.or.j<=1.or.j>=ly)then
        !  node(k)%hold=3
        !  node(node(k)%altre)%hold=3
        !end if
        !placode
        if((i>px1.and.i<px2))then
          gex(k,4)=1d0 ; gex(node(k)%altre,4)=1d0
          gex(k,13)=1d0 ; gex(node(k)%altre,13)=1d0
        else
          gex(k,5)=1d0 ; gex(node(k)%altre,5)=1d0
        end if
        !signalling centre
        if(i>sx1.and.i<sx2)then
          gex(k,11)=1d0 ; gex(node(k)%altre,11)=1d0
          gex(k,22)=1d0 ; gex(node(k)%altre,22)=1d0
          gex(k,13)=0d0 ; gex(node(k)%altre,13)=0d0
        end if
        !!AP proliferation supressor
        !if((i>apx1.and.i<apx2).and.(j>apy1.and.j<apy2))then
        !  gex(k,23)=1d0
        !end if


        !mesench
        do ii=1,layer
          k=cels(cell_grid_mes(i,j,ii))%node(1)

          gex(k,3)=1d0
          gex(k,10)=1d0
          !borders
          !if(i<=1.or.i>=lx.or.j<=1.or.j>=ly)then
          !  node(k)%hold=3
          !end if
          !placode
          !if((i>px1.and.i<px2))then
            !gex(k,6)=1d0
            gex(k,15)=1d0
          !else
            !gex(k,7)=1d0
          !end if
        end do


        !supra
        do ii=1,xlayer
          k=cels(cell_grid_supra(i,j,ii))%node(1)

          gex(k,2)=1d0
          gex(k,9)=1d0
          !borders
          !if(i<=1.or.i>=lx.or.j<=1.or.j>=ly)then
          !  node(k)%hold=3
          !end if
          !placode
          if((i>px1.and.i<px2).and.ii==1)then
            gex(k,6)=1d0
            gex(k,14)=1d0
          else
            gex(k,7)=1d0
          end if
        end do



      end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!
    end if

  if(geometry==1)then  
    !ring=2  !scale
    ring=1 !tooth

    ring2=3

    polarized=0 !tooth
    !polarized=1 !scale

    do i=1,nd
      if(node(i)%tipus<3)then
        !if(node(i)%tipus==2)then
          !epithelial basal layer, basal side
          gex(i,1)=1d0  !receptors
          node(i)%marge=0
          d=sqrt(node(i)%x**2+node(i)%y**2)
          !setting expression on placode and interplacode
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,4)=1d0
            if(polarized==0.or.node(i)%x<=0.15)then
              gex(i,13)=1d0 !; gex(node(i)%altre,13)=1d0
            else
              !gex(i,23)=1d0
              gex(i,13)=1d0 !; gex(node(i)%altre,13)=1d0
            end if
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,5)=1d0
          end if
          !!!!!!!!!!!!!
          !signalling centre
          if(d-2*node(i)%req*(ring-1)<epsilod)then
            node(i)%marge=0 !; node(node(i)%altre)%marge=0
            gex(i,11)=1d0 !; gex(node(i)%altre,11)=1d0
            gex(i,22)=1d0 !; gex(node(i)%altre,22)=1d0
            gex(i,13)=0d0 !; gex(node(i)%altre,13)=0d0
          end if

        !end if
        !epithelial basal layer, whole layer       
        gex(i,8)=1d0

      end if
      if(node(i)%tipus==3)then
        if(node(i)%z<=0)then
          !mesenchymal layer
          gex(i,3)=1d0
          gex(i,10)=1d0

          !setting expression on placode and interplacode
          d=sqrt(node(i)%x**2+node(i)%y**2)
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,6)=1d0
            gex(i,15)=1d0
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,7)=1d0
          end if
          !!!!!!!!!!!!!

        else
          !epithelial suprabasal layer
          gex(i,2)=1d0
          gex(i,9)=1d0

          !setting expression on placode and interplacode
          d=sqrt(node(i)%x**2+node(i)%y**2)
          if(d-2*node(i)%req*(ring2-1)<epsilod)then
            gex(i,4)=1d0
            if(polarized==0.or.node(i)%x<=0.15)then
              gex(i,14)=1d0
            else
              gex(i,23)=1d0
              gex(i,14)=1d0
            end if
          end if
          if(d-2*node(i)%req*(ring2-1)>=epsilod)then
            gex(i,5)=1d0
          end if
          !!!!!!!!!!!!!
          !signalling centre
          !if(d-2*node(i)%req*(ring-1)<epsilod)then
          !  gex(i,11)=1d0
          !  gex(i,22)=1d0
          !  gex(i,14)=0d0 !; gex(node(i)%altre,13)=0d0
          !end if


        end if
      end if


    end do
    end if

    call update_npag

    node(:)%talone=0.0d0

    !ramax=maxval(node(:)%da)*3

    
end subroutine

















end module
