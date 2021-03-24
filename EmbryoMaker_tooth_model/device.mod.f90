!  afer:

!-opcio de fastmath i O2 al pgfortran
!-la part genetica i lo de la alocatacio
module device
use general
use genetic
use neighboring
use cudafor

integer,public,parameter :: max_ww=5, max_prepost=5  !max sizes for the gene type intrinsic matrices
integer, public, parameter     ::ngmax=20       ! max number of genes, for static allocatation of matrices !>>Miquel27-11-14

integer, public :: nblocks,nthreads_per_block
!integer, device, public :: d_nthreads_per_block

!!!!!!THE GEN VARIABLE TYPE CANNOT HAVE ALLOCATABLE ARRAYS WHEN DECLARED IN THE GPU,!!!!!!
!!!!!!THUS WE MAKE A SLIGHTLY DIFFERENT TYPE OF GEN TO BE DECLARED IN THE GPU!!!!!!!!!


  !real*8,allocatable,device   ::gen_w(:,:)      ! interaction strenghts of the genes regulating the transcription of gene i
  !integer,allocatable,device  ::gen_nww(:)       ! number of reactions catalyzed by i
  !real*8,allocatable,device   ::ww(:,:)   ! which reactions i regulates, dim1 index of the reaction(as in nww), dim2(1) the pre form, dim2(2) the post reaction, dim2(3) the w strength
!  real*8               ::diffu     ! extracellular diffusivity of gene proteins or protein forms
!  real*8               ::idiffu    ! intracellular diffusivity of gene proteins or protein forms
!  real*8               ::mu        ! degradation of gene proteins or protein forms
!  integer              ::npre      ! the number of pre form that gives rise to this form 
!  integer,allocatable,managed  ::pre(:)    ! the pre forms of a gene: the forms that give rise to it, NOTICE the reverse reaction only occurs if the post of a gene is its pre
!                                   ! if that does not occur it means that the reaction is effectivelly irreversible
!  integer              ::npost     ! number of post forms of a gene
!  integer,allocatable,managed  ::post(:)   ! the forms that arises from this form
!  real*8               ::kindof    ! 1 transcribed/translated protein or RNA without further modification; 
!  real*8, allocatable,managed  ::wa(:)    ! effect of a gene on each n








type, public :: dev_genes

  !of genes on genes
  !character*50         ::name  !!!ECTO MOD
  real*8,allocatable,managed   ::w(:)      ! interaction strenghts of the genes regulating the transcription of gene i
  integer              ::nww       ! number of reactions catalyzed by i
  real*8,allocatable,managed   ::ww(:,:)   ! which reactions i regulates, dim1 index of the reaction(as in nww), dim2(1) the pre form, dim2(2) the post reaction, dim2(3) the w strength
  real*8               ::diffu     ! extracellular diffusivity of gene proteins or protein forms
  real*8               ::idiffu    ! intracellular diffusivity of gene proteins or protein forms
  real*8               ::mu        ! degradation of gene proteins or protein forms
  integer              ::npre      ! the number of pre form that gives rise to this form 
  integer,allocatable,managed  ::pre(:)    ! the pre forms of a gene: the forms that give rise to it, NOTICE the reverse reaction only occurs if the post of a gene is its pre
                                   ! if that does not occur it means that the reaction is effectivelly irreversible
  integer              ::npost     ! number of post forms of a gene
  integer,allocatable,managed  ::post(:)   ! the forms that arises from this form
  real*8               ::kindof    ! 1 transcribed/translated protein or RNA without further modification; 
                                   ! 2 transcribed/translated protein or RNA without that can be modified postraductionally
                                   ! 3 comes from a previous form   
                                   ! 4 extracelular signal: NOTICE THE MUST NOT HAVE A PRE: since they get secreted immediately and thus is like an irreversible reaction
                                   ! 5 form that is intrinsically transported to the apical side of an epithelial cell (through kinesins)  !in these two cases, transport depends only
                                   ! 6 form that is intrinsically transported to the basal side of an epithelial cell (through dineins)    !on the absolute concentration on one node,
                                   ! 7 membrane-bound signal/receptor: like notch-delta ; it requires a pre since 7 is only the activated form of the receptor
                                   ! 8 bound receptor of extracellular signal, it has to have one (and only one) pre, and promote its own change to the 
                                                                         !pre form (unbounding of the ligand)
                                   ! notice:                                                                                               !so the %diffus probably should be smaller.
				   !          extracellular signals are automatically secreted and can only affect genes in a node through receptors
				   !

!  character*100        ::label      !this is intended to be used as a label to identify and track the gene (trhough evolution for example) !>>Miquel8-8-14


  !of genes on cell behaviours or properties
  real*8, allocatable,managed  ::wa(:)    ! effect of a gene on each node parameter. These are ordered by number; after the number of the node parameter
                                  ! come each of the cell behaviours
                                  ! nparam_per_node+1 = cell growth
                                  ! nparam_per_node+2 = cell cycle increase, when 1 the cell can divide if it has the right size
                                  ! nparam_per_node+3 = apoptosis
                                  ! nparam_per_node+4 = this gene promotes secretion at this rate 
                                  ! nparam_per_node+5 = 1 means this gene is secretable 
                                  !                     there is no parameter for adh because that comes from the genes expressed there
                                  ! nparam_per_node+6 = repcel of the secreted node
                                  ! nparam_per_node+7 = da (as prportion of reqcel) of the secreted node
                                  ! nparam_per_node+8 = this gene is polarizing the cell [affecting pol vector of cel]
                                  ! nparam_per_node+9 = this gene is telling the cells that they should grow in the polarized direction with
                                  !                     a noise that is 1-(this value)
                                  ! nparam_per_node+10= this affects the cell property minsize_for_div, the number of nodes required for a cell
                                  !                     to divide
                                  ! nparam_per_node+11= dependence of the plane of division relative to chemical gradients over the Hertwig vector
                                  ! nparam_per_node+12= activates assymetric division
                                  ! nparam_per_node+13= activates epitelial-to-mesenchymal transition (EMT)
                                  ! nparam_per_node+14= ECM proteolisis (like apoptosis for ECM nodes)
                                  ! nparam_per_node+15= changes maxsize_for_div the max number of nodes that allows division
                                  ! nparam_per_node+16= noise is biased by polarization vector
end type












!type, public :: dev_genes
!
!  !of genes on genes
!  !character*50         ::name  !!!ECTO MOD
!  real*8,dimension(:)   ::w(ngmax)      ! interaction strenghts of the genes regulating the transcription of gene i
!  integer              ::nww       ! number of reactions catalyzed by i
!  real*8,dimension(:,:)   ::ww(max_ww,3)   ! which reactions i regulates, dim1 index of the reaction(as in nww), dim2(1) the pre form, dim2(2) the post reaction, dim2(3) the w strength
!  real*8               ::diffu     ! extracellular diffusivity of gene proteins or protein forms
!  real*8               ::idiffu    ! intracellular diffusivity of gene proteins or protein forms
!  real*8               ::mu        ! degradation of gene proteins or protein forms
!  integer              ::npre      ! the number of pre form that gives rise to this form 
!  integer,dimension(:)  ::pre(max_prepost)    ! the pre forms of a gene: the forms that give rise to it, NOTICE the reverse reaction only occurs if the post of a gene is its pre
!                                   ! if that does not occur it means that the reaction is effectivelly irreversible
!  integer              ::npost     ! number of post forms of a gene
!  integer,dimension(:)  ::post(max_prepost)   ! the forms that arises from this form
!  real*8               ::kindof    ! 1 transcribed/translated protein or RNA without further modification; 
!                                   ! 2 transcribed/translated protein or RNA without that can be modified postraductionally
!                                   ! 3 comes from a previous form   
!                                   ! 4 extracelular signal: NOTICE THE MUST NOT HAVE A PRE: since they get secreted immediately and thus is like an irreversible reaction
!                                   ! 5 form that is intrinsically transported to the apical side of an epithelial cell (through kinesins)  !in these two cases, transport depends only
!                                   ! 6 form that is intrinsically transported to the basal side of an epithelial cell (through dineins)    !on the absolute concentration on one node,
!                                   ! 7 membrane-bound signal/receptor: like notch-delta ; it requires a pre since 7 is only the activated form of the receptor
!                                   ! 8 bound receptor of extracellular signal, it has to have one (and only one) pre, and promote its own change to the 
!                                                                         !pre form (unbounding of the ligand)
!                                   ! notice:                                                                                               !so the %diffus probably should be smaller.
!				   !          extracellular signals are automatically secreted and can only affect genes in a node through receptors
!				   !
!
!!  character*100        ::label      !this is intended to be used as a label to identify and track the gene (trhough evolution for example) !>>Miquel8-8-14
!
!
!  !of genes on cell behaviours or properties
!  real*8, dimension(:)  ::wa(nga)    ! effect of a gene on each node parameter. These are ordered by number; after the number of the node parameter
!                                  ! come each of the cell behaviours
!                                  ! nparam_per_node+1 = cell growth
!                                  ! nparam_per_node+2 = cell cycle increase, when 1 the cell can divide if it has the right size
!                                  ! nparam_per_node+3 = apoptosis
!                                  ! nparam_per_node+4 = this gene promotes secretion at this rate 
!                                  ! nparam_per_node+5 = 1 means this gene is secretable 
!                                  !                     there is no parameter for adh because that comes from the genes expressed there
!                                  ! nparam_per_node+6 = repcel of the secreted node
!                                  ! nparam_per_node+7 = da (as prportion of reqcel) of the secreted node
!                                  ! nparam_per_node+8 = this gene is polarizing the cell [affecting pol vector of cel]
!                                  ! nparam_per_node+9 = this gene is telling the cells that they should grow in the polarized direction with
!                                  !                     a noise that is 1-(this value)
!                                  ! nparam_per_node+10= this affects the cell property minsize_for_div, the number of nodes required for a cell
!                                  !                     to divide
!                                  ! nparam_per_node+11= dependence of the plane of division relative to chemical gradients over the Hertwig vector
!                                  ! nparam_per_node+12= activates assymetric division
!                                  ! nparam_per_node+13= activates epitelial-to-mesenchymal transition (EMT)
!                                  ! nparam_per_node+14= ECM proteolisis (like apoptosis for ECM nodes)
!                                  ! nparam_per_node+15= changes maxsize_for_div the max number of nodes that allows division
!                                  ! nparam_per_node+16= noise is biased by polarization vector
!end type

!!!!!NODE VARIABLES!!!!!!!!
type(nod),device,allocatable,dimension(:) :: d_node
integer, device ::d_nd,d_mnn

!real*8,device,allocatable,dimension(:,:) :: d_pxpypz !d_dex,d_px,d_py,d_pz,d_fmeanl


!!!!!!GENETIC VARIABLES!!!!!!
type(dev_genes),managed,allocatable,dimension(:) ::d_gen
real*8,device,allocatable,dimension(:,:) :: d_gex
real*8,device,allocatable, dimension(:,:):: d_kadh


integer, device, allocatable,dimension(:) :: d_npag   !number of genes that have an effect in that cellular parameter
integer, device, allocatable,dimension(:,:) :: d_whonpag !list of the indexes of the genes affecting each cellular parameter
!real*8,device,allocatable, dimension(:,:):: d_gext
!real*8,device,allocatable, dimension(:,:):: d_dgex
!integer, device, allocatable,dimension(:) :: d_checkn
integer, device :: d_ng,d_nga


! derived matrices that do not change and are just for making the code run faster

integer,device              ::d_nkindof(9)   ! number of each kindof gene
integer,device,allocatable  ::d_wkindof(:,:) ! which are those for each kindof

real*8 , device, allocatable :: d_nw(:),d_w(:,:)
integer, device, allocatable :: d_nwpre(:,:)    ! number of forms catalyzing the reaction pre j of i 
real*8 , device, allocatable :: d_wpre(:,:,:,:) ! dim1 the form i, dim2 index of the reactions leading to i from its pres (as in npre)
                                               ! dim3(1) the list of enzymes, dim4(1) their gene indices dim4(2) their ws  
integer, device, allocatable :: d_nwpost(:,:)   ! number of forms catalyzing the reaction post j of i 
real*8 , device, allocatable :: d_wpost(:,:,:,:)! dim1 the form i, dim2 index of the reactions leading to i from its posts (as in npost)
                                               ! dim3(1) the list of enzymes, dim4(1) their gene indices dim4(2) their ws


!!!!!!NEIGHBORING VARIABLES!!!!!!

integer, device, allocatable,dimension(:) :: d_nneigh  !number of neighbors
integer, device, allocatable,dimension(:,:) :: d_neigh  !neighbor matrix
real*8, device, allocatable,dimension(:,:) :: d_dneigh  !neighbor matrix distances
integer, device,allocatable,dimension(:,:) :: d_trans_neigh(:,:)  !neighbor matrix
real*8, device,allocatable,dimension(:,:) :: d_trans_dneigh(:,:)  !neighbor matrix distances

integer, device, allocatable,dimension(:,:) :: d_borders_neigh  !neighbor matrix

integer, device, allocatable,dimension(:,:,:):: d_boxes !number of elements per grid box
integer, device ,allocatable,dimension(:) :: d_list !neighbor building grid matrix
integer, device ::d_nboxes,d_fscreen
real*8, device ::d_screen_radius
real*8, device ::d_urv,d_rdiffmax
real*8, device :: d_maxlen
real*8 :: maxlen
integer,device::d_nparam_per_node

!!!!!GENERAL PARAMETERS!!!!!!
real*8, device ::d_dmax,d_resmax,d_delta,d_epsilod,d_angletor
integer, device ,allocatable,dimension(:) :: d_ffu
!integer, device:: d_n(5000),d_nn(5000),d_nnn(5000),d_nnnn(5000)
!integer :: n(5000),nn(5000),nnn(5000),nnnn(5000)

!integer::seldev

contains

!******************************************************************************************************************************

subroutine init_dev(seldev)
integer:: i,istat,ndev
integer,value::seldev
type(cudadeviceprop):: prop



!***********SET THE DEVICE*************!

  istat=cudaGetDeviceCount(ndev)
  print*,"******number of devices=",ndev
  print*,"******listing the devices"
  do i=0,ndev-1
    istat=cudaGetDeviceProperties( prop, i )
    !print*,i,"properties",prop

    print *,"index",i,"name",prop%name
    print *,"shared memory per block in bytes",prop%sharedMemPerBlock
    print *,"registers per block in bytes",prop%regsPerBlock
    !print *,prop%major,prop%minor,"compute capability"
    !print *,prop%multiprocessorcount,"number of multiprocessors"
    !print *,prop%clockrate,"clock rate"
    !print *,prop%computemode,"compute mode"
    !print *,prop%concurrentkernels,"concurrent kernels"
    !print *,prop%pcideviceid,"device pci id"
    !print *,prop%kernelexectimeoutenabled,"kernelexectimeoutenabled"
    !print *,prop%maxthreadspermultiprocessor,"maxthreadspermultiprocessor"

  end do

  !seldev=0  ! it should be from 0 to 2
  print*,"******setting device number",seldev
  istat=cudaSetDevice(seldev)
  istat=cudaGetDeviceProperties( prop, seldev )
  print *,"device selected",prop%name

  !print *,istat,"istat ty"

!initialize all device variables, now used for forces an gene_stuff



!!!!!NODE INITIALIZATIONS!!!!!!!


  if (allocated(d_node)) deallocate(d_node)
  allocate(d_node(nda))


  d_node(1:nda)=node(1:nda)


  d_nd=nd 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!NEIGHBORING INITIALIZATIONS!!!!!!!
  mnn_dynam=mnn
  omnn=0
  !print*,"mnn_dynam",mnn_dynam
  !print*,"debug:: before allocating d_neigh",nd,mnn_dynam

  if(allocated(d_neigh)) deallocate(d_neigh) !only for use within the device
  allocate(d_neigh(nda,mnn_dynam))
  d_neigh=0
  !print*,"debug:: after allocating d_neigh"
  if(allocated(d_nneigh)) deallocate(d_nneigh)
  allocate(d_nneigh(nda))
  d_nneigh=0
  !print*,"debug:: after allocating d_nneigh"

  if(allocated(d_dneigh)) deallocate(d_dneigh) !only for use within the device
  allocate(d_dneigh(nda,mnn_dynam))
  d_dneigh=0
  !print*,"debug:: after allocating d_dneigh"

  !if(allocated(d_trans_neigh)) deallocate(d_trans_neigh) !only for use within the device
  !allocate(d_neigh(nda,mnn_dynam))
  !d_trans_neigh=0
  !
  !if(allocated(d_trans_dneigh)) deallocate(d_trans_dneigh) !only for use within the device
  !allocate(d_trans_dneigh(nda,mnn_dynam))
  !d_trans_dneigh=0d0
  
  d_nboxes=nboxes
  if(ffu(3)==1)then
    d_fscreen=1
    d_screen_radius=screen_radius
  else
    d_fscreen=0
  end if

  !if(prop_noise>epsilod)then
    if(allocated(neigh)) deallocate(neigh) !only for use within the device
    allocate(neigh(nda,mnn_dynam))
    neigh=0
    if(allocated(nneigh)) deallocate(nneigh)
    allocate(nneigh(nda))
    nneigh=0
    if(allocated(dneigh)) deallocate(dneigh) !only for use within the device
    allocate(dneigh(nda,mnn_dynam))
    dneigh=0
  !end if  
  
  if(allocated(d_boxes)) deallocate(d_boxes)
  allocate(d_boxes(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))
  d_boxes=boxes
  !print*,"debug:: after initializing d_boxes"

  if(allocated(d_list)) deallocate(d_list)
  allocate(d_list(nda))
  !d_list(1:nd)=list(1:nd)
  !istat=cudaMalloc(d_list,nd)
  
  
  !if(ffu(24)==1)then
  !  allocate(d_borders_neigh(nborders,4))
  !  d_borders_neigh=borders_neigh
  !end if
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!GENETIC INITIALIZATIONS!!!!!!!!!!

  if (allocated(d_gen)) deallocate(d_gen)
  allocate(d_gen(ng))


  do i=1,ng   !ARRAYS WITHIN TYPE ARRAYS MUST BE FIXED SIZE  

    !w
    !if (allocated(d_gen(i)%w)) deallocate(d_gen(i)%w)
    allocate(d_gen(i)%w(ng))
    istat = cudaDeviceSynchronize()
    d_gen(i)%w=gen(i)%w

    !d_gen(i)%w(1:ng)=gen(i)%w(1:ng)

    allocate(d_gen(i)%ww(ng*ng,3))
    istat = cudaDeviceSynchronize()

    d_gen(i)%nww=gen(i)%nww
    d_gen(i)%ww=gen(i)%ww
    !d_gen(i)%ww(1:gen(i)%nww,1:3)=gen(i)%ww(1:gen(i)%nww,1:3)

    !wa
    !if (allocated(d_gen(i)%wa)) deallocate(d_gen(i)%wa)
    allocate(d_gen(i)%wa(nga))
    istat = cudaDeviceSynchronize()
    d_gen(i)%wa=gen(i)%wa
    !d_gen(i)%wa(1:nga)=gen(i)%wa(1:nga)


    d_gen(i)%npre=gen(i)%npre
    d_gen(i)%npost=gen(i)%npost
    d_gen(i)%diffu=gen(i)%diffu
    d_gen(i)%kindof=gen(i)%kindof
    d_gen(i)%mu=gen(i)%mu


    if(gen(i)%npost>0)then
      !if (allocated(d_gen(i)%post)) deallocate(d_gen(i)%post)
      allocate(d_gen(i)%post(gen(i)%npost))
      istat = cudaDeviceSynchronize()
      d_gen(i)%post=gen(i)%post

      !d_gen(i)%post(1:gen(i)%npost)=gen(i)%post(1:gen(i)%npost)
       !print*,i,"d_post1",d_gen(i)%post(1)

    end if

    if(gen(i)%npre>0)then
      !if (allocated(d_gen(i)%pre)) deallocate(d_gen(i)%pre)
      allocate(d_gen(i)%pre(gen(i)%npre))
      istat = cudaDeviceSynchronize()
      d_gen(i)%pre=gen(i)%pre

      !d_gen(i)%pre(1:gen(i)%npre)=gen(i)%pre(1:gen(i)%npre)
    end if


    
  end do


 !do i=1,ng
 !  print*,i,"npost",gen(i)%npost
 !end do
 !
 !
 ! d_gen(1:ng)%npre=gen(1:ng)%npre
 ! d_gen(1:ng)%npost=gen(1:ng)%npost
 ! d_gen(1:ng)%diffu=gen(1:ng)%diffu
 ! d_gen(1:ng)%kindof=gen(1:ng)%kindof
 ! d_gen(1:ng)%mu=gen(1:ng)%mu
 !
 !do i=1,ng
 !  print*,i,"npost",d_gen(i)%npost
 !end do

  !do i=1,ng
  !  if(gen(i)%npost>0)then
  !    !if (allocated(d_gen(i)%post)) deallocate(d_gen(i)%post)
  !    !allocate(d_gen(i)%post(gen(i)%npost))
  !
  !     print*,i,"post1",gen(i)%post(1)
  !
  !    d_gen(i)%post(1:gen(i)%npost)=gen(i)%post(1:gen(i)%npost)
  !     print*,i,"d_post1",d_gen(i)%post(1)
  !
  !  end if
  !
  !  if(gen(i)%npre>0)then
  !    !if (allocated(d_gen(i)%pre)) deallocate(d_gen(i)%pre)
  !    !allocate(d_gen(i)%pre(gen(i)%npre))
  !    d_gen(i)%pre(1:gen(i)%npre)=gen(i)%pre(1:gen(i)%npre)
  !  end if
  !end do


  if (allocated(d_npag)) deallocate(d_npag)
  allocate(d_npag(nga))
  d_npag(1:nga)=npag(1:nga)

  if (allocated(d_whonpag)) deallocate(d_whonpag)  
  allocate(d_whonpag(nga,ng))
  d_whonpag(1:nga,1:ng)=whonpag(1:nga,1:ng)

  if (allocated(d_gex)) deallocate(d_gex)  
  allocate(d_gex(nda,ng))
  d_gex(1:nd,1:ng)=gex(1:nd,1:ng)


  !if (allocated(d_dgex)) deallocate(d_dgex)  
  !allocate(d_dgex(nd,ng))
  !d_dgex=0
  
  !if (allocated(d_checkn)) deallocate(d_checkn)  
  !allocate(d_checkn(nd))
  !d_checkn=0  


  d_ng=ng ; d_nga=nga


  if(ntipusadh>0)then
    if (allocated(d_kadh)) deallocate(d_kadh)  
    allocate(d_kadh(ntipusadh,ntipusadh))
    d_kadh=kadh
  end if

  
  ! derived utilitary matrices that do not change and are just for making the code run faster
  if (allocated(d_nw)) deallocate(d_nw)  
  allocate(d_nw(ng))
  d_nw=nw
  if (allocated(d_w)) deallocate(d_w)  
  allocate(d_w(ng,ng))
  d_w=w
  !print*,"ng*****************",ng
  !do i=1,ng
  !print*,"w",w(i,:)
  !end do
  
  k=maxval(gen(1:ng)%npre)
  if(k>0)then
    if (allocated(d_nwpre)) deallocate(d_nwpre)  
    allocate(d_nwpre(ng,k))
    d_nwpre=nwpre
    if (allocated(d_wpre)) deallocate(d_wpre)  
    allocate(d_wpre(ng,k,ng,2))
    d_wpre=wpre
  end if
  k=maxval(gen(1:ng)%npost)
  if(k>0)then
    k=maxval(gen(1:ng)%npost)
    if (allocated(d_nwpost)) deallocate(d_nwpost)  
    allocate(d_nwpost(ng,k))
    d_nwpost=nwpost
    if (allocated(d_wpost)) deallocate(d_wpost)  
    allocate(d_wpost(ng,k,ng,2))
    d_wpost=wpost
  end if
  
  
  if (allocated(d_wkindof)) deallocate(d_wkindof)  
  allocate(d_wkindof(9,ng))
  d_wkindof=wkindof
  d_nkindof=nkindof

  
  !!!!FORCES INITIALIZATIONS!!!!!
  !if (allocated(d_pxpypz)) deallocate(d_pxpypz)  
  !allocate(d_pxpypz(nd,5))
  !d_pxpypz(1:nd,:)=0d0

  !if (allocated(d_px)) deallocate(d_px)  
  !allocate(d_px(nd))
  !d_px(1:nd)=0d0
  !istat=cudaMalloc(d_px,nd)

  !if (allocated(d_py)) deallocate(d_py)  
  !allocate(d_py(nd))
  !d_py(1:nd)=0d0
  !istat=cudaMalloc(d_py,nd)
  
  !if (allocated(d_pz)) deallocate(d_pz)  
  !allocate(d_pz(nd))
  !d_pz(1:nd)=0d0
  !istat=cudaMalloc(d_pz,nd)

  !if (allocated(d_dex)) deallocate(d_dex)  
  !allocate(d_dex(nd))
  !d_dex(1:nd)=0d0
  !istat=cudaMalloc(d_dex,nd)

  !if (allocated(d_fmeanl)) deallocate(d_fmeanl)  
  !allocate(d_fmeanl(nd))
  !d_fmeanl(1:nd)=0d0
  !istat=cudaMalloc(d_fmeanl,nd)
  
  !!!!GLOBAL PARAMETERS!!!!!!!
  d_dmax=dmax
  d_urv=urv
  d_epsilod=epsilod
  d_angletor=angletor
  
  allocate(d_ffu(nfu))
  d_ffu=ffu
  d_mnn=mnn

  !nblocks=nd/1024+1
  !nthreads_per_block=nd/nblocks+1
  !nthreads_per_block=1024
  !nblocks=nd/nthreads_per_block+1

  if(nda>256)then
    nthreads_per_block=256
    nblocks=nda/nthreads_per_block+1
  else  
    nthreads_per_block=nda
    nblocks=1
  end if
  
  d_nparam_per_node=nparam_per_node

  print*,"getot",getot,"nd",nd
  print*,"nthreads_per_block",nthreads_per_block,"nblocks",nblocks

  !d_nthreads_per_block=nthreads_per_block

!print *,nblocks,nthreads_per_block,d_nthreads_per_block
end subroutine init_dev
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_device_pre_iteration !do this every iteration before launching any kernel
integer::i,ivar
integer::old_mnn_dynam,new_nd
!updates matrix sizes for gex and makes a copy to the device

  !print*,"1.1"
  new_nd=d_nd
  !temp2=d_ng
  !print*,"1.2"
  !print*,"nd",nd,"d_nd",temp
  !print*,"new_nd",new_nd,"nd",nd
  if(new_nd /= nd)then
!  if(temp /= nda)then
  !d_nda=nda
  !ivar=45896

    !nblocks=nd/1024+1
    !nthreads_per_block=nd/nblocks
    !d_nthreads_per_block=nthreads_per_block
    if(nda>256)then
      nthreads_per_block=256
      nblocks=nda/nthreads_per_block+1
    else  
      nthreads_per_block=nda
      nblocks=1
    end if
    !print*,"getot",getot,"nd",nd
    !print*,"nthreads_per_block",nthreads_per_block,"nblocks",nblocks


    !print*,"hi ha canvi"
    deallocate(d_node) ; allocate(d_node(nda),STAT=ivar)
    !print*,"debug:: no es d_node",ivar
    !print*,"ng",ng
    deallocate(d_gex) ; allocate(d_gex(nda,ng),STAT=ivar)
    !print*,"debug:: no es d_gex",ivar

    !deallocate(d_dgex) ; allocate(d_dgex(nd,ng),STAT=ivar)
    !print*,"debug:: no es d_dgex",ivar

    !deallocate(d_checkn) ; allocate(d_checkn(nd),STAT=ivar)

    !print*,"debug:: no es d_checkn",ivar
    deallocate(d_list) ; allocate(d_list(nda))
    !istat=cudaFree(d_list)
    !istat=cudaMalloc(d_list, nd)

    !print*,"debug:: canvi en nd"

    deallocate(d_neigh,d_dneigh) ; allocate(d_neigh(nda,omnn),d_dneigh(nda,omnn))
    if(allocated(d_nneigh)) deallocate(d_nneigh) ; allocate(d_nneigh(nda))
    !if(prop_noise>epsilod)then
      deallocate(neigh,dneigh) ; allocate(neigh(nda,omnn),dneigh(nda,omnn))
      deallocate(nneigh) ; allocate(nneigh(nda))
    !end if

    !deallocate(d_trans_neigh,d_trans_dneigh) ; allocate(d_trans_neigh(nda,mnn_dynam),d_trans_dneigh(nda,mnn_dynam))
    
    !print*,"debug:: no es veins"

    ! deallocate(d_px) ; allocate(d_px(nd))
    !!print*,"debug:: no es d_px"

    ! deallocate(d_py) ; allocate(d_py(nd))
    !!print*,"debug:: no es d_py"

    ! deallocate(d_pz) ; allocate(d_pz(nd))
    !print*,"debug:: no es d_pz"

    ! deallocate(d_dex) ; allocate(d_dex(nd))
    !print*,"debug:: no es d_dex"
    
    ! deallocate(d_fmeanl) ; allocate(d_fmeanl(nd))
    !print*,"debug:: no es d_fmeanl"
    
    !deallocate(d_pxpypz) ; allocate(d_pxpypz(nd,5))
    !print*,"debug:: no es d_fmeanl"
    d_nd=nd
  end if

  

  
  !if(mnn_dynam/=old_mnn_dynam.or.new_nd/=nd)then
  !  !print*,"hi ha canvi a veins"
  !  if(allocated(d_neigh)) deallocate(d_neigh)
  !      !print*,"deallocated"
  !  allocate(d_neigh(1:nd,1:mnn_dynam))
  !      !print*,"allocated d_neigh"
  !      !print*,"deallocated"  
  !  if(allocated(d_dneigh)) deallocate(d_dneigh)
  !  allocate(d_dneigh(1:nd,1:mnn_dynam))
  !      !print*,"allocated d_dneigh"
  !  !if(prop_noise>epsilod)then
  !  !  if(allocated(neigh)) deallocate(neigh)
  !  !  allocate(neigh(1:nd,1:mnn_dynam))
  !  !  if(allocated(dneigh)) deallocate(dneigh)
  !  !  allocate(dneigh(1:nd,1:mnn_dynam))
  !  !end if
  !end if

  
  d_node(1:nda)=node(1:nda)
  !print*,"1.4"
  !d_nneigh(1:nd)=0d0 ; d_neigh(1:nd,1:mnn)=0d0 ; d_dneigh(1:nd,1:mnn)=0d0
  !nneigh=0 ; neigh=0 ; dneigh=0d0
  !d_nneigh=0 ; d_neigh=0 ; d_dneigh=0d0

  !print*,"1.5"
  d_gex(1:nda,1:ng)=gex(1:nda,1:ng)
  !do i=1,nd
  !  !print*,i,"gex",gex(i,1:13)
  !end do
  !print*,"1.5.5"

  !d_dgex(1:nd,1:ng)=0

  !print*,"1.6"
  !d_dgex(1:nd,1:ng)=0
  !print*,"1.7"

  !print*,"1.9"
end subroutine set_device_pre_iteration

!**************************************************************************************************************
!!!!!NEIGHBOR BUILD

subroutine dev_neighboring(nblocks,nthreads_per_block,nd,nda)
integer,value::nblocks,nthreads_per_block,nd,nda
integer :: trans_nneigh(nda)  !number of neighbors
!integer, device,allocatable:: d_trans_neigh(:,:)  !neighbor matrix
!real*8, device, allocatable:: d_trans_dneigh(:,:)  !neighbor matrix distances
integer, device,dimension(:,:):: d_trans_neigh(nda,500)  !neighbor matrix
real*8, device,dimension(:,:):: d_trans_dneigh(nda,500)  !neighbor matrix distances

!integer,allocatable :: trans_neigh(:,:)  !neighbor matrix
!real*8,allocatable:: trans_dneigh(:,:)  !neighbor matrix distances
integer::old_mnn_dynam,new_omnn,grid_length,old_nboxes,shared_allocation
type(dim3)::block_grid

  trans_nneigh=0
  !d_trans_neigh=0
  !d_trans_dneigh=0d0

  if(prop_noise>epsilod .or. ffu(1)==1)then
    if(size(neigh,dim=2)/=size(d_neigh,dim=2))then
      !print*,"ADJUST******"
      deallocate(d_neigh,d_dneigh)
      allocate(d_neigh(nda,omnn),d_dneigh(nda,omnn))
    end if
  end if
  
  !REMAKE GRID BOXES!!!!!
  rv=2*maxval(node(:nd)%da);urv=1d0/rv !;print*,"RV",rv
  rdiffmax=2*maxval(node(:nd)%da)*dmax
  d_rdiffmax=rdiffmax

  !print*,"1.2.9"  
  old_nboxes=nboxes
  call iniboxes
  !d_maxlen=maxlen
  d_rv=rv ; d_urv=urv

  if(old_nboxes /= nboxes)then
    deallocate(d_boxes) ; allocate(d_boxes(-nboxes:nboxes,-nboxes:nboxes,-nboxes:nboxes))
  end if  
  d_list=list
  !print*,"1.7.5"
  d_boxes=boxes
  !print*,"1.8"
  !d_nd=nd !; d_nda=nda
  !d_ncels=ncels ; d_ncals=ncals
  d_nboxes=nboxes
  
  
  !print*,"1.3"

  
  !rdiffmax=2*maxval(node(:nd)%da)*dmax  
  old_mnn_dynam=mnn_dynam
  grid_length=2*nint(rdiffmax*urv)+1
  !print*,"debug :: rdiffmax",rdiffmax,"urv",urv

  mnn_dynam=mnn_dyn*grid_length**3 !setting dynamic neighbor matrix width !>>Miquel22-1-15

  !print*,"debug :: grid length",grid_length,"mnn_dyn",mnn_dyn,"mnn_dynam",mnn_dynam

  !print*,"mnn_dynam nou",mnn_dynam,"old_mnn_dynam",old_mnn_dynam
  !if(old_mnn_dynam/=mnn_dynam)then
  !    if (allocated(d_trans_neigh)) deallocate(d_trans_neigh)
  !    if (allocated(d_trans_dneigh)) deallocate(d_trans_dneigh)
  !    allocate(d_trans_neigh(nda,mnn_dynam),d_trans_dneigh(nda,mnn_dynam))
  !    d_trans_neigh=0 ; d_trans_dneigh=0d0
  !end if

  !allocate(trans_neigh(nda,mnn_dynam),trans_dneigh(nda,mnn_dynam))

  !allocate(d_trans_neigh(nda,mnn_dynam),d_trans_dneigh(nda,mnn_dynam))
  d_trans_neigh=0;d_trans_dneigh=0d0

  
  !print*,"shared allocation",shared_allocation,sizeof(a),sizeof(i)
  !print*,"neighboring 1.2"
  !print*,"pre_kernel new_omnn",new_omnn,"old omnn",omnn,"mnn_dynam",mnn_dynam,"nd",nd,"nda",nda,"nboxes",nboxes  

  !block_grid=dim3(grid_length,grid_length,grid_length)
  !shared_allocation=(4+8)*grid_length**3*(mnn_dyn+1)+4*grid_length**3+(4+8)*mnn_dynam  !memory space for dynamic shared arrays in the kernel
  !call neighbor_gather_kernel<<<nd,block_grid,shared_allocation>>>(d_trans_neigh,d_trans_dneigh,mnn_dynam,mnn_dyn,nd,maxlen) !>>>>>>GPU ADD 28-11-14

  !istat=cudaDeviceSynchronize

  !print*,"gather DONE",nd

!  allocate(trans_neigh(nda,mnn_dynam))
!
!  trans_nneigh(1:nda)=d_nneigh(1:nda)
!  !trans_neigh=d_trans_neigh
!
!  !trans_dneigh=d_trans_dneigh
! do i=1,nd
!  print*,"nneigh1",trans_nneigh(i)
!  print*,"neigh1",trans_neigh(i,1:trans_nneigh(i))
!end do
!  deallocate(trans_neigh)

  !print*,""
  !print*,"post_kernel ","old omnn",omnn,"mnn_dynam",mnn_dynam,"nd",nd,"nda",nda,"nboxes",nboxes

  
  !if(ffu(3)==1)then
  !  call neighbor_gabriel_kernel<<<nblocks,nthreads_per_block>>>(d_trans_neigh,d_trans_dneigh,mnn_dynam,screen_radius,nd) !>>>>>>GPU ADD 28-11-14
  !  istat=cudaDeviceSynchronize
  !  !print*,"passa?"
  !end if
  !print*,"gabriel DONE",nd

  !print*,"neighboring 1.1"

  call neighbor_build_kernel<<<nblocks,nthreads_per_block>>>(d_trans_neigh,d_trans_dneigh,d_nneigh,mnn_dynam,mnn_dyn,nd) !>>>>>>GPU ADD 28-11-14

  
    !allocate(trans_neigh(nda,mnn_dynam))

  trans_nneigh(1:nda)=d_nneigh(1:nda)
  !trans_neigh=d_trans_neigh

  !trans_dneigh=d_trans_dneigh
 !do i=1,nd
 ! print*,"nneigh1",trans_nneigh(i)
 ! print*,"neigh1",trans_neigh(i,1:trans_nneigh(i))
!end do
!  deallocate(trans_neigh)

  
  
  trans_nneigh(1:nda)=d_nneigh(1:nda)
  new_omnn=maxval(trans_nneigh(1:nd))
  !print*,"new_omnn",new_omnn
  !new_omnn=maxval(d_nneigh(1:nd))+10
  !print*,"nneigh1",trans_nneigh(1)
  !trans_neigh=d_trans_neigh
  !trans_dneigh=d_trans_dneigh
  !print*,"neigh1",trans_neigh(1,1:trans_nneigh(1))
  !print*,"new omnn",new_omnn
  !!print*,"neighboring 1.2"
  !  print*,"post_kernel new_omnn",new_omnn,"old omnn",omnn,"mnn_dynam",mnn_dynam,"nd",nd,"nda",nda,"nboxes",nboxes
    !print*,"sizeofneigh",size(neigh,dim=2)
  if(new_omnn/=omnn)then
    omnn=new_omnn
    !print*,"omnn",omnn
    if(allocated(d_neigh)) deallocate(d_neigh)
    if(allocated(d_dneigh)) deallocate(d_dneigh)

    allocate(d_neigh(nda,omnn),d_dneigh(nda,omnn))

    if(allocated(neigh)) deallocate(neigh)
    if(allocated(dneigh)) deallocate(dneigh)

    allocate(neigh(nda,omnn),dneigh(nda,omnn))
  end if
  !print*,"neighboring 1.5"

  call neighbor_transfer_kernel<<<nblocks,nthreads_per_block>>>(d_trans_neigh,d_trans_dneigh,omnn,nd) !>>>>>>GPU ADD 28-11-14
    !print*,"neighboring 1.6",size(neigh,dim=2),size(d_neigh,dim=2)
    
  if(prop_noise>epsilod .or. ffu(1)==1)then !needed for noise and single_node mode
    neigh=d_neigh
    !print*,"neighboring 1.6.1"
    dneigh=d_dneigh
    !print*,"neighboring 1.6.2"
    nneigh=d_nneigh
    !print*,"neighboring 1.6.3"
  end if
    !print*,"neighboring 1.7"
  istat=cudaDeviceSynchronize

  !deallocate(d_trans_neigh,d_trans_dneigh)
    
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!attributes(global) subroutine neighbor_gather_kernel(neigh,dneigh,mnn_dynam,maxbox,nd,maxlen)
!
!integer:: ivv,nbo,iii1,iii2,iii3
!real*8:: ix,iy,iz,dist,nbh
!integer::i,j,ii,jj,k,kk,iii,jjj,kkk
!real*8::a,b,c,d
!
!integer::tipi  !>>Miquel31-12-14
!real*8::dai !,maxlen    !>>Miquel31-12-14
!
!real*8::cx,cy,cz  !>>Miquel6-3-14
!!integer::o                 !>>Miquel6-3-14
!
!integer :: neigh(:,:)!,nneigh(:)  !neighbor matrix
!real*8 :: dneigh(:,:)  !neighbor matrix distances
!integer,value::mnn_dynam,maxbox,nd
!real*8,value::maxlen
!
!
!integer:: lnneigh
!!integer :: lneigh(mnn_dynam) ! ACHTUNG ACHTUNG! ACHTUNG! ACHTUNG! ACHTUNG! ACHTUNG! ACHTUNG! ACHTUNG
!!real*8  :: ldneigh(mnn_dynam)   
!
!integer,shared,dimension(:,:,:,:)::stack(blockDim%x,blockDim%y,blockDim%z,maxbox+1)
!real*8,shared,dimension(:,:,:,:)::dstack(blockDim%x,blockDim%y,blockDim%z,maxbox+1)
!integer,shared,dimension(:,:,:)::nstack(blockDim%x,blockDim%y,blockDim%z)
!integer,shared,dimension(:)::lneigh(500)
!real*8,shared,dimension(:)::ldneigh(500)
!
!integer::stack_count
!integer::nod
!
!  !lneigh=0
!  !lnneigh=0
!
!  !print*,"tall0"  
!  !here each block explores the boxes for one node, each thread explores one box
!  i=threadIdx%x  !thread to stack mapping
!  j=threadIdx%y
!  k=threadIdx%z
!
!  iii=int(real(blockdim%x)*0.5d0)+1
!  ii=i-iii   !setting the index of the neighboring box
!  jj=j-iii
!  kk=k-iii
!  nod=blockIdx%x
!
!  !print*,"blockDims",blockDim%x,blockDim%y,blockDim%z,"bock id",blockIdx%x
!  !print*,"threads",threadIdx%x,threadIdx%y,threadIdx%z
!  !print*,""
!  !if(nod==1)print*,"nod",nod,"ijk",i,j,k,"tall1"
!  !print*,"nod",nod,"ii,jj,kk",ii,jj,kk
!  !print*,"sizeof stack1",size(stack,dim=1)
!  !print*,"sizeof stack2",size(stack,dim=2)
!  !print*,"sizeof stack3",size(stack,dim=3)
!  !print*,"sizeof stack4",size(stack,dim=4)
!  
!  
!  !if (i>nd) goto 567
!
!  !neigh(i,:)=0
!  !dneigh(i,:)=0d0
!  
!  
!  !if (d_rdiffmax<2*d_node(i)%da)then
!  !  nbh=2*d_node(i)%da  !the neighbor search range for diffusion is the same as for node interactions
!  !else
!  !  nbh=d_rdiffmax
!  !end if
!  ix=d_node(nod)%x     ; iy=d_node(nod)%y     ; iz=d_node(nod)%z   
!  iii=nint(iz*d_urv) ; jjj=nint(iy*d_urv) ; kkk=nint(ix*d_urv)
!  ivv=d_node(nod)%altre
!  !print*,"nod",nod,"ijk",i,j,k,"tall1"
!  stack_count=0
!
!  tipi=d_node(nod)%tipus
!  dai=d_node(nod)%da
!  !if(nod==1)print*,"nod",nod,"ijk",i,j,k,"tall2"
!  
!  !if(i==1 .and. j==1 .and. k==1)then !this thread performs this arbitrarily, only needs to be done once per block
!  !  if(tipi<3)then
!  !    stack_count=1
!  !    !if(nod==1)print*,"nod",nod,"ijk",i,j,k,"stack_count ivv",stack_count
!  !    stack(i,j,k,stack_count)=ivv
!  !    dstack(i,j,k,stack_count)=sqrt((d_node(ivv)%x-ix)**2+(d_node(ivv)%y-iy)**2+(d_node(ivv)%z-iz)**2)
!  !  end if
!  !end if
!  !if(nod==1)print*,"nod",nod,"ijk",i,j,k,"tall3"
!
!
!  iii1=iii+ii
!  iii2=jjj+jj
!  iii3=kkk+kk
!  !print*,"iii1",iii1,"iii2",iii2,"iii3",iii3
!  ie=d_boxes(iii1,iii2,iii3)
!  !if(nod==1)print*,"nod",nod,"ijk",i,j,k,"ie first",ie,"stack_count",stack_count
!  do while(ie.ne.0)
!    !print*,"ie",ie
!    if (ie==nod) then ; ie=d_list(ie) ; cycle ; end if
!    if (ie==ivv) then ; ie=d_list(ie) ; cycle ; end if
!    dist=sqrt((d_node(ie)%x-ix)**2+(d_node(ie)%y-iy)**2+(d_node(ie)%z-iz)**2)
!    !>>Miquel31-12-14
!    a=dai+d_node(ie)%da
!    if(tipi>=3)then  
!      if(d_node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
!        if(dist>a)then ; ie=d_list(ie) ; cycle ; end if
!      else  !mesench/ECM vs epithelial
!        b=sqrt(2*(a**2))
!        if(dist>b)then ; ie=d_list(ie) ; cycle ; end if
!      end if
!    else !epithelial vs epithelial
!      b=sqrt(2*(a**2))
!      if(tipi==d_node(ie)%tipus)then !same face epithelials
!        if(b<maxlen) b=maxlen
!        if(dist>b)then ; ie=d_list(ie) ; cycle ; end if
!      else
!        if(dist>b)then ; ie=d_list(ie) ; cycle ; end if
!      end if
!    end if  
!    stack_count=stack_count+1
!    !if(nod==1)print*,"nod",nod,"ijk",i,j,k,"stack_count",stack_count
!    stack(i,j,k,stack_count)=ie
!    dstack(i,j,k,stack_count)=dist!sqrt((d_node(ie)%x-ix)**2+(d_node(ie)%y-iy)**2+(d_node(ie)%z-iz)**2)
!    !lnneigh=lnneigh+1
!    !neigh(lnneigh)=ie
!    !ldneigh(lnneigh)=dist
!    ie=d_list(ie)
!    !if(nod==1)print*,ie,"end of cycle"
!  end do
!  nstack(i,j,k)=stack_count
!  !if(nod==1)print*,"nod",nod,"ijk",i,j,k,"tall4","nstack",nstack(i,j,k),stack_count
!  call syncthreads ()
!
!  if(i==blockDim%x.and.j==blockDim%y.and.k==blockDim%z)then !this thread will put everything together
!    !if(nod==1)print*,nod,"tall5"
!    lnneigh=0
!    if(tipi<3)then
!      lnneigh=1
!      lneigh(1)=ivv
!      ldneigh(1)=sqrt((d_node(ivv)%x-ix)**2+(d_node(ivv)%y-iy)**2+(d_node(ivv)%z-iz)**2)
!    end if
!    do ii=1,blockDim%x
!      do jj=1,blockDim%y
!        do kk=1,blockDim%z
!         !if(nod==1) print*,"nstack",nstack(i,j,k)
!          do iii=1,nstack(ii,jj,kk)
!            ie=stack(ii,jj,kk,iii)
!            !if(nod==1) print*,"ii,jj,kk",ii,jj,kk,"iii",iii,"ie",ie
!            !if(ie==0) exit
!            lnneigh=lnneigh+1
!            lneigh(lnneigh)=ie
!            ldneigh(lnneigh)=dstack(ii,jj,kk,iii)
!          end do
!        end do
!      end do
!    end do
!    d_nneigh(nod)=lnneigh
!    neigh(nod,1:lnneigh)=lneigh(1:lnneigh)
!    dneigh(nod,1:lnneigh)=ldneigh(1:lnneigh)
!    !if(nod==1) print*,nod,"tall6"
!  end if
!end subroutine neighbor_gather_kernel
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!attributes(global) subroutine neighbor_gabriel_kernel(neigh,dneigh,mnn_dynam,screen_radius,nd)
!
!integer:: ivv,ij,ik,jk,iv
!real*8:: ix,iy,iz,dist,udist,nbh
!integer::i,j,ii,jj,k,kk
!real*8::a,b,c,d
!
!integer::tipi  !>>Miquel31-12-14
!real*8::dai !,maxlen    !>>Miquel31-12-14
!
!real*8::jx,jy,jz,cx,cy,cz  !>>Miquel6-3-14
!integer::o                 !>>Miquel6-3-14
!
!integer :: neigh(:,:)!,nneigh(:)  !neighbor matrix
!real*8 :: dneigh(:,:)  !neighbor matrix distances
!integer,value::mnn_dynam,nd
!real*8,value::screen_radius
!
!integer:: lnneigh
!integer :: lneigh(500)!,trans_lneigh(mnn_dynam) ! ACHTUNG ACHTUNG! ACHTUNG! ACHTUNG! ACHTUNG! ACHTUNG! ACHTUNG! ACHTUNG
!real*8  :: ldneigh(500)!,trans_ldneigh(mnn_dynam)
!
!
!
!  !lneigh=0
!  !lnneigh=0
!
!  !ivv=d_node(i)%altre
!  i=blockDim%x*(blockidx%x-1)+threadidx%x  !here each node will be a block
!  if (i>nd) goto 567
!
!
!  ivv=d_node(i)%altre
!  lnneigh=d_nneigh(i)
!  !lneigh(1:mnn_dynam)=neigh(i,1:mnn_dynam)
!  !ldneigh(1:mnn_dynam)=dneigh(i,1:mnn_dynam)
!  !screening by Gabriel graph !>>Miquel6-3-14
!  !a neighbor connection is deleted if the sphere which diameter is the vector connecting the two nodes contains any other node
!  do j=1,lnneigh-1  !ordering the neighors with increasing distance
!    b=dneigh(i,j)
!    ii=0
!    do k=j+1,lnneigh
!      c=dneigh(i,k)
!      if(b>c)then
!        ii=k ; b=dneigh(i,k)
!      end if
!    end do
!    if(ii/=0)then
!      !jj=neigh(i,j)
!      kk=neigh(i,ii) ; c=dneigh(i,ii) !the swap
!      neigh(i,ii)=neigh(i,j) ; dneigh(i,ii)=dneigh(i,j)
!      neigh(i,j)=kk ; dneigh(i,j)=c
!    end if
!  end do
!
!  ix=d_node(i)%x ; iy=d_node(i)%y ; iz=d_node(i)%z
!  !the screening
!  !ii=0 !the number of eliminated connections
!  do j=lnneigh,1,-1
!    jj=neigh(i,j)
!    if(jj==ivv) cycle
!    a=dneigh(i,j)*0.5*screen_radius !the radius of the sphere !>>Miquel28-7-14
!!      jx=d_node(jj)%x ; jy=d_node(jj)%y ; jz=d_node(jj)%z
!    cx=(ix+d_node(jj)%x)*0.5 ; cy=(iy+d_node(jj)%y)*0.5 ; cz=(iz+d_node(jj)%z)*0.5 !the midpoint
!    do k=j-1,1,-1
!      kk=neigh(i,k)
!      d=sqrt((cx-d_node(kk)%x)**2+(cy-d_node(kk)%y)**2+(cz-d_node(kk)%z)**2)
!      if(d-a<d_epsilod)then !there is one node within the sphere, we must delete this connection
!        neigh(i,j)=0
!        !do l=j,lnneigh-ii-1
!        !  neigh(i,l)=neigh(i,l+1)
!        !  dneigh(i,l)=dneigh(i,l+1)
!        !end do
!        !neigh(i,lnneigh-ii)=0
!        !dneigh(i,lnneigh-ii)=0
!        !ii=ii+1
!        exit
!      end if
!    end do
!  end do
!  ii=0
!  do j=1,lnneigh
!    jj=neigh(i,j)
!    if(jj/=0)then
!      ii=ii+1
!      lneigh(ii)=jj ; ldneigh(ii)=dneigh(i,j)
!    end if
!  end do
!  !lnneigh=ii
!
!
!  !if(i==nd) print*,"i",i,"lnneigh",lnneigh,"ii",ii
!  !lnneigh=snneigh
!  !if(i==nd) print*,"i",i,"lnneigh post",lnneigh
!  
!  d_nneigh(i)=ii
!  neigh(i,1:ii)=lneigh(1:ii)
!  dneigh(i,1:ii)=ldneigh(1:ii)
!567 continue
!  
!end subroutine neighbor_gabriel_kernel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



attributes(global) subroutine neighbor_build_kernel(neigh,dneigh,nneigh,mnn_dynam,maxbox,nd)

integer:: ivv,ii1,ii2,ii3,nbo,iii1,iii2,iii3,ie,ierr,ij,ik,jk,iv
real*8:: ix,iy,iz,dist,udist,nbh
integer::i,j,ii,jj,k,kk
real*8::a,b,c,d

integer::tipi  !>>Miquel31-12-14
real*8::dai !,maxlen    !>>Miquel31-12-14

!!!triangulation variables
!integer:: npt !number of points
!integer:: sizht !size of hash table
!integer:: maxbf,maxfc,nbf,nfc,nface,ntetra   !size of arrays
!integer,allocatable :: vecinod(:,:),vecic(:)!miguel
!real*8,allocatable :: dvecinod(:,:)
!real*8,allocatable :: vcl(:,:) !point coordinates
!integer,allocatable :: vm(:),ht(:),bf(:,:),fc(:,:) !point indices (the algorithm reorders)
!integer,dimension(:) :: border(nd)
real*8::jx,jy,jz,cx,cy,cz  !>>Miquel6-3-14
integer::o                 !>>Miquel6-3-14

integer :: neigh(:,:),nneigh(:)  !neighbor matrix
real*8 :: dneigh(:,:)  !neighbor matrix distances
integer,value::mnn_dynam,maxbox,nd,mnn


integer:: lnneigh
!integer :: lneigh(mnn_dynam) ! ACHTUNG ACHTUNG! ACHTUNG! ACHTUNG! ACHTUNG! ACHTUNG! ACHTUNG! ACHTUNG
!real*8  :: ldneigh(mnn_dynam)   
integer :: lneigh(500) ! ACHTUNG ACHTUNG! ACHTUNG! ACHTUNG! ACHTUNG! ACHTUNG! ACHTUNG! ACHTUNG
real*8  :: ldneigh(500)   



  lneigh=0
  lnneigh=0

  i=blockDim%x*(blockidx%x-1)+threadidx%x  !here each node will be a block
  if (i>nd) goto 567

  !neigh(i,:)=0
  !dneigh(i,:)=0d0
  
  
  if (d_rdiffmax<2*d_node(i)%da)then
    nbh=2*d_node(i)%da  !the neighbor search range for diffusion is the same as for node interactions
  else
    nbh=d_rdiffmax
  end if
  ix=d_node(i)%x     ; iy=d_node(i)%y     ; iz=d_node(i)%z   
  ii1=nint(iz*d_urv) ; ii2=nint(iy*d_urv) ; ii3=nint(ix*d_urv)
  ivv=d_node(i)%altre

  nbo=nint(nbh*d_urv) !;print*,"nbo",nbo,"rdiffmax",rdiffmax
  lnneigh=0 !; dif_nneigh(i)=0

  if(d_node(i)%tipus<3)then
    lnneigh=1
    lneigh(1)=ivv
    ldneigh(1)=sqrt((d_node(ivv)%x-ix)**2+(d_node(ivv)%y-iy)**2+(d_node(ivv)%z-iz)**2)
  end if

  tipi=d_node(i)%tipus
  dai=d_node(i)%da
  do i1=-nbo,nbo 
    iii1=ii1+i1
    do i2=-nbo,nbo
      iii2=ii2+i2
      do i3=-nbo,nbo
        iii3=ii3+i3
        ie=d_boxes(iii3,iii2,iii1)
        do while(ie.ne.0)
          if (ie==i) then ; ie=d_list(ie) ; cycle ; end if
          if (ie==ivv) then ; ie=d_list(ie) ; cycle ; end if
          dist=sqrt((d_node(ie)%x-ix)**2+(d_node(ie)%y-iy)**2+(d_node(ie)%z-iz)**2)
          !>>Miquel31-12-14
          a=dai+d_node(ie)%da
          if(tipi>=3)then  
            if(d_node(ie)%tipus>=3)then !mesench/ECM vs mesench/ECM
              if(dist>a)then ; ie=d_list(ie) ; cycle ; end if
            else  !mesench/ECM vs epithelial
              b=sqrt(2*(a**2))
              if(dist>b)then ; ie=d_list(ie) ; cycle ; end if
            end if
          else !epithelial vs epithelial
            b=sqrt(2*(a**2))
            if(tipi==d_node(ie)%tipus)then !same face epithelials
              if(b<d_maxlen) b=d_maxlen
              if(dist>b)then ; ie=d_list(ie) ; cycle ; end if
            else
              if(dist>b)then ; ie=d_list(ie) ; cycle ; end if
            end if
          end if  
          !>>Miquel31-12-14
          lnneigh=lnneigh+1
          lneigh(lnneigh)=ie
          ldneigh(lnneigh)=dist
          ie=d_list(ie)
        end do
      end do
    end do
  end do

  if(d_ffu(3)==1)then
    !screening by Gabriel graph !>>Miquel6-3-14
    !a neighbor connection is deleted if the sphere which diameter is the vector connecting the two nodes contains any other node
    do j=1,lnneigh-1  !ordering the neighors with increasing distance
      b=ldneigh(j)
      ii=0
      do k=j+1,lnneigh
        c=ldneigh(k)
        if(b>c)then
          ii=k ; b=ldneigh(k)
        end if
      end do
      if(ii/=0)then
        !jj=neigh(i,j)
        kk=lneigh(ii) ; c=ldneigh(ii) !the swap
        lneigh(ii)=lneigh(j) ; ldneigh(ii)=ldneigh(j)
        lneigh(j)=kk ; ldneigh(j)=c
      end if
    end do
  
    !the screening
    ii=0 !the number of eliminated connections
    do j=lnneigh,1,-1
      jj=lneigh(j)
      if(jj==ivv) cycle
      a=ldneigh(j)*0.5*d_screen_radius !the radius of the sphere !>>Miquel28-7-14
!      jx=d_node(jj)%x ; jy=d_node(jj)%y ; jz=d_node(jj)%z
      cx=(ix+d_node(jj)%x)*0.5 ; cy=(iy+d_node(jj)%y)*0.5 ; cz=(iz+d_node(jj)%z)*0.5 !the midpoint
      do k=j-1,1,-1
        kk=lneigh(k)
        d=sqrt((cx-d_node(kk)%x)**2+(cy-d_node(kk)%y)**2+(cz-d_node(kk)%z)**2)
        if(d-a<d_epsilod)then !there is one node within the sphere, we must delete this connection
          do l=j,lnneigh-ii-1
            lneigh(l)=lneigh(l+1)
            ldneigh(l)=ldneigh(l+1)
          end do
          lneigh(lnneigh-ii)=0
          ldneigh(lnneigh-ii)=0
          ii=ii+1
          exit
        end if
      end do
    end do
    lnneigh=lnneigh-ii
  end if

  
  nneigh(i)=lnneigh
  neigh(i,1:lnneigh)=lneigh(1:lnneigh)
  dneigh(i,1:lnneigh)=ldneigh(1:lnneigh)
567 continue
  
end subroutine neighbor_build_kernel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

attributes(global) subroutine neighbor_transfer_kernel(neigh1,dneigh1,omnn,nd)
integer :: neigh1(:,:)!,neigh2(:,:)  !neighbor matrix
real*8 :: dneigh1(:,:)!,dneigh2(:,:)  !neighbor matrix distances
integer,value::omnn,nd
integer::i,j


  i=blockDim%x*(blockidx%x-1)+threadidx%x  !here each node will be a block
  if (i>nd) goto 567
  d_neigh(i,1:omnn)=neigh1(i,1:omnn)
  d_dneigh(i,1:omnn)=dneigh1(i,1:omnn)
567 continue
end subroutine

!!!!!FORCES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine dev_iterdiferencial
integer::nodmo,i,j,k,ii
real*8::a,b,c
real*8::ox,oy,oz !miguel4-1-13
integer::istat
real*8,device,dimension(:,:)::d_pxpypz(nd,5)
real*8,pinned,allocatable::pxpypz(:,:)
real*8,pinned,allocatable::cel_pol(:,:)
real*8,device,dimension(:,:)::d_cel_pol(nd,3)

  !CALCULATING FORCES AND MOVEMENT VECTORS

  !print*,"debug:: pre calling forces_kernel",getot
  d_pxpypz=0d0
  
  if(npag(nparam_per_node+16)>0) then
    allocate(cel_pol(nd,3))
    do i=1,nd
      k=node(i)%icel
      if(k>0)then
        cel_pol(i,1)=cels(k)%polx
        cel_pol(i,2)=cels(k)%poly
        cel_pol(i,3)=cels(k)%polz
      end if
    end do
    d_cel_pol=cel_pol
  end if
  
  
  call forces_kernel<<<nblocks,nthreads_per_block>>>(d_pxpypz,k_bu,d_cel_pol)
  if(npag(nparam_per_node+16)>0) then
    deallocate(cel_pol)
  end if

!  call forces_kernel<<<nd,1>>>
  !call cudaDeviceSynchronize
  !istat=cudaDeviceSynchronize
  !print*,"debug:: forces_kernel done, now copying back",getot
  allocate(pxpypz(nd,5))
  pxpypz=d_pxpypz
  !print*,"debug:: pxpypz copyed to host",getot
  px(1:nd)=pxpypz(1:nd,1)
  py(1:nd)=pxpypz(1:nd,2)
  pz(1:nd)=pxpypz(1:nd,3)
  dex(1:nd)=pxpypz(1:nd,4)
  fmeanl(1:nd)=pxpypz(1:nd,5)
  !print*,getot,"maxval fmeanl",maxval(fmeanl(1:nd))
  !print*,"check*** p18",px(18),py(18),pz(18)
  !print*,"check*** p19",px(19),py(19),pz(19)

  !print*,"check*** p2",px(2),py(2),pz(2)
  !print*,"check*** p3",px(3),py(3),pz(3)
  !print*,"check*** dex",dex(1)
  !print*,"check*** fmeanl",fmeanl(1)

  
  
  
  
  !print*,"debug:: forces data copyed to host",getot
  a=0.0d0
  do i=1,nd
    if (dex(i)>a) then ; a=dex(i) ; ii=i ; end if
  end do
  
  if (a<epsilod) then ; rtime=rtime+delta ; delta=deltamin; return ; end if  !>>> Is 29-8-13
  if(ffu(12)==0)then
    delta=resmax/a      !delta is adjusted so the node that has the greatest force
                      !will always move the same absolute distance, and the rest will
                      !move proportionally to that force
  else
    delta=deltamin
!    delta=deltamax 
  end if
                      
  if (delta>deltamax) delta=deltamax
  if (delta<deltamin) delta=deltamin
  
  deallocate(pxpypz)

end subroutine dev_iterdiferencial

!*************************************************************************

attributes(global) subroutine forces_kernel (pxpypz,k_bu,cel_pol)

real*8   ::ix,iy,iz
real*8   ::a,b,c,d,e,f,g
integer  ::celi,celj,nod
real*8   ::youe,repe,adhe,adho,repcele,deqe,ideqe
real*8   ::younod,repnod,adhnod,repcelnod,reqnod,tornod,stornod    !>>>> Miquel 16-8-13
real*8   ::ax,ay,az,bx,by,bz
real*8   ::ud,udd,uddd,dd,ddd
real*8   ::cx,cy,cz,ccx,ccy,ccz,dotp,pesco
real*8   ::icx,icy,icz,idd,iudd,id !>>>>>>>> MIQUEL 4-3-13
real*8   ::mcx,mcy,mcz,md,umd					!>>>>>>>> MIQUEL 30-4-13
real*8   ::nodda,posca
integer  ::ivv		!>>>>>>>> MIQUEL 4-3-13
integer  ::i,j,k,ii,jj,kk,ic,iii,jjj,kkk,iiii,jjjj,kkkk,iv,kjjj,jkkk
integer  ::nuve,inuve
integer  ::tipi,tipic																!>>>>>>>>>>>>>>>>>>>>>Miquel 23-4-13
integer  ::switch,twoep
real*8   ::aa

integer  ::twomes !!!!ECTO MOD !>>Miquel26-5-15

integer  ::lateral,vertical !flags that tell if there is a lateral or vertical component to take into account !>>Miquel28-1-14

real*8   ::rvx,rvy,rvz   !the resulting force vector
real*8   ::uvx,uvy,uvz   !unit vector
real*8   ::pox,poy,poz   !polarisation vector (from the cell)

real*8   ::ad,fd  !>>Miquel28-1-14

!real*8   ::upr(d_nd)
!integer  ::iupr(d_nd)
!real*8   ::r(d_nd),er(d_nd)

real*8   ::rcilx,rcily,rcilz
real*8   ::rtorx,rtory,rtorz
real*8   ::rstorx,rstory,rstorz
real*8   ::rsprx,rspry,rsprz

real*8   ::pxpypz(:,:)
real*8   ::cel_pol(:,:)
real*8,value::k_bu

!integer, parameter :: d_mnn=200  !declared in general module !>>Miquel24-2-14

                              ! ACHTUNG, WE ASSUME THAT THE MAXIMAL NUMBER OF NODES INTERACTING WITH A NODE IS d_mnn, otherwise these 
                              ! matrices become unbareably large for large systems
!integer  ::vei(d_mnn,d_nd)        ! force storing arrays for optimization    >>>>>>Miquel 7-8-16
!real*8   ::vforce(d_mnn,d_nd)     ! so we don't calculate each interaction twice
!real*8   ::vuvx(d_mnn,d_nd)       ! (for now we do it only for cylinder and normal interactions,
!real*8   ::vuvy(d_mnn,d_nd)       ! not torsion nor springs)
!real*8   ::vuvz(d_mnn,d_nd)       !
!real*8   ::vtorforce(d_mnn,d_nd)  ! NOTE THE NON-CONVENTIONAL FORM OF THE ARRAYS:
!real*8   ::vtoruvx(d_mnn,d_nd)    ! DIM 1 IS THE LISTS THE NEIGHBORS,
!real*8   ::vtoruvy(d_mnn,d_nd)    ! DIM 2 LISTS THE NODE ID
!real*8   ::vtoruvz(d_mnn,d_nd)    !
!real*8   ::vstorforce(d_mnn,d_nd) !
!real*8   ::vstoruvx(d_mnn,d_nd)   !
!real*8   ::vstoruvy(d_mnn,d_nd)   !
!real*8   ::vstoruvz(d_mnn,d_nd)   !
!integer ::nveins(d_nd)         !
integer  ::epinveins!(d_nd)      ! we store how many same-side epithelial neighbors an epithelial node has !>>>Miquel4-4-14
integer  ::alone              ! 0 if the node is really alone

!  vcilx=0 ; vcily=0 ; vcilz=0 ; vsprx=0     !force vectors storage matrices, for different components
!  vtorx=0 ; vtory=0 ; vtorz=0 ; vspry=0     !of the resulting force
!  vstorx=0 ; vstory=0 ; vstorz=0 ; vsprz=0
!  vei=0 ; vforce=0 ; vuvx=0 ; vuvy=0 ; vuvz=0 ; nveins=0  !>>>>>>Miquel 7-8-16
!  vtorforce=0 ; vtoruvx=0 ; vtoruvy=0 ; vtoruvz=0         !
!  vstorforce=0 ; vstoruvx=0 ; vstoruvy=0 ; vstoruvz=0     !

  rcilx=0d0  ; rcily=0d0  ; rcilz=0d0  !they store the force components for all the nodes !>>Miquel4-4-14
  rtorx=0d0  ; rtory=0d0  ; rtorz=0d0
  rstorx=0d0 ; rstory=0d0 ; rstorz=0d0
  epinveins=0
  
    !nod=blockidx%x  !here each node will be a block
    nod=blockDim%x*(blockidx%x-1)+threadidx%x
    if (nod>d_nd) goto 987

    !d_fmeanl(nod)=0 !; d_fmeanv=0  !storage vector that makes the balance between compressive and tensile forces within a node    !>>>Miquel23-1-14
 
    if (d_node(nod)%hold==2) then                ! >>> Is 30-6-14
      pxpypz(nod,:)=0  !module of the force vector ! >>> Is 30-6-14
      !d_px(nod)=0 ; d_py(nod)=0 ; d_pz(nod)=0       ! >>> Is 30-6-14
      goto 987 !>>GPU ADD 2-12-14
    end if                                    ! >>> Is 30-6-14

    rsprx=0.0d0  ; rspry=0.0d0; rsprz=0.0d0
    ix=d_node(nod)%x ; iy=d_node(nod)%y ; iz=d_node(nod)%z
    tipi=d_node(nod)%tipus ; celi=d_node(nod)%icel
    !iii=nint(ix*urv)    ; jjj=nint(iy*urv)   ; kkk=nint(iz*urv)	!>>GPU REMOVE 2-12-14
    rvx=0d0    ; rvy=0d0    ; rvz=0d0
    nuve=0!nuve=nveins(nod) !>>>>>Miquel 7-8-13 : this is not always zero since we fill this matrix from its neighbours
    switch=0
    alone=0

    !SPRINGS
    if(tipi<3)then
      iv=d_node(nod)%altre
      ax=d_node(iv)%x   ; ay=d_node(iv)%y    ; az=d_node(iv)%z
      cx=ax-ix        ; cy=ay-iy         ; cz=az-iz
      dd=sqrt(cx**2+cy**2+cz**2)
      udd=1d0/dd
      !if(d_ffu(1)==1)then
        uvx=cx*udd ; uvy=cy*udd ; uvz=cz*udd  !the unit vector
        ddd=dd-d_node(nod)%reqs-d_node(iv)%reqs  !>>Miquel5-2-14
        f=2*d_node(nod)%ke*ddd    !the force
        rsprx=f*uvx ; rspry=f*uvy ; rsprz=f*uvz
        !fmeanv(nod)=ddd !>>>Miquel5-2-14
      !end if
    else
      iv=0    !>>>>Miquel17-1-14
    end if
    younod=d_node(nod)%you                                            !>>>>>Miquel 16-8-13
    repnod=d_node(nod)%rep                                            !
    adhnod=d_node(nod)%adh !default inespecific adhesion of node nodmo!
    repcelnod=d_node(nod)%repcel                                      !
    tornod=d_node(nod)%tor                                            !
    stornod=d_node(nod)%stor                                          !
    reqnod=d_node(nod)%req                                            !
    nodda=d_node(nod)%da

    !NODE's REPULSIONS AND ADHESIONS
    !if(nod==275) print*,"nneigh",d_nneigh(nod)

    do i=1,d_nneigh(nod)
      ic=d_neigh(nod,i)
      !if(nod==275) print*,"neigh de 275",d_neigh(nod,i)
      !if(ic<nod)then !we have calculated that interaction already !>>Miquel4-4-14
      !  alone=1
      !  cycle
      !end if

      !so it turns out that the nod-ic interactions has not been calculated before

      bx=d_node(ic)%x   ; by=d_node(ic)%y    ; bz=d_node(ic)%z
      ccx=bx-ix       ; ccy=by-iy        ; ccz=bz-iz		
      !d=sqrt(ccx**2+ccy**2+ccz**2)
      d=d_dneigh(nod,i)
      ud=1d0/d
      tipic=d_node(ic)%tipus
      twoep=0
      ad=0 ; fd=0 ; ddd=0     !>>Miquel28-1-14
      
      twomes=0 !!ECTO MOD

      alone=1      !this is crappy but fast, it makes that lonely nodes are eliminated in squares
      if (tipi<3)then
        if(tipic<3)then
          ivv=d_node(ic)%altre
          icx=d_node(ivv)%x-bx ; icy=d_node(ivv)%y-by ; icz=d_node(ivv)%z-bz 
          idd=sqrt(icx**2+icy**2+icz**2)      ; iudd=1d0/idd	!ic's spring vector	
          posca=icx*cx+icy*cy+icz*cz
          if (tipi==tipic) then           ! BOTH NODES ARE EPITHELIAL IN THE SAME SIDE we are equal so we must have lateral adhesion
            if (posca>d_epsilod) then    !SO WE THE NEIGHBOUR IS IN A CONTIGUOUS PART OF THE EPITHELIUM
              mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=1d0/sqrt(mcx**2+mcy**2+mcz**2)  !we take as spring vector the sum of ic's and nod's spring vectors
              dotp=(mcx*ccx+mcy*ccy+mcz*ccz)
              ddd=abs(dotp)*md  !vertical component
              !a=nodda+d_node(ic)%da
              ad=d**2-ddd**2 ; if(ad<d_epsilod) ad=d_epsilod
              ad=sqrt(ad)              !lateral component
              if (ad-nodda-d_node(ic)%da>d_epsilod) cycle
              pesco=dotp*md**2
              uvx=(ccx-mcx*pesco) ; uvy=(ccy-mcy*pesco) ; uvz=(ccz-mcz*pesco)  !unit vector of the within cilinder force !el modul és el mateix que el vector c
              !a=1d0/sqrt(uvx**2+uvy**2+uvz**2)
              a=1/ad
              uvx=uvx*a ; uvy=uvy*a ; uvz=uvz*a
              fd=ad  !fd as the distance we will use to assert the force modulus  !>>Miquel28-1-14
              twoep=1
              epinveins=epinveins+1 !>>>GPU CHANGED !>>Miquel2-12-14
              !epinveins(ic)=epinveins(ic)+1
            else
              ! apical/basal contact, two epithelia from the same side
              if(d_ffu(16)==0)then
                mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=sqrt(mcx**2+mcy**2+mcz**2)  !we take as spring vector the sum of ic's and nod's spring vectors
                if(md<d_epsilod)then !the two cylinders are parallel, the vector used is the spring !>>Miquel20-8-14
                  dotp=(cx*ccx+cy*ccy+cz*ccz)
                  if (dotp<d_epsilod) then         ! projection of the vector from nod to ic into the vector from nod to iv
                    ddd=-dotp*udd                 ! that is the distance UP in the direction of altre
                    a=nodda+d_node(ic)%da
                    if (ddd-a<d_epsilod) then
                      ddd=d**2-ad**2 ;if(ddd<d_epsilod) ddd=d_epsilod
                      ddd=sqrt(ddd)              !lateral component
                      if (ad-a>d_epsilod) cycle
                      uvx=-cx*udd ; uvy=-cy*udd ; uvz=-cz*udd  !unit vector vertical
                      fd=ddd     !distance used, vertical component
                      twoep=2
                      goto 300
                    end if
                    cycle 
                  end if  
                else  !proper apical cylindric interface applied !>>Miquel20-8-14
                  md=1d0/md
                  dotp=(mcx*ccx+mcy*ccy+mcz*ccz)
                  ad=abs(dotp)*md  !vertical component
                  a=nodda+d_node(ic)%da
                  if(ad-a<d_epsilod)then
                    ddd=d**2-ad**2 ;if(ddd<d_epsilod) ddd=d_epsilod
                    ddd=sqrt(ddd)              !lateral component
                    if (ddd-a>d_epsilod) cycle
                    pesco=dotp*md**2
                    uvx=(ccx-mcx*pesco) ; uvy=(ccy-mcy*pesco) ; uvz=(ccz-mcz*pesco)  !unit vector of the within cilinder force !el modul és el mateix que el vector c
                    a=1/ddd
                    uvx=uvx*a ; uvy=uvy*a ; uvz=uvz*a
                    fd=ddd     !distance used, vertical component
                    twoep=2
                    goto 300
                  end if
                  cycle 
                end if  
              else  !apical-apical contact from the same face, sphere-sphere interface used !>>Miquel20-8-14
                if(d-nodda-d_node(ic)%da<d_epsilod)then
                  uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud
                  fd=d
                  twoep=2
                  goto 300
                end if
                cycle 
              end if  
            end if
          else
            if (posca<0.0d0) then           
              ! we are from contrary sides so we must have face adhesion
              ! this is adhesion and repulsion between epithelial sides
              dotp=(cx*ccx+cy*ccy+cz*ccz)                     
              if (dotp<d_epsilod) then         ! projection of the vector from nod to ic into the vector from nod to iv
                ddd=-dotp*udd                 ! that is the distance UP in the direction of altre
                a=nodda+d_node(ic)%da
                if (ddd-a<d_epsilod) then
                  ad=d**2-ddd**2 ;if(ad<d_epsilod) ad=d_epsilod
                  ad=sqrt(ad)              !lateral component
                  if (ad-a>d_epsilod) cycle
                  uvx=cx*udd ; uvy=cy*udd ; uvz=cz*udd  !unit vector vertical
                  fd=ddd     !distance used, vertical component
                  twoep=2
                  goto 300
                end if
              end if
            end if
            cycle
          end if
        else
          ! nod IS EPITHELIAL and ic is not
          ! we check the distance to radial distance in the plane of the ic cylinder 
          dotp=(cx*ccx+cy*ccy+cz*ccz)!;print*,"dotp epi-mes",dotp,nod,ic !hi ha un vector del revés, per tant això està al revés també
          if (dotp<0.0) then
            ddd=abs(dotp*udd)    !vertical component
            a=nodda+d_node(ic)%da
            if (ddd-a<d_epsilod) then
              ad=d**2-ddd**2 ;if(ad<d_epsilod) ad=d_epsilod
              ad=sqrt(ad)              !lateral component
              if (ad-a>d_epsilod) cycle
              !uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud  
              !fd=d  !now interactions between epith. and mesench. work in a sphere-sphere manner !>>Miquel21-2-14
              
              uvx=-cx*udd ; uvy=-cy*udd ; uvz=-cz*udd  !unit vector vertical
              fd=ddd     !distance used, vertical component
              twomes=1
            else
              cycle
            end if
          else
            cycle
          end if
        end if
      else
        if(tipic<3)then          ! IC IS EPITHELIAL and nod is not
          ivv=d_node(ic)%altre     ! we check the distance to radial distance in the plane of the ic cylinder 
          icx=d_node(ivv)%x-bx ; icy=d_node(ivv)%y-by ; icz=d_node(ivv)%z-bz   
          idd=1d0/sqrt(icx**2+icy**2+icz**2)
          dotp=(icx*ccx+icy*ccy+icz*ccz)  !hi ha un vector del revés, per tant això està al revés també
          if (dotp>0.0) then
            ddd=dotp*idd        !vertical component
            a=nodda+d_node(ic)%da
            if (ddd-a<0.0) then
              ad=d**2-ddd**2 ;if(ad<d_epsilod) ad=d_epsilod
              ad=sqrt(ad)              !lateral component
              if (ad-a>d_epsilod) cycle
              !uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud
              !fd=d  !now interactions between epith. and mesench. work in a sphere-sphere manner !>>Miquel21-2-14
              
              uvx=icx*idd ; uvy=icy*idd ; uvz=icz*idd  !unit vector vertical
              fd=ddd     !distance used, vertical component
              twomes=1
            else
              cycle
            end if
          else
            cycle
          end if
        else
          fd=d !BOTH NODES ARE NON-EPITHELIAL: we just consider the interactions between nodes as such
          if (fd-nodda-d_node(ic)%da>d_epsilod) cycle
          uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud  !unit vector of the within cilinder force !el modul és el mateix que el vector c
          twomes=2
        end if
      end if

300   nuve=nuve+1

      !if (nuve>d_mnn) then; 
      !  print *,"PANIC!!! PANIC!!!PANIC!!!; this is going to crash because there is too many " ; 
      !  print *,"neighbors per a node to interact, and vuvx matrix and those are too small"
      !end if 
  
       !ALL THAT WAS JUST TO CALCULATE THE RIGHT DISTANCE BETWEEN NODES, ddd, NOW WE CALCULATE THE ACTUAL ENERGIES
      !if (d_ffu(2)==1) then IS 25-12-13
        if(d_node(nod)%icel==d_node(ic)%icel)then
          deqe=(reqnod+d_node(ic)%req)             !>>>> Miquel 16-8-13
          if(fd-deqe<-d_epsilod)then 				
            repe=0.5d0*(repnod+d_node(ic)%rep)     !>>>> Miquel 16-8-13
            f=2*repe*(fd-deqe)
          else
            youe=(younod+d_node(ic)%you)          !>>>> Miquel 16-8-13
            f=youe*(fd-deqe)
          end if
        else
          deqe=(reqnod+d_node(ic)%req)          !>>>> Miquel 16-8-13
          if(fd-deqe<-d_epsilod)then 
            repcele=(repcelnod+d_node(ic)%repcel) !>>>> Miquel 16-8-13
            f=repcele*(fd-deqe)       !in fact that is the derivative of the energy
          else
            !if(twoep/=2)then   !!SCALE MOD. prevent apical sides of the epith. to stick when the epith. is folded
              adhe=0.5d0*(adhnod+d_node(ic)%adh)          !>>>> Miquel 16-8-13
              if (d_npag(1)>0) then ! we have adhesion molecules
                do j=1,d_npag(1)
                  k=d_whonpag(1,j)
                  do kjjj=1,d_npag(1)
                    jkkk=d_whonpag(1,kjjj)
                    if (d_gex(nod,k)>0.0d0.and.d_gex(ic,jkkk)>0.0d0) then     
                      adhe=adhe+d_gex(nod,k)*d_gex(ic,jkkk)*d_kadh(int(d_gen(k)%wa(1)),int(d_gen(jkkk)%wa(1))) ! >>> Is 7-6-14    !this is specific adhesion
                    end if
                  end do
                end do
              end if
              f=2*adhe*(fd-deqe)          !in fact that is the derivative of the energy
            !end if
          end if
        end if !;print*,"f",f,"fd",fd
        !vforce(nuve,nod)=f     !>>>>>Miquel 7-8-13
        !if (twoep/=2) vforce(jjjj,ic)=f      !>>>>>Miquel 7-8-13
        rcilx=rcilx+f*uvx ; rcily=rcily+f*uvy ; rcilz=rcilz+f*uvz  !>>GPU CHANGED 2-12-14
        !rcilx(ic)=rcilx(ic)-f*uvx ; rcily(ic)=rcily(ic)-f*uvy ; rcilz(ic)=rcilz(ic)-f*uvz  !>>>Miquel 4-4-14
        !if(tipi<3 .and. tipic<3)then 
           
        pxpypz(nod,5)=pxpypz(nod,5)+(fd-deqe) !; fmeanl(ic)=fmeanl(ic)+(fd-deqe) ! >>> Is 21-6-14 !>>GPU CHANGED 2-12-14

        !endif  !tensile forces positive !>>Miquel23-1-14
        !er(nuve)=f
      !end if


      !TORSION
      if(d_ffu(4)==0 .and.twoep==1 ) then !it is only between epithelial nodes

        !surface tension-like torsion (original)
        mcx=icx+cx; mcy=icy+cy; mcz=icz+cz; md=1d0/sqrt(mcx**2+mcy**2+mcz**2)  !we take as spring vector the sum of ic's and nod's spring vectors
        dotp=(mcx*ccx+mcy*ccy+mcz*ccz)*md !vertical projection, more stable than the angle
        !cal fer la mitja de les 2 molles, la d'ic i la de nod, perquè despres no hi hagi assimetries
        if(abs(dotp)-d_angletor*d>d_epsilod)then
          uvx=mcx*md ; uvy=mcy*md ; uvz=mcz*md
          f=(tornod+d_node(ic)%tor)*dotp         !>>>>>Miquel 7-8-13
          rstorx=rstorx+f*uvx ; rstory=rstory+f*uvy ; rstorz=rstorz+f*uvz !>>GPU CHANGED 2-12-14

          !parallel springs torsion   !new  !>>Miquel14-3-14
          dotp=(cx*ccx+cy*ccy+cz*ccz)*udd !vertical projection, more stable than the angle
          f=(tornod+d_node(ic)%tor)*dotp
          uvx=ccx*ud ; uvy=ccy*ud ; uvz=ccz*ud
          rtorx=rtorx+f*uvx ; rtory=rtory+f*uvy ; rtorz=rtorz+f*uvz
        end if
      end if
      
 

      !!! ECTO MOD: implementation of single node cell migration (TOOOOOORQUE)
      ii=d_nparam_per_node+16
      if(d_npag(ii)>0 .and. twomes>0 ) then !it is only between epithelial nodes
        !nod crawls over ic
        if(tipi/=3) goto 587
        mcx=cel_pol(nod,1) ; mcy=cel_pol(nod,2) ; mcz=cel_pol(nod,3) !unit vector
        !if(nod==15) print*,"nod",nod,"pol vec",mcx,mcy,mcz
        if(mcx==0 .and. mcy==0 .and. mcz==0) goto 587 !no polarization, no migration
        dotp=ccx*mcx+ccy*mcy+ccz*mcz
        if(dotp>d_epsilod)then !angle between neighbor and polarization vector < 90º
          pesco=dotp*ud**2
          a=0d0
          do k=1,d_npag(ii) ; 
            kk=d_whonpag(ii,k) 
            if (d_gex(nod,kk)>0.0d0) then ; 
              a=a+d_gex(nod,kk)*d_gen(kk)%wa(ii)
            endif 
          enddo
          f=dotp*ud*a
          uvx=mcx-ccx*pesco ; uvy=mcy-ccy*pesco ; uvz=mcz-ccz*pesco
          aa=sqrt(uvx**2+uvy**2+uvz**2)
          if(aa>epsilod)then ; b=1/sqrt(aa) ; else ; b=0d0 ; end if
          uvx=f*uvx*b ; uvy=f*uvy*b ; uvz=f*uvz*b
          !if(nod==15) print*,"nod",nod,"uv",uvx,uvy,uvz,"f",f,"dotp",dotp
          !print*,"nod",nod,"ic",ic,"dotp",dotp
          rtorx=rtorx+uvx ; rtory=rtory+uvy ; rtorz=rtorz+uvz
          !rtorx(ic)=rtorx(ic)-uvx ; rtory(ic)=rtory(ic)-uvy ; rtorz(ic)=rtorz(ic)-uvz
        end if
        !ic crawls over nod
587     if(tipic/=3) goto 324
        mcx=cel_pol(ic,1) ; mcy=cel_pol(ic,2) ; mcz=cel_pol(ic,3) !unit vector
        if(mcx==0 .and. mcy==0 .and. mcz==0) goto 324 !no polarization, no migration
        dotp=-ccx*mcx-ccy*mcy-ccz*mcz
        if(dotp>d_epsilod)then !angle between neighbor and polarization vector < 90º
          pesco=dotp*ud**2
          a=0d0
          do k=1,d_npag(ii) ; 
            kk=d_whonpag(ii,k) 
            if (d_gex(ic,kk)>0.0d0) then ; 
              a=a+d_gex(ic,kk)*d_gen(kk)%wa(ii)
            endif 
          enddo
          f=dotp*ud*a
          uvx=mcx+ccx*pesco ; uvy=mcy+ccy*pesco ; uvz=mcz+ccz*pesco
          aa=sqrt(uvx**2+uvy**2+uvz**2)
          if(aa>epsilod)then ; b=1/sqrt(aa) ; else ; b=0d0 ; end if
          uvx=f*uvx*b ; uvy=f*uvy*b ; uvz=f*uvz*b
          !print*,"nod",nod,"ic",ic,"dotp",dotp
          !rtorx(ic)=rtorx(ic)+uvx ; rtory(ic)=rtory(ic)+uvy ; rtorz(ic)=rtorz(ic)+uvz
          rtorx=rtorx-uvx ; rtory=rtory-uvy ; rtorz=rtorz-uvz
        end if
      end if




324   continue
    end do

    if(tipi<3)then
      if(epinveins>0)then                  !>>>Miquel24-3-14
        !fmeanl(nod)=0
      !else
        pxpypz(nod,5)=pxpypz(nod,5)/epinveins
      end if
    else
      if(nuve>0)then                  !>>>Miquel24-3-14
        !fmeanl(nod)=0
      !else
        pxpypz(nod,5)=pxpypz(nod,5)/nuve
      end if
    end if

789 if (d_ffu(8)==1) then   !this is kind of crappy in the sense that 
      if (alone==0) then
        d_node(nod)%talone=d_node(nod)%talone+1
      else
        d_node(nod)%talone=0
      end if
    end if

    !if (aut==0) then ! if aut==0 there is not display   
    !  vsprx(nod)=rsprx ; vspry(nod)=rspry ; vsprz(nod)=rsprz                     !putting the force components into storage vectors for display
    !  vcilx(nod)=rcilx(nod)    ; vcily(nod)=rcily(nod)    ; vcilz(nod)=rcilz(nod)         !this part should be turned down
    !  vtorx(nod)=rtorx(nod)    ; vtory(nod)=rtory(nod)    ; vtorz(nod)=rtorz(nod)         !when there is no display
    !  vstorx(nod)=rstorx(nod)  ; vstory(nod)=rstory(nod)  ; vstorz(nod)=rstorz(nod)!
    !end if

    if (epinveins>0) then  ! >>> Is 7-6-14
      a=epinveins ! >>> Is 7-6-14
      a=1.0d0/a        ! >>> Is 7-6-14
      rvx=rsprx+rcilx+(rtorx+rstorx)*a  ! >>> Is 7-6-14 !summing the force components into a net force vector
      rvy=rspry+rcily+(rtory+rstory)*a  ! >>> Is 7-6-14
      rvz=rsprz+rcilz+(rtorz+rstorz)*a  ! >>> Is 7-6-14
    else
      !rvx=rsprx+rcilx  ! >>> Is 7-6-14 !summing the force components into a net force vector
      !rvy=rspry+rcily  ! >>> Is 7-6-14
      !rvz=rsprz+rcilz  ! >>> Is 7-6-14
      rvx=rsprx+rcilx+rtorx  !ECTO MOD
      rvy=rspry+rcily+rtory  !ECTO MOD
      rvz=rsprz+rcilz+rtorz  !ECTO MOD
    end if

    !Physical boundaries force  !>>>>>Miquel9-1-14
    if(d_node(nod)%hold==1)then
      uvx=d_node(nod)%orix-ix
      uvy=d_node(nod)%oriy-iy
      uvz=d_node(nod)%oriz-iz
      d=uvx**2+uvy**2+uvz**2 ; if(d>d_epsilod)then ; ud=1d0/sqrt(d);else; ud=0;end if
      a=d_node(nod)%khold
      rvx=rvx+uvx*a*ud    !now the force is constant        !no need to calculate the unit vector, because here you have the product of the unit vector and the distance,
      rvy=rvy+uvy*a*ud    !to make it a spring remove ud    !wich is the original vector
      rvz=rvz+uvz*a*ud
    end if
    
    if(d_node(nod)%hold==3)then
      uvx=0d0!d_node(nod)%orix-ix
      uvy=0d0!d_node(nod)%oriy-iy
      uvz=d_node(nod)%oriz-iz
      d=uvx**2+uvy**2+uvz**2 ; if(d>d_epsilod)then ; ud=1d0/sqrt(d);else;ud=0;end if
      a=d_node(nod)%khold
      !rvx=rvx+uvx*a*ud    !now the force is constant        !no need to calculate the unit vector, because here you have the product of the unit vector and the distance,
      !rvy=rvy+uvy*a*ud    !to make it a spring remove ud    !wich is the original vector
      rvz=rvz+uvz*a*ud
    end if
    
        !HAIRMOD !!!!!!!!!!!!!!!!!!!!!!!!!!!!>>Miquel10-11-14
    !Buoyancy force on one epithelial side  !>>>>>Miquel13-10-14
    if(d_ffu(6)==1.and. tipi==2)then
      uvx=cx
      uvy=cy
      uvz=cz
      !5d-1!buoyancy constant
      rvx=rvx+uvx*k_bu*udd*(d_gex(nod,2)*d_gen(2)%wa(2))    !now the force is constant        !no need to calculate the unit vector, because here you have the product of the unit vector and the distance,
      rvy=rvy+uvy*k_bu*udd*(d_gex(nod,2)*d_gen(2)%wa(2))    !to make it a spring remove ud    !wich is the original vector
      rvz=rvz+uvz*k_bu*udd*(d_gex(nod,2)*d_gen(2)%wa(2))
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !ii=d_nparam_per_node+16
    !!print*,"ii",ii,"d_nparam_per_node",d_nparam_per_node
    !if(tipi==3 .and. d_npag(ii)>0) then
    !  a=0d0
    !  do k=1,d_npag(ii)
    !    kk=d_whonpag(ii,k)
    !    if (d_gex(nod,kk)>0.0d0) then
    !      a=a+d_gex(nod,kk)*d_gen(kk)%wa(ii)  !wa in units of probability
    !      print*,nod,"rvpre",rvx,rvy,rvz
    !      rvx=rvx+k_bu*cel_pol(nod,1)*a
    !      rvy=rvy+k_bu*cel_pol(nod,2)*a
    !      rvz=rvz+k_bu*cel_pol(nod,3)*a
    !      print*,nod,"polaritzacio a",a,"cel_pol",cel_pol(nod,1),cel_pol(nod,2),cel_pol(nod,3)
    !      print*,nod,"rv",rvx,rvy,rvz
    !
    !    end if
    !  end do
    !end if    
    
    
    pxpypz(nod,3)=sqrt(rvx**2+rvy**2+rvz**2)  !module of the force vector
    pxpypz(nod,1)=rvx ; pxpypz(nod,2)=rvy ; pxpypz(nod,3)=rvz

  !print*,"nod",nod,"d_px",d_px(nod),"d_py",d_py(nod),"d_pz",d_pz(nod),"d_dex",d_dex(nod)
  
987 continue !>>GPU ADD 2-12-14

end subroutine forces_kernel

!!!!!!!!!!!********************************************************************************

!!!!!genestep

subroutine dev_gene_stuff ! gene interactions and gene diffusion
real*8, dimension (:,:) :: digex(nd,ng) ! the increment: only for calculations
integer :: checkn(nd)
integer::istat
real*8,device, dimension (:,:) :: d_digex(nd,ng) ! the increment: only for calculations
integer,device :: d_checkn(nd)
integer::shared_allocation
!real*8,device, dimension (:,:) :: d_dif_dif(nd,ng) ! diffusion differential

  !ldi=epsilod !!!!!!!!!!!!!!!!!!!!!!!!!!! ACHTUNG

  !print*,"debug:: calling genestep_kernel",getot

  !print*,"gex pre",gex(1,1:nd)

  !print*,"w",d_w

  d_digex=0d0 ; d_checkn=0
  
  call genestep_kernel<<<nblocks,nthreads_per_block>>>(d_digex,d_checkn,ng)

  !block_grid=dim3(omnn,1,1)
  
  !shared_allocation=8*ng
  
  !call genestep_diffusion_kernel<<<nd,omnn,shared_allocation>>>(d_dif_dif,d_checkn,ng)
  !istat=cudaDeviceSynchronize
  !call genestep_reaction_kernel<<<blocks,threads_per_block>>>(d_dif_dif,d_digex,ng)

  
!  call genestep_kernel<<<nd,1>>>
  !call cudaDeviceSynchronize
  !istat=cudaDeviceSynchronize

  !print*,"debug:: copying back genetic data",getot

  digex(1:nd,1:ng)=d_digex(1:nd,1:ng)
  !print*,"digex",digex(1,:)
  checkn(1:nd)=d_checkn(1:nd)
  
  !print*,"debug:: updating agex matrix",getot
  
  if (ffu(9)/=2.and.ffu(19)/=2) then
    agex(:nd,1:ng) = gex(:nd,1:ng) + delta*digex(:nd,1:ng) ;        
  else
print *,"RUNGE-KUTTA DOES NOT WORK"
  call exit(12)
!    call rungekutta4gex
  end if

  !do i=1,nd
  !  do j=1,ng
  !   if (abs(agex(i,j))<ldi) gex(i,j)=0
  !  end do
  !end do

  do i=1,nd        
    if (checkn(i)==0) then ! disruption of lateral signalling
      do jj=1,nkindof(7)
        k=wkindof(7,jj)
        if (gen(k)%npre>0) then ! >>> Is 9-5-14 
          kkk=gen(k)%pre(1)   !each notch molecule should have one and only one pre to make anything
          agex(i,kkk)=agex(i,kkk)+gex(i,k)! we lose all the bound forms !!!! RZ 4-3-14
          agex(i,k)=0.0d0     ! we lose all the bound forms
        end if  ! >>> Is 9-5-14
         ! mind that we do not consider partial disruptions here
      end do
    end if
  end do
  !print*,"debug:: dev_gene_stuff done",getot
end subroutine dev_gene_stuff



!*******************************************************************************************

attributes(global) subroutine genestep_kernel(digex,checkn,ng)

integer :: i, k, kk ! RZ 17-11-14 ! added k, kk
integer :: ii1,i1,ii2,i2,ii3,i3,iii1,iii2,iii3
integer :: j,jj,jjj,kkk,iii
real*8  :: a,b
integer :: kki,ivv
real*8  :: sm,smd,smp
integer :: ie,tipi,celi
real*8  :: ix,iy,iz
integer,value::ng
real*8  :: gext(ng)!,agex(nda,ng)
real*8  :: dist,nbh,udist,udistn
integer :: nv,nbo
real*8 :: did(ng)
!integer :: checkn(nd)
!integer :: c_nneigh(nd)
real*8  ::krev,kbound,aa
integer::ii
integer::recep !ECTO MOD

real*8,dimension(:,:) :: digex ! the increment: only for calculations
integer,dimension(:) :: checkn

       !if (allocated(dgex)) deallocate(dgex)
       !allocate(dgex(nd,ng))
       !dgex=0.0d0
       !if (allocated(checkn)) deallocate(checkn)
       !allocate(checkn(nd))
       !checkn=0.0d0
!      rdiffmax=2*maxval(node(:nd)%da)*dmax !now it's in the iteracio subroutine !>>Miquel27-2-14

        !i=blockidx%x  !here each node will be a block
        i=blockDim%x*(blockidx%x-1)+threadidx%x
        if (i>d_nd) goto 457

!print*,"debug: genestep_kernel tall0",i
        did(1:ng)=0.0d0
        gext(1:ng)=0d0  !>>Miquel22-8-14
!print*,"debug: ng",ng,gext(1:ng)

        tipi=d_node(i)%tipus
        celi=d_node(i)%icel !>>Miquel31-3-14
        !if (rdiffmax<2*node(i)%da)   then ;print *,"RRRRRRRRRRRRR"; nbh=2*node(i)%da   ; else ; nbh=rdiffmax ; end if  ! this may make nbh different for different nodes 
        !!! and then non-symmetric diffusion 
        nbh=d_rdiffmax          
        nv=0
        !ix=d_node(i)%x     ; iy=d_node(i)%y     ; iz=d_node(i)%z   
        !ii1=nint(iz*d_urv) ; ii2=nint(iy*d_urv) ; ii3=nint(ix*d_urv)
        ivv=d_node(i)%altre

        nbo=nint(nbh*d_urv)
        !nneigh(i)=0 
        checkn(i)=0 ! RZ 4-3-14

!print*,"debug: genestep_kernel tall1",i,"nneigh",d_nneigh(i)
              !print*,"diffusion HOLA"
        ! ************************* STANDARD MINIMAL LOOP (1,2,3,4) **********************************
        ! DIFFUSION
        
        !if(d_node(i)%border/=0)then
        !  kk=d_node(i)%border
        !  do i1=1,4
        !    ie=d_borders_neigh(kk,i1)
        !    if(ie==0) exit
        !    dist=d_node(i)%req+d_node(ie)%req
        !    !udist=1.0d0/(dist)
        !    udist=1.0d0/(1d0+dist)            
        !    if (d_nkindof(4)>0) then
        !      if (dist<nbh) then
        !        do jjj=1,d_nkindof(4) ! kindof 4 diffuses within the cell but it is only in its surface (so it is not within the cell really)
        !          j=d_wkindof(4,jjj)
        !          if ((d_node(ie)%tipus==1 .and. tipi==2).or.(d_node(ie)%tipus==2 .and. tipi==1)) then 
        !            !did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))/(dist+1) no diffusion across epitelia for signaling molecules
        !            !did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))*udist 
        !          else 
        !            did(j)=did(j)+d_gen(j)%diffu*(d_gex(ie,j)-d_gex(i,j))*udist !;print*,"i,ie",i,ie,"udist",udist
        !          end if
        !        end do
        !      end if
        !    end if
        !  end do
        !end if
            
        do i1=1,d_nneigh(i)
          ie=d_neigh(i,i1)
          dist=d_dneigh(i,i1)
          !udist=1.0d0/(dist) ! >>> Is 25-5-14
          udist=1.0d0/(1d0+dist) ! >>> Is 25-5-14
     !print*,"debug: genestep_kernel tall2",i,i1,ie  
          ! now in case we have kindof 5,6 and 7: NOTICE that these are membrane molecules so they do not diffuse within the cell
          !     then in order to reach the membrane they have to come through a previous form of a lower kindof 
          do jj=1,3
            if (d_nkindof(jj)>0) then
              if (ie==ivv) then   !diffusion with the lower part of the epithelium is always on
                if(tipi<3)then
                  do jjj=1,d_nkindof(jj)
                    j=d_wkindof(jj,jjj) 
                    did(j)=did(j)+d_gen(j)%diffu*(d_gex(ie,j)-d_gex(i,j))*udist
                  end do
                  cycle
                end if
              elseif (celi==d_node(ie)%icel) then
              !print*,"diffusion i",i,"ie",ie,"jj",jj,"dist",dist
                if(dist.lt.nbh) then ! make sure neighbours within given distance are taken into account   
                  do jjj=1,d_nkindof(jj)
                    j=d_wkindof(jj,jjj) !; print*,jj,"j diffu",j,"gex",gex(i,j)
                    did(j)=did(j)+d_gen(j)%diffu*(d_gex(ie,j)-d_gex(i,j))*udist
                    !print*,"j",j,"did",did(j),"diffu",d_gen(j)%diffu,"ie",ie,"d_gex(i)",d_gex(i,j),"d_gex(ie)",d_gex(ie,j),"udist",udist
                    !print*,"nbh",nbh
                  end do
                end if
              end if
            end if
          end do

          !print*,"d_nkindof(4)",d_nkindof(4),"dist",dist,"nbh",nbh
          if (d_nkindof(4)>0) then
            if (dist<nbh) then
              do jjj=1,d_nkindof(4) ! kindof 4 diffuses within the cell but it is only in its surface (so it is not within the cell really)
                j=d_wkindof(4,jjj)
!ECTO MOD       !if ((d_node(ie)%tipus==1 .and. tipi==2).or.(d_node(ie)%tipus==2 .and. tipi==1)) then 
                  !did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))/(dist+1) no diffusion across epitelia for signaling molecules
                  !did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))*udist 
!ECTO MOD       !else 
                  did(j)=did(j)+d_gen(j)%diffu*(d_gex(ie,j)-d_gex(i,j))*udist !;print*,"i,ie",i,ie,"udist",udist
                    !print*,"j",j,"did",did(j),"diffu",d_gen(j)%diffu,"ie",ie,"d_gex(i)",d_gex(i,j),"d_gex(ie)",d_gex(ie,j),"udist",udist
                    !print*,"nbh",nbh
                  !print*,"diffusion????? diffu",d_gen(j)%diffu,"d_gex",d_gex(ie,j),"udist",udist
!ECTO MOD       !end if
              end do
            end if
          end if
          if (d_nkindof(9)>0) then
            if (dist<nbh) then
              do jjj=1,d_nkindof(9) ! kindof 4 diffuses within the cell but it is only in its surface (so it is not within the cell really)
                j=d_wkindof(9,jjj)
!ECTO MOD       !if ((d_node(ie)%tipus==1 .and. tipi==2).or.(d_node(ie)%tipus==2 .and. tipi==1)) then 
                  !did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))/(dist+1) no diffusion across epitelia for signaling molecules
                  !did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))*udist 
!ECTO MOD       !else 
                  did(j)=did(j)+d_gen(j)%diffu*(d_gex(ie,j)-d_gex(i,j))*udist !;print*,"i,ie",i,ie,"udist",udist
                    !print*,"j",j,"did",did(j),"diffu",d_gen(j)%diffu,"ie",ie,"d_gex(i)",d_gex(i,j),"d_gex(ie)",d_gex(ie,j),"udist",udist
                    !print*,"nbh",nbh

!ECTO MOD       !end if
              end do
            end if
          end if
  
          ! kindof 5
          if (d_nkindof(5)>0) then
            if (tipi<3) then
              if (ie==ivv) then   !diffusion with the lower part of the epithelium is always on
                do jjj=1,d_nkindof(5)  !active apical-basal diffusion due to microtubule transport
                  j=d_wkindof(5,jjj)
                  if (tipi==2) then
                    did(j)=did(j)-d_gen(j)%diffu*d_gex(i,j)*udist !loss due to kinesin transport, here diffu is like the transport rate by kinesin !>>>Miquel16-12-13
                  else if(tipi==1)then
                    did(j)=did(j)+d_gen(j)%diffu*d_gex(ie,j)*udist !gain due to kinesin transport, here diffu is like the transport rate by kinesin !>>>Miquel16-12-13
                  end if
                end do
              else
                if(dist.lt.nbh) then ! make sure neighbours within given distance are taken into account   
                  if(tipi==d_node(ie)%tipus)then
                    do jjj=1,d_nkindof(5)
                      j=d_wkindof(5,jjj)
                      did(j)=did(j)+d_gen(j)%diffu*(d_gex(ie,j)-d_gex(i,j))*udist
                    end do
                  end if
                end if
              end if
            end if
          end if
          if (d_nkindof(6)>0) then
            if (tipi<3) then
              if (ie==ivv) then   !diffusion with the lower part of the epithelium is always on
                do jjj=1,d_nkindof(6)  !active apical-basal diffusion due to microtubule transport
                  j=d_wkindof(6,jjj)
                  if(tipi==2)then
                    did(j)=did(j)-d_gen(j)%diffu*d_gex(i,j)*udist !loss due to kinesin transport, here diffu is like the transport rate by kinesin !>>>Miquel16-12-13
                  else if(tipi==1)then
                    did(j)=did(j)+d_gen(j)%diffu*d_gex(ie,j)*udist !gain due to kinesin transport, here diffu is like the transport rate by kinesin !>>>Miquel16-12-13
                  end if
                end do
              else
                if(dist.lt.nbh) then ! make sure neighbours within given distance are taken into account  
                  if(tipi==d_node(ie)%tipus)then
                    do jjj=1,d_nkindof(6)
                      j=d_wkindof(6,jjj)
                      did(j)=did(j)+d_gen(j)%diffu*(d_gex(ie,j)-d_gex(i,j))*udist
                    end do
                  end if
                end if
              end if
            end if
          end if
          if (d_nkindof(7)>0) then !membrane-bound notch-delta kind of interaction >>> Is 9-2-14 
            if(tipi<4) then ! there is no reason to exclude mesenchymal cells RZ 4-3-14
              if (d_node(ie)%icel.ne.celi) then   ! ONLY from different cells
                if(tipi==d_node(ie)%tipus)then !>>Miquel2-10-14
                  if(dist.le.(d_node(i)%da+d_node(ie)%da)) then     !new kinetics !>>Miquel22-8-14
                    udistn=(d_node(i)%da+d_node(ie)%da)/(d_node(i)%da+d_node(ie)%da+dist) !notice this is the inverse Is
                    checkn(i)=1 !nodes close enough for binding
                    do jj=1,d_nkindof(7)
                      j=d_wkindof(7,jj)
                      if (d_gen(j)%npost>0) then !the dissociation reaction
                        do kki=1,d_gen(j)%npost   
                          !smd=0.0d0
                          kkk=d_gen(j)%post(kki) !the post form (either ligand or receptor)
                          do k=1,d_nwpost(j,kki)
                            kk=d_gen(int(d_wpost(j,kki,k,1)))%post(1) !the other bound form  ! RZ 17-11-14 INT
                            krev= d_wpost(j,kki,k,2) !dissociation constant
                          end do                           
                          b=krev*d_gex(i,j)*d_gex(ie,kk)   !loss of bound form by dissociation
                          gext(kkk) = gext(kkk) + b !this is the gain of ligand or receptor by dissociation of bound form
                          gext(j)   = gext(j)   - b !this is the loss of bound form by dissociation (it loses 1 for every couple of ligand-receptor released)
                        end do
                      end if
                      if (d_gen(j)%npre>0) then !the binding reaction
                        do kki=1,d_gen(j)%npre   
                          kkk=d_gen(j)%pre(kki) !the pre form (either ligand or receptor)
                          do k=1,d_nwpre(j,kki)
                            kk=d_wpre(j,kki,k,1) !catalist (either ligand or receptor) and also the one
                            kbound=d_wpre(j,kki,k,2) !binding constant
                            smd=d_gex(i,kkk)*d_gex(ie,kk) !productory of the concentrations of ligand and receptor (kk and kkk)
                          end do
                          b=kbound*smd         ! gain of bound form by binding of receptor and ligand
                          gext(kkk) = gext(kkk) - b !this is the loss of dissociated form due to binding
                          gext(j)   = gext(j)   + b !this is the gain of bound form
                        end do
                      end if
                    end do
                  end if
                end if
              else
                if(dist.lt.nbh) then ! make sure neighbours within given distance are taken into account  
                  if(tipi==d_node(ie)%tipus)then
                    do jjj=1,d_nkindof(7)
                      j=d_wkindof(7,jjj)
                      did(j)=did(j)+d_gen(j)%diffu*(d_gex(ie,j)-d_gex(i,j))*udist
                    end do
                  end if
                end if
              end if
            end if
          end if
          if (d_nkindof(8)>0) then  ! receptors of extracellular diffusible molecules ONLY THE ACTIVE FORM, the inactive is kindof2 or 3
            if (celi==d_node(ie)%icel) then
              if (tipi==d_node(ie)%tipus) then ! the active form diffuses only in the membrane (the inactive everywhere)
                if(dist.lt.nbh) then ! make sure neighbours within given distance are taken into account   
                  do jjj=1,d_nkindof(8)
                    j=d_wkindof(8,jjj)
                    did(j)=did(j)+d_gen(j)%diffu*(d_gex(ie,j)-d_gex(i,j))*udist
                  end do
                end if
              end if
            end if
          end if
        end do

        !loss by diffusion when there are borders
        if(d_node(i)%border==1 .and. (d_nkindof(4)>0 .or. d_nkindof(9)>0))then
          do jjj=1,d_nkindof(4) ! kindof 4 diffuses within the cell but it is only in its surface (so it is not within the cell really)
            j=d_wkindof(4,jjj)
            did(j)=did(j)-d_gen(j)%diffu*d_gex(i,j)
          end do
          do jjj=1,d_nkindof(9) ! kindof 4 diffuses within the cell but it is only in its surface (so it is not within the cell really)
            j=d_wkindof(9,jjj)
            did(j)=did(j)-d_gen(j)%diffu*d_gex(i,j)
          end do
        end if

!print*,"debug: genestep_kernel tall3",i  


        !LIMIT TO DIFFUSION, so that diffusion does not get to infinite distance with small amounts
        !do j=1,ng
        !  if (abs(did(j))<ldi) did(j)=0.0d0
        !end do

        !REACTION 
        if (d_node(i)%marge==0.and.d_node(i)%tipus<4) then ! IN THE NUCLEUS (THERE WILL BE TRANSCRIPTION)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NEWEST VERSIONS
          do jjj=1,d_nkindof(1)           !*****************PRODUCTION OF PRIMARY FORMS (WITHOUT FURTHER FORMS)
            j=d_wkindof(1,jjj)     
            a=d_gex(i,j)
            sm=0.0d0
            do k=1,int(d_nw(j))
              kk=d_w(j,k)
              !crappy implementation of signal receptors   !ECTO MOD
              if((d_gen(kk)%kindof==4 .or. d_gen(kk)%kindof==9).and. d_gen(kk)%npre>0)then
                recep=0
                do kki=1,d_gen(kk)%npre   
                  kkk=d_gen(kk)%pre(kki)
                  if(d_gex(i,kkk)>d_epsilod) recep=recep+1
                end do
                if(recep==0) cycle
              end if
              !crappy implementation of signal receptors              
              sm = sm + d_gen(j)%w(kk)*d_gex(i,kk)  !the amount of TFs
            end do
            if (sm<0.0) sm=0.0
            gext(j) = sm/(1+sm) - d_gen(j)%mu*a + did(j) !gain of product by transcription (Michaelis-Menten) !>>>> Miquel 20-12-13
          end do
          do jjj=1,d_nkindof(4)           !*****************PRODUCTION OF PRIMARY FORMS (WITHOUT FURTHER FORMS)
            j=d_wkindof(4,jjj)     
            a=d_gex(i,j)
            sm=0.0d0
            do k=1,int(d_nw(j))
              kk=d_w(j,k)
              !crappy implementation of signal receptors   !ECTO MOD
              if((d_gen(kk)%kindof==4 .or. d_gen(kk)%kindof==9).and. d_gen(kk)%npre>0)then
                recep=0
                do kki=1,d_gen(kk)%npre   
                  kkk=d_gen(kk)%pre(kki)
                  if(d_gex(i,kkk)>d_epsilod) recep=recep+1
                end do
                if(recep==0) cycle
              end if
              !crappy implementation of signal receptors              
              sm = sm + d_gen(j)%w(kk)*d_gex(i,kk)  !the amount of TFs
            end do
            if (sm<0.0) sm=0.0
            gext(j) = sm/(1+sm) !- d_gen(j)%mu*a + did(j) !gain of product by transcription (Michaelis-Menten) !>>>> Miquel 20-12-13
          end do
          do jjj=1,d_nkindof(9)           !*****************PRODUCTION OF PRIMARY FORMS (WITHOUT FURTHER FORMS)
            j=d_wkindof(9,jjj)
            a=d_gex(i,j)
            sm=0.0d0 ; smn=0.0d0
            do k=1,int(d_nw(j))
              kk=d_w(j,k)
              !crappy implementation of signal receptors
              if(d_gen(kk)%kindof==9 .and. d_gen(kk)%npre>0)then
                recep=0
                do kki=1,d_gen(kk)%npre   
                  kkk=d_gen(kk)%pre(kki)
                  if(d_gex(i,kkk)>d_epsilod) recep=recep+1
                end do
                if(recep==0) cycle
              end if
              !crappy implementation of signal receptors
              !sm=sm+d_gen(j)%w(kk)*d_gex(i,kk)
              aa=d_gen(j)%w(kk)
              if(aa>=d_epsilod)then
                sm = sm + aa*(d_gex(i,kk)**2)  !the amount of TFs
              !  sm = sm + aa*d_gex(i,kk)  !the amount of TFs
              else
                smn = smn + abs(aa)*d_gex(i,kk)  !the amount of TFs
              !  smn = smn + abs(aa)*d_gex(i,kk)  !the amount of TFs
              end if
            end do
            if (sm<-d_epsilod) sm=0.0
            gext(j) = gext(j) + sm/(1+smn) !- d_gen(j)%mu*a + did(j) !gain of product by transcription (Michaelis-Menten) !>>>> Miquel 20-1
            !gext(j) = gext(j) + sm - d_gen(j)%mu*a**2 + did(j) !gain of product by transcription (Michaelis-Menten) !>>>> Miquel 20-1
          end do
          do jjj=1,d_nkindof(2)           !*****************PRODUCTION OF PRIMARY FORMS (WITH FURTHER FORMS)
            j=d_wkindof(2,jjj)
            a=d_gex(i,j)
            sm=0.0d0
            do k=1,int(d_nw(j))
              kk=d_w(j,k)
              !crappy implementation of signal receptors
              if((d_gen(kk)%kindof==4 .or. d_gen(kk)%kindof==9).and. d_gen(kk)%npre>0)then
                recep=0
                do kki=1,d_gen(kk)%npre   
                  kkk=d_gen(kk)%pre(kki)
                  if(d_gex(i,kkk)>d_epsilod) recep=recep+1
                end do
                if(recep==0) cycle
              end if
              !crappy implementation of signal receptors              
              sm = sm + d_gen(j)%w(kk)*d_gex(i,kk)  !the amount of TFs
            end do
            if (sm<0.0) sm=0.0
            if (d_gen(j)%npost>0) then 
              do kki=1,d_gen(j)%npost   
                kkk=d_gen(j)%post(kki); if(d_gen(kkk)%kindof==8) cycle !>>Miquel21-8-14
                smd=0.0d0
                do k=1,d_nwpost(j,kki)
                  kk=int(d_wpost(j,kki,k,1)) ! RZ 17-11-14 INT
                  smd=smd + d_wpost(j,kki,k,2)*d_gex(i,kk)   !loss of product by transforming into another form !>>> Miquel20-12-13
                end do
                if (smd<0.0) smd=0.0   ! Is 9-2-14
                b=smd*a/(1+a)         ! Michaelis-Menten
                gext(kkk) = gext(kkk) + b ! that is the gain by the post due to its transformation from j
                gext(j)   = gext(j)   - b ! that is the lost j has due to its transformation to j
              end do
              gext(j) = gext(j) + sm/(1+sm) - d_gen(j)%mu*a + did(j) !>>> Is 9-2-14
            else
              gext(j) = sm/(1+sm) - d_gen(j)%mu*a + did(j)  !>>> Miquel20-12-13
            end if
          end do
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
          !do jjj=1,d_nkindof(9)           !*****************PRODUCTION OF PRIMARY FORMS (WITHOUT FURTHER FORMS)
          !  j=d_wkindof(9,jjj)     
          !  a=d_gex(i,j)
          !  sm=0.0d0 ; smn=0.0d0
          !  !print*,i,"j",j,"a",a
          !  do k=1,int(d_nw(j))
          !    kk=d_w(j,k)
          !    aa=d_gen(j)%w(kk)
          !    !print*,i,"j",j,"kk",kk,"npost",d_gen(kk)%npost
          !    if(d_gen(kk)%kindof==4)then
          !      if(d_gen(kk)%npre>0)then
          !        do kki=1,d_gen(kk)%npre  
          !          kkk=d_gen(kk)%pre(kki)
          !          !print*,"entra"
          !          !print*,i,"j",j,"kk",kk!,"post1",d_gen(kk)%post(1)
          !          !print*,i,"j",j,"kk",kk,"kkk",kkk
          !          !print*,i,"j",j,"kk",kk,"dgex(i,kkk)",d_gex(i,kkk)
          !          if(d_gex(i,kkk)>d_epsilod)then
          !            if(aa>=d_epsilod)then
          !              sm = sm + aa*(d_gex(i,kk)**2)  !the amount of TFs
          !            else
          !              smn = smn + abs(aa)*d_gex(i,kk)  !the amount of TFs
          !            end if
          !            exit
          !          end if
          !        end do
          !      end if
          !    else
          !      if(aa>=d_epsilod)then
          !        sm = sm + aa*(d_gex(i,kk)**2)  !the amount of TFs
          !      else
          !        smn = smn + abs(aa)*d_gex(i,kk)  !the amount of TFs
          !      end if
          !    endif
          !  end do
          !  !if (sm<0.0) sm=0.0
          !  gext(j) = gext(j) + sm/(1+smn) - d_gen(j)%mu*a + did(j) !gain of product by transcription (Michaelis-Menten) !>>>> Miquel 20-12-13
          !end do

          !!!!!!!!!!!!!!!!!!!!!!! OLD VERSION
          !do jjj=1,d_nkindof(1)           !*****************PRODUCTION OF PRIMARY FORMS (WITHOUT FURTHER FORMS)
          !  j=d_wkindof(1,jjj)     
          !  a=d_gex(i,j)
          !  sm=0.0d0
          !  do k=1,int(d_nw(j))
          !    kk=d_w(j,k)
          !    sm = sm + d_gen(j)%w(kk)*d_gex(i,kk)  !the amount of TFs
          !  end do
          !  if (sm<0.0) sm=0.0
          !  gext(j) = sm/(1+sm) - d_gen(j)%mu*a + did(j) !gain of product by transcription (Michaelis-Menten) !>>>> Miquel 20-12-13
          !end do
          !do jjj=1,d_nkindof(4)           !*****************PRODUCTION OF PRIMARY FORMS (WITHOUT FURTHER FORMS)
          !  j=d_wkindof(4,jjj)
          !  a=d_gex(i,j)
          !  sm=0.0d0 ; smn=0.0d0
          !  do k=1,int(d_nw(j))
          !    kk=d_w(j,k)
          !    aa=d_gen(j)%w(kk)
          !    if(aa>=d_epsilod)then
          !      sm = sm + aa*(d_gex(i,kk)**2)  !the amount of TFs
          !    else
          !      smn = smn + abs(aa)*d_gex(i,kk)  !the amount of TFs
          !    end if
          !  end do
          !  !if (sm<0.0) sm=0.0
          !  gext(j) = gext(j) + sm/(1+smn) - d_gen(j)%mu*a + did(j) !gain of product by transcription (Michaelis-Menten) !>>>> Miquel 20-1
          !end do           
          !do jjj=1,d_nkindof(2)           !*****************PRODUCTION OF PRIMARY FORMS (WITH FURTHER FORMS)
          !  j=d_wkindof(2,jjj)
          !  a=d_gex(i,j)
          !  sm=0.0d0
          !  do k=1,int(d_nw(j))
          !    kk=d_w(j,k)
          !    sm = sm + d_gen(j)%w(kk)*d_gex(i,kk)  !the amount of TFs
          !  end do
          !  if (sm<0.0) sm=0.0
          !  if (d_gen(j)%npost>0) then 
          !    do kki=1,d_gen(j)%npost   
          !      kkk=d_gen(j)%post(kki); if(d_gen(kkk)%kindof==8) cycle !>>Miquel21-8-14
          !      smd=0.0d0
          !      do k=1,d_nwpost(j,kki)
          !        kk=int(d_wpost(j,kki,k,1)) ! RZ 17-11-14 INT
          !        smd=smd + d_wpost(j,kki,k,2)*d_gex(i,kk)   !loss of product by transforming into another form !>>> Miquel20-12-13
          !      end do
          !      if (smd<0.0) smd=0.0   ! Is 9-2-14
          !      b=smd*a/(1+a)         ! Michaelis-Menten
          !      gext(kkk) = gext(kkk) + b ! that is the gain by the post due to its transformation from j
          !      gext(j)   = gext(j)   - b ! that is the lost j has due to its transformation to j
          !    end do
          !    gext(j) = gext(j) + sm/(1+sm) - d_gen(j)%mu*a + did(j) !>>> Is 9-2-14
          !  else
          !    gext(j) = sm/(1+sm) - d_gen(j)%mu*a + did(j)  !>>> Miquel20-12-13
          !  end if
          !end do
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
        else  !****************** NOT IN THE NUCLEUS (NO TRANSCRIPTION)

          do jjj=1,d_nkindof(1)           !*****************PRODUCTION OF PRIMARY FORMS (WITHOUT FURTHER FORMS)
            j=d_wkindof(1,jjj)  
            gext(j) = - d_gen(j)%mu*d_gex(i,j) + did(j) !gain of product by transcription (Michaelis-Menten) !>>>> Miquel 20-12-13
          end do 
          do jjj=1,d_nkindof(2)           !*****************PRODUCTION OF PRIMARY FORMS (WITH FURTHER FORMS)
            j=d_wkindof(2,jjj)     
            a=d_gex(i,j)
            if (d_gen(j)%npost>0) then !
              do kki=1,d_gen(j)%npost
                kkk=d_gen(j)%post(kki); if(d_gen(kkk)%kindof==8) cycle !>>Miquel21-8-14
                smd=0.0d0
                do k=1,d_nwpost(j,kki)
                  kk=int(d_wpost(j,kki,k,1)) ! RZ 17-11-14 INT
                  smd=smd + d_wpost(j,kki,k,2)*d_gex(i,kk)   !loss of product by transforming into another form !>>> Miquel20-12-13
                end do                           
                if (smd<0.0) smd=0.0   ! Is 9-2-14
                b=smd*a/(1+a)         ! Michaelis-Menten
                gext(kkk) = gext(kkk) + b ! that is the gain by the post due to its transformation from j
                gext(j)   = gext(j)   - b ! that is the lost j has due to its transformation to j
              end do
              gext(j) = gext(j) - d_gen(j)%mu*a + did(j) !>>> Is 9-2-14
            else
              gext(j) = - d_gen(j)%mu*a + did(j)  !>>> Miquel20-12-13
            end if
          end do
        end if
        
        !********PRODUCTION OF SECONDARY FORMS: they are never transcripts : so it is the same whether we are in the nucleus or not
        do jjj=1,d_nkindof(3)           ! THIS IS IMPORTANT, IF YOU WANT TO MAKE A NOTCH YOU SHOULD MAKE ITS PREVIOUS FORM AND TRANSCRIBE IT
          j=d_wkindof(3,jjj)
          a=d_gex(i,j)
          if (d_gen(j)%npost>0) then 
            do kki=1,d_gen(j)%npost
              kkk=d_gen(j)%post(kki); if(d_gen(kkk)%kindof==8) cycle !>>Miquel21-8-14
              smd=0.0d0
              do k=1,d_nwpost(j,kki)
                kk=int(d_wpost(j,kki,k,1))  ! RZ 17-11-14 INT
                smd=smd + d_wpost(j,kki,k,2)*d_gex(i,kk)   !loss of product by transforming into another form !>>> Miquel20-12-13
              end do                           
              if (smd<0.0) smd=0.0   ! Is 9-2-14
              b=smd*a/(1+a)         ! Michaelis-Menten
              gext(kkk) = gext(kkk) + b ! that is the gain by the post due to its transformation from j
              gext(j)   = gext(j)   - b ! that is the lost j has due to its transformation to j
            end do
            !gext(j) = gext(j) - gen(j)%mu*a + did(j) !>>> Is 9-2-14
          end if
          gext(j) = gext(j)  - d_gen(j)%mu*a + did(j)  ! Is 9-2-14
        end do 

        !********PRODUCTION OF SECONDARY FORMS: they are never transcripts : so it is the same whether we are in the nucleus or not
        do jjj=1,d_nkindof(4)           ! THIS IS IMPORTANT, IF YOU WANT TO MAKE A NOTCH YOU SHOULD MAKE ITS PREVIOUS FORM AND TRANSCRIBE IT
          j=d_wkindof(4,jjj)     
          a=d_gex(i,j)
          if (d_gen(j)%npost>0) then 
            do kki=1,d_gen(j)%npost   
              kkk=d_gen(j)%post(kki); if(d_gen(kkk)%kindof==8) cycle !>>Miquel21-8-14
              smd=0.0d0
              do k=1,d_nwpost(j,kki)
                kk=int(d_wpost(j,kki,k,1)) ! RZ 17-11-14 INT
                smd=smd + d_wpost(j,kki,k,2)*d_gex(i,kk)   !loss of product by transforming into another form !>>> Miquel20-12-13
              end do                           
              if (smd<0.0) smd=0.0   ! Is 9-2-14
              b=smd*a/(1+a)         ! Michaelis-Menten
              gext(kkk) = gext(kkk) + b ! that is the gain by the post due to its transformation from j
              gext(j)   = gext(j)   - b ! that is the lost j has due to its transformation to j
            end do
            !gext(j) = gext(j) - gen(j)%mu*a + did(j) !>>> Is 9-2-14
          end if

          ! NOTICE: growth factors do not return any thing to pre forms, there is not back reaction, because they get secreted as they get produced
          gext(j) = gext(j)  - d_gen(j)%mu*a + did(j)  ! Is 9-2-14
        end do
        
        do jjj=1,d_nkindof(9)           ! THIS IS IMPORTANT, IF YOU WANT TO MAKE A NOTCH YOU SHOULD MAKE ITS PREVIOUS FORM AND TRANSCRIBE IT
          j=d_wkindof(9,jjj)     
          a=d_gex(i,j)
          if (d_gen(j)%npost>0) then 
            do kki=1,d_gen(j)%npost   
              kkk=d_gen(j)%post(kki); if(d_gen(kkk)%kindof==8) cycle !>>Miquel21-8-14
              smd=0.0d0
              do k=1,d_nwpost(j,kki)
                kk=int(d_wpost(j,kki,k,1)) ! RZ 17-11-14 INT
                smd=smd + d_wpost(j,kki,k,2)*d_gex(i,kk)   !loss of product by transforming into another form !>>> Miquel20-12-13
              end do                           
              if (smd<0.0) smd=0.0   ! Is 9-2-14
              b=smd*a/(1+a)         ! Michaelis-Menten
              gext(kkk) = gext(kkk) + b ! that is the gain by the post due to its transformation from j
              gext(j)   = gext(j)   - b ! that is the lost j has due to its transformation to j
            end do
            !gext(j) = gext(j) - gen(j)%mu*a + did(j) !>>> Is 9-2-14
          end if
        
          ! NOTICE: growth factors do not return any thing to pre forms, there is not back reaction, because they get secreted as they get produced
          gext(j) = gext(j)  - d_gen(j)%mu*a + did(j)  ! Is 9-2-14
        end do
        
        do jj=5,7 !********PRODUCTION OF OTHER SECONDARY FORMS: they are never transcripts : so it is the same whether we are in the nucleus or not
        !  if (jj==7) cycle   !>>> Is 1-3-14 ! outcommented RZ 4-3-14
          if (jj==7) then ! >>> RZ 4-3-14
            do jjj=1,d_nkindof(jj)
              j=d_wkindof(jj,jjj)     
              a=d_gex(i,j)
              gext(j) = gext(j)  - d_gen(j)%mu*a + did(j)  ! Is 9-2-14
            end do
          else   ! <<< RZ 4-3-14
            do jjj=1,d_nkindof(jj)           ! THIS IS IMPORTANT, IF YOU WANT TO MAKE A NOTCH YOU SHOULD MAKE ITS PREVIOUS FORM AND TRANSCRIBE IT
              j=d_wkindof(jj,jjj)     
              a=d_gex(i,j)
              if (d_gen(j)%npost>0) then 
                do kki=1,d_gen(j)%npost   
                  kkk=d_gen(j)%post(kki); if(d_gen(kkk)%kindof==8) cycle !>>Miquel21-8-14
                  smd=0.0d0
                  do k=1,d_nwpost(j,kki)
                    kk=int(d_wpost(j,kki,k,1)) ! RZ 17-11-14 INT
                    smd=smd + d_wpost(j,kki,k,2)*d_gex(i,kk)   !loss of product by transforming into another form !>>> Miquel20-12-13
                  end do                           
                  if (smd<0.0) smd=0.0   ! Is 9-2-14
                  b=smd*a/(1+a)         ! Michaelis-Menten
                  gext(kkk) = gext(kkk) + b ! that is the gain by the post due to its transformation from j
                  gext(j)   = gext(j)   - b ! that is the lost j has due to its transformation to j
                end do
                !gext(j) = gext(j) - gen(j)%mu*a + did(j) !>>> Is 9-2-14
              end if
              gext(j) = gext(j)  - d_gen(j)%mu*a + did(j)  ! Is 9-2-14
            end do 
          end if  ! RZ 4-3-14
        end do

        !bound receptors !>>Miquel18-8-14
        do jjj=1,d_nkindof(8)
          j=d_wkindof(8,jjj)
          if (d_gen(j)%npost<1) cycle !they should though ! Is 1-10-14
          a=d_gex(i,j)
          !the dissociation reaction
          krev=0d0
          kkk=d_gen(j)%post(1) !the post form (either ligand or receptor)
          iii=d_gen(j)%post(2) !the post form (either ligand or receptor)
          krev= d_wpost(j,1,1,2) !dissociation constant
          b=krev*a       !loss of bound form by dissociation
          gext(kkk) = gext(kkk) + b !this is the gain of ligand or receptor by dissociation of bound form
          gext(iii) = gext(iii) + b !this is the gain of ligand or receptor by dissociation of bound form

          gext(j)   = gext(j)   - b !this is the loss of bound form by dissociation (it loses 1 for every couple of ligand-receptor released)

          !the binding reaction
          if (d_gen(j)%npre<1) cycle !they should though !>>Miquel1-10-14
          kbound=0d0
          kkk=d_gen(j)%pre(1) !the pre form (either ligand or receptor)
          iii=d_gen(j)%pre(2) !the pre form (either ligand or receptor)
          kbound=d_wpre(j,1,1,2) !binding constant
          smd=d_gex(i,iii)*d_gex(i,kkk) !productory of the concentrations of ligand and receptor
          b=kbound*smd         ! gain of bound form by binding of receptor and ligand

          gext(kkk) = gext(kkk) - b !this is the loss of dissociated form due to binding
          gext(iii) = gext(iii) - b !this is the loss of dissociated form due to binding

          gext(j)   = gext(j)   + b !this is the gain of bound form due to binding
          gext(j) = gext(j)  - d_gen(j)%mu*a + did(j)  ! Is 9-2-14
        end do

!print*,"debug: genestep_kernel tall4",i 
        !do j=1,ng
        !print*,"debug: genestep_kernel tall4.",j,gext(j)
        !  d_dgex(i,j)=gext(j)
        !end do
        digex(i,1:ng)=gext(1:ng)  ! this is the increment in the gene this step
!print*,"debug: genestep_kernel tall5",i 
457 continue
end subroutine genestep_kernel

!*******************************************************************************************

!attributes(global) subroutine genestep_diffusion_kernel(dif_dif,checkn,ng)
!
!integer :: i, k, kk ! RZ 17-11-14 ! added k, kk
!!integer :: ii1,i1,ii2,i2,ii3,i3,iii1,iii2,iii3
!integer :: j,jj,jjj,kkk,iii
!real*8  :: a,b
!integer :: kki,ivv
!real*8  :: sm,smd,smp
!integer :: ie,tipi,celi
!real*8  :: ix,iy,iz
!integer,value::ng
!real*8  :: gext(ng)!,agex(nda,ng)
!real*8  :: dist,nbh,udist,udistn
!integer :: nbo
!!real*8 :: did(ng)
!!integer :: checkn(nd)
!!integer :: c_nneigh(nd)
!real*8  ::krev,kbound
!
!!real*8,dimension(:,:) :: digex ! the increment: only for calculations
!real*8,dimension(:,:) :: dif_dif ! the increment: only for calculations
!integer,dimension(:) :: checkn
!integer::nod
!integer::lnneigh
!real*8,shared :: did(ng)
!
!       !if (allocated(dgex)) deallocate(dgex)
!       !allocate(dgex(nd,ng))
!       !dgex=0.0d0
!       !if (allocated(checkn)) deallocate(checkn)
!       !allocate(checkn(nd))
!       !checkn=0.0d0
!!      rdiffmax=2*maxval(node(:nd)%da)*dmax !now it's in the iteracio subroutine !>>Miquel27-2-14
!
!        !i=blockidx%x  !here each node will be a block
!        nod=blockIdx%x    !each block is a node
!        i=threadidx%x
!        
!        lnneigh=d_nneigh(nod)
!        
!        if (i>lnneigh) goto 457
!
!!print*,"debug: genestep_kernel tall0",i
!        did(1:ng)=0.0d0
!        !gext(1:ng)=0d0  !>>Miquel22-8-14
!!print*,"debug: ng",ng,gext(1:ng)
!
!        tipi=d_node(i)%tipus
!        celi=d_node(i)%icel !>>Miquel31-3-14
!        !if (rdiffmax<2*node(i)%da)   then ;print *,"RRRRRRRRRRRRR"; nbh=2*node(i)%da   ; else ; nbh=rdiffmax ; end if  ! this may make nbh different for different nodes 
!        !!! and then non-symmetric diffusion 
!        nbh=d_rdiffmax
!        !nv=0
!        ix=d_node(i)%x     ; iy=d_node(i)%y     ; iz=d_node(i)%z   
!        !ii1=nint(iz*d_urv) ; ii2=nint(iy*d_urv) ; ii3=nint(ix*d_urv)
!        ivv=d_node(i)%altre
!
!        !nbo=nint(nbh*d_urv)
!        !nneigh(i)=0 
!        if(i==1) checkn(nod)=0 ! RZ 4-3-14
!
!!print*,"debug: genestep_kernel tall1",i,"nneigh",d_nneigh(i)
!
!        ! ************************* STANDARD MINIMAL LOOP (1,2,3,4) **********************************
!        ! DIFFUSION
!        !do i1=1,d_nneigh(i)
!          ie=d_neigh(nod,i)
!          dist=d_dneigh(nod,i)
!          udist=1.0d0/(1d0+dist) ! >>> Is 25-5-14
!!print*,"debug: genestep_kernel tall2",i,i1,ie  
!          ! now in case we have kindof 5,6 and 7: NOTICE that these are membrane molecules so they do not diffuse within the cell
!          !     then in order to reach the membrane they have to come through a previous form of a lower kindof 
!          do jj=1,3
!            if (d_nkindof(jj)>0) then
!              if (ie==ivv) then   !diffusion with the lower part of the epithelium is always on
!                if(tipi<3)then
!                  do jjj=1,d_nkindof(jj)
!                    j=d_wkindof(jj,jjj) 
!                    did(j)=did(j)+d_gen(j)%diffu*(d_gex(ie,j)-d_gex(i,j))*udist
!                  end do
!                  cycle
!                end if
!              elseif (celi==d_node(ie)%icel) then
!                if(dist.lt.nbh) then ! make sure neighbours within given distance are taken into account   
!                  do jjj=1,d_nkindof(jj)
!                    j=d_wkindof(jj,jjj) !; print*,jj,"j diffu",j,"gex",gex(i,j)
!                    did(j)=did(j)+d_gen(j)%diffu*(d_gex(ie,j)-d_gex(i,j))*udist
!                  end do
!                end if
!              end if
!            end if
!          end do
!
!          !print*,"d_nkindof(4)",d_nkindof(4),"dist",dist,"nbh",nbh
!          if (d_nkindof(4)>0) then
!            if (dist<nbh) then
!              do jjj=1,d_nkindof(4) ! kindof 4 diffuses within the cell but it is only in its surface (so it is not within the cell really)
!                j=d_wkindof(4,jjj)
!                if ((d_node(ie)%tipus==1 .and. tipi==2).or.(d_node(ie)%tipus==2 .and. tipi==1)) then 
!                  !did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))/(dist+1) no diffusion across epitelia for signaling molecules
!                  !did(j)=did(j)+gen(j)%diffu*(gex(ie,j)-gex(i,j))*udist 
!                else 
!                  did(j)=did(j)+d_gen(j)%diffu*(d_gex(ie,j)-d_gex(i,j))*udist !;print*,"i,ie",i,ie,"udist",udist
!                  !print*,"diffusion????? diffu",d_gen(j)%diffu,"d_gex",d_gex(ie,j),"udist",udist
!                end if
!              end do
!            end if
!          end if
!  
!          ! kindof 5
!          if (d_nkindof(5)>0) then
!            if (tipi<3) then
!              if (ie==ivv) then   !diffusion with the lower part of the epithelium is always on
!                do jjj=1,d_nkindof(5)  !active apical-basal diffusion due to microtubule transport
!                  j=d_wkindof(5,jjj)
!                  if (tipi==2) then
!                    did(j)=did(j)-d_gen(j)%diffu*d_gex(i,j)*udist !loss due to kinesin transport, here diffu is like the transport rate by kinesin !>>>Miquel16-12-13
!                  else if(tipi==1)then
!                    did(j)=did(j)+d_gen(j)%diffu*d_gex(ie,j)*udist !gain due to kinesin transport, here diffu is like the transport rate by kinesin !>>>Miquel16-12-13
!                  end if
!                end do
!              else
!                if(dist.lt.nbh) then ! make sure neighbours within given distance are taken into account   
!                  if(tipi==d_node(ie)%tipus)then
!                    do jjj=1,d_nkindof(5)
!                      j=d_wkindof(5,jjj)
!                      did(j)=did(j)+d_gen(j)%diffu*(d_gex(ie,j)-d_gex(i,j))*udist
!                    end do
!                  end if
!                end if
!              end if
!            end if
!          end if
!          if (d_nkindof(6)>0) then
!            if (tipi<3) then
!              if (ie==ivv) then   !diffusion with the lower part of the epithelium is always on
!                do jjj=1,d_nkindof(6)  !active apical-basal diffusion due to microtubule transport
!                  j=d_wkindof(6,jjj)
!                  if(tipi==2)then
!                    did(j)=did(j)-d_gen(j)%diffu*d_gex(i,j)*udist !loss due to kinesin transport, here diffu is like the transport rate by kinesin !>>>Miquel16-12-13
!                  else if(tipi==1)then
!                    did(j)=did(j)+d_gen(j)%diffu*d_gex(ie,j)*udist !gain due to kinesin transport, here diffu is like the transport rate by kinesin !>>>Miquel16-12-13
!                  end if
!                end do
!              else
!                if(dist.lt.nbh) then ! make sure neighbours within given distance are taken into account  
!                  if(tipi==d_node(ie)%tipus)then
!                    do jjj=1,d_nkindof(6)
!                      j=d_wkindof(6,jjj)
!                      did(j)=did(j)+d_gen(j)%diffu*(d_gex(ie,j)-d_gex(i,j))*udist
!                    end do
!                  end if
!                end if
!              end if
!            end if
!          end if
!          if (d_nkindof(7)>0) then !membrane-bound notch-delta kind of interaction >>> Is 9-2-14 
!            if(tipi<4) then ! there is no reason to exclude mesenchymal cells RZ 4-3-14
!              if (d_node(ie)%icel.ne.celi) then   ! ONLY from different cells
!                if(tipi==d_node(ie)%tipus)then !>>Miquel2-10-14
!                  if(dist.le.(d_node(i)%da+d_node(ie)%da)) then     !new kinetics !>>Miquel22-8-14
!                    udistn=(d_node(i)%da+d_node(ie)%da)/(d_node(i)%da+d_node(ie)%da+dist) !notice this is the inverse Is
!                    checkn(nod)=1 !nodes close enough for binding
!                    do jj=1,d_nkindof(7)
!                      j=d_wkindof(7,jj)
!                      if (d_gen(j)%npost>0) then !the dissociation reaction
!                        do kki=1,d_gen(j)%npost   
!                          !smd=0.0d0
!                          kkk=d_gen(j)%post(kki) !the post form (either ligand or receptor)
!                          do k=1,d_nwpost(j,kki)
!                            kk=d_gen(int(d_wpost(j,kki,k,1)))%post(1) !the other bound form  ! RZ 17-11-14 INT
!                            krev= d_wpost(j,kki,k,2) !dissociation constant
!                          end do                           
!                          b=krev*d_gex(i,j)*d_gex(ie,kk)   !loss of bound form by dissociation
!                          did(kkk) = did(kkk) + b !this is the gain of ligand or receptor by dissociation of bound form
!                          did(j)   = did(j)   - b !this is the loss of bound form by dissociation (it loses 1 for every couple of ligand-receptor released)
!                        end do
!                      end if
!                      if (d_gen(j)%npre>0) then !the binding reaction
!                        do kki=1,d_gen(j)%npre   
!                          kkk=d_gen(j)%pre(kki) !the pre form (either ligand or receptor)
!                          do k=1,d_nwpre(j,kki)
!                            kk=d_wpre(j,kki,k,1) !catalist (either ligand or receptor) and also the one
!                            kbound=d_wpre(j,kki,k,2) !binding constant
!                            smd=d_gex(i,kkk)*d_gex(ie,kk) !productory of the concentrations of ligand and receptor (kk and kkk)
!                          end do
!                          b=kbound*smd         ! gain of bound form by binding of receptor and ligand
!                          did(kkk) = did(kkk) - b !this is the loss of dissociated form due to binding
!                          did(j)   = did(j)   + b !this is the gain of bound form
!                        end do
!                      end if
!                    end do
!                  end if
!                end if
!              else
!                if(dist.lt.nbh) then ! make sure neighbours within given distance are taken into account  
!                  if(tipi==d_node(ie)%tipus)then
!                    do jjj=1,d_nkindof(7)
!                      j=d_wkindof(7,jjj)
!                      did(j)=did(j)+d_gen(j)%diffu*(d_gex(ie,j)-d_gex(i,j))*udist
!                    end do
!                  end if
!                end if
!              end if
!            end if
!          end if
!          if (d_nkindof(8)>0) then  ! receptors of extracellular diffusible molecules ONLY THE ACTIVE FORM, the inactive is kindof2 or 3
!            if (celi==d_node(ie)%icel) then
!              if (tipi==d_node(ie)%tipus) then ! the active form diffuses only in the membrane (the inactive everywhere)
!                if(dist.lt.nbh) then ! make sure neighbours within given distance are taken into account   
!                  do jjj=1,d_nkindof(8)
!                    j=d_wkindof(8,jjj)
!                    did(j)=did(j)+d_gen(j)%diffu*(d_gex(ie,j)-d_gex(i,j))*udist
!                  end do
!                end if
!              end if
!            end if
!          end if
!        !end do
!        
!        if(i==1)then !this need only be done by one thread
!        !loss by diffusion when there are borders
!          if(d_node(i)%border==1.and. d_nkindof(4)>0)then
!            do jjj=1,d_nkindof(4) ! kindof 4 diffuses within the cell but it is only in its surface (so it is not within the cell really)
!              j=d_wkindof(4,jjj)
!              did(j)=did(j)-d_gen(j)%diffu*d_gex(i,j)
!            end do
!          end if
!        end if
!        
!        call syncthreads()
!        if(i==1) dif_dif(nod,1:ng)=did(1:ng) 
!
!457 continue
!
!end subroutine
!
!
!!*******************************************************************************************
!
!attributes(global) subroutine genestep_reaction_kernel(dif_dif,digex,ng)
!
!integer :: i, k, kk ! RZ 17-11-14 ! added k, kk
!!integer :: ii1,i1,ii2,i2,ii3,i3,iii1,iii2,iii3
!integer :: j,jj,jjj,kkk,iii
!real*8  :: a,b
!integer :: kki,ivv
!real*8  :: sm,smd,smp
!integer :: ie,tipi,celi
!!real*8  :: ix,iy,iz
!integer,value::ng
!real*8  :: gext(ng)!,agex(nda,ng)
!!real*8  :: dist,nbh,udist,udistn
!!integer :: nv,nbo
!real*8 :: did(ng)
!!integer :: checkn(nd)
!!integer :: c_nneigh(nd)
!real*8  ::krev,kbound
!
!real*8,dimension(:,:) :: dif_dif ! the increment: only for calculations
!real*8,dimension(:,:) :: digex ! the increment: only for calculations
!!integer,dimension(:) :: checkn
!
!       !if (allocated(dgex)) deallocate(dgex)
!       !allocate(dgex(nd,ng))
!       !dgex=0.0d0
!       !if (allocated(checkn)) deallocate(checkn)
!       !allocate(checkn(nd))
!       !checkn=0.0d0
!!      rdiffmax=2*maxval(node(:nd)%da)*dmax !now it's in the iteracio subroutine !>>Miquel27-2-14
!
!        !i=blockidx%x  !here each node will be a block
!        i=blockDim%x*(blockidx%x-1)+threadidx%x
!        if (i>d_nd) goto 457
!
!!print*,"debug: genestep_kernel tall0",i
!        !did(1:ng)=0.0d0
!        gext(1:ng)=0d0  !>>Miquel22-8-14
!
!        did(1:ng)=dif_dif(i,1:ng)
!!print*,"debug: ng",ng,gext(1:ng)
!
!        tipi=d_node(i)%tipus
!        celi=d_node(i)%icel !>>Miquel31-3-14
!        !if (rdiffmax<2*node(i)%da)   then ;print *,"RRRRRRRRRRRRR"; nbh=2*node(i)%da   ; else ; nbh=rdiffmax ; end if  ! this may make nbh different for different nodes 
!        !!! and then non-symmetric diffusion 
!        !nbh=d_rdiffmax          
!        !nv=0
!        !ix=d_node(i)%x     ; iy=d_node(i)%y     ; iz=d_node(i)%z   
!        !ii1=nint(iz*d_urv) ; ii2=nint(iy*d_urv) ; ii3=nint(ix*d_urv)
!        !ivv=d_node(i)%altre
!
!        !nbo=nint(nbh*d_urv)
!        !nneigh(i)=0 
!        !checkn(i)=0 ! RZ 4-3-14
!
!!print*,"debug: genestep_kernel tall1",i,"nneigh",d_nneigh(i)
!
!        ! ************************* STANDARD MINIMAL LOOP (1,2,3,4) **********************************
!        ! DIFFUSION
!
!        !REACTION 
!        if (d_node(i)%marge==0 .and. d_node(i)%tipus<4) then ! IN THE NUCLEUS (THERE WILL BE TRANSCRIPTION)
!          do jjj=1,d_nkindof(1)           !*****************PRODUCTION OF PRIMARY FORMS (WITHOUT FURTHER FORMS)
!            j=d_wkindof(1,jjj)     
!            a=d_gex(i,j)
!            sm=0.0d0
!            do k=1,int(d_nw(j))
!              kk=d_w(j,k)
!              sm = sm + d_gen(j)%w(kk)*d_gex(i,kk)  !the amount of TFs
!            end do
!            if (sm<0.0) sm=0.0
!            gext(j) = sm/(1+sm) - d_gen(j)%mu*a + did(j) !gain of product by transcription (Michaelis-Menten) !>>>> Miquel 20-12-13
!          end do 
!          do jjj=1,d_nkindof(2)           !*****************PRODUCTION OF PRIMARY FORMS (WITH FURTHER FORMS)
!            j=d_wkindof(2,jjj)
!            a=d_gex(i,j)
!            sm=0.0d0
!            do k=1,int(d_nw(j))
!              kk=d_w(j,k)
!              sm = sm + d_gen(j)%w(kk)*d_gex(i,kk)  !the amount of TFs
!            end do
!            if (sm<0.0) sm=0.0
!            if (d_gen(j)%npost>0) then 
!              do kki=1,d_gen(j)%npost   
!                kkk=d_gen(j)%post(kki); if(d_gen(kkk)%kindof==8) cycle !>>Miquel21-8-14
!                smd=0.0d0
!                do k=1,d_nwpost(j,kki)
!                  kk=int(d_wpost(j,kki,k,1)) ! RZ 17-11-14 INT
!                  smd=smd + d_wpost(j,kki,k,2)*d_gex(i,kk)   !loss of product by transforming into another form !>>> Miquel20-12-13
!                end do
!                if (smd<0.0) smd=0.0   ! Is 9-2-14
!                b=smd*a/(1+a)         ! Michaelis-Menten
!                gext(kkk) = gext(kkk) + b ! that is the gain by the post due to its transformation from j
!                gext(j)   = gext(j)   - b ! that is the lost j has due to its transformation to j
!              end do
!              gext(j) = gext(j) + sm/(1+sm) - d_gen(j)%mu*a + did(j) !>>> Is 9-2-14
!            else
!              gext(j) = sm/(1+sm) - d_gen(j)%mu*a + did(j)  !>>> Miquel20-12-13
!            end if
!          end do
!        else  !****************** NOT IN THE NUCLEUS (NO TRANSCRIPTION)
!
!          do jjj=1,d_nkindof(1)           !*****************PRODUCTION OF PRIMARY FORMS (WITHOUT FURTHER FORMS)
!            j=d_wkindof(1,jjj)  
!            gext(j) = - d_gen(j)%mu*d_gex(i,j) + did(j) !gain of product by transcription (Michaelis-Menten) !>>>> Miquel 20-12-13
!          end do 
!          do jjj=1,d_nkindof(2)           !*****************PRODUCTION OF PRIMARY FORMS (WITH FURTHER FORMS)
!            j=d_wkindof(2,jjj)     
!            a=d_gex(i,j)
!            if (d_gen(j)%npost>0) then !
!              do kki=1,d_gen(j)%npost
!                kkk=d_gen(j)%post(kki); if(d_gen(kkk)%kindof==8) cycle !>>Miquel21-8-14
!                smd=0.0d0
!                do k=1,d_nwpost(j,kki)
!                  kk=int(d_wpost(j,kki,k,1)) ! RZ 17-11-14 INT
!                  smd=smd + d_wpost(j,kki,k,2)*d_gex(i,kk)   !loss of product by transforming into another form !>>> Miquel20-12-13
!                end do                           
!                if (smd<0.0) smd=0.0   ! Is 9-2-14
!                b=smd*a/(1+a)         ! Michaelis-Menten
!                gext(kkk) = gext(kkk) + b ! that is the gain by the post due to its transformation from j
!                gext(j)   = gext(j)   - b ! that is the lost j has due to its transformation to j
!              end do
!              gext(j) = gext(j) - d_gen(j)%mu*a + did(j) !>>> Is 9-2-14
!            else
!              gext(j) = - d_gen(j)%mu*a + did(j)  !>>> Miquel20-12-13
!            end if
!          end do
!        end if
!        
!        !********PRODUCTION OF SECONDARY FORMS: they are never transcripts : so it is the same whether we are in the nucleus or not
!        do jjj=1,d_nkindof(3)           ! THIS IS IMPORTANT, IF YOU WANT TO MAKE A NOTCH YOU SHOULD MAKE ITS PREVIOUS FORM AND TRANSCRIBE IT
!          j=d_wkindof(3,jjj)
!          a=d_gex(i,j)
!          if (d_gen(j)%npost>0) then 
!            do kki=1,d_gen(j)%npost
!              kkk=d_gen(j)%post(kki); if(d_gen(kkk)%kindof==8) cycle !>>Miquel21-8-14
!              smd=0.0d0
!              do k=1,d_nwpost(j,kki)
!                kk=int(d_wpost(j,kki,k,1))  ! RZ 17-11-14 INT
!                smd=smd + d_wpost(j,kki,k,2)*d_gex(i,kk)   !loss of product by transforming into another form !>>> Miquel20-12-13
!              end do                           
!              if (smd<0.0) smd=0.0   ! Is 9-2-14
!              b=smd*a/(1+a)         ! Michaelis-Menten
!              gext(kkk) = gext(kkk) + b ! that is the gain by the post due to its transformation from j
!              gext(j)   = gext(j)   - b ! that is the lost j has due to its transformation to j
!            end do
!            !gext(j) = gext(j) - gen(j)%mu*a + did(j) !>>> Is 9-2-14
!          end if
!          gext(j) = gext(j)  - d_gen(j)%mu*a + did(j)  ! Is 9-2-14
!        end do 
!
!        !********PRODUCTION OF SECONDARY FORMS: they are never transcripts : so it is the same whether we are in the nucleus or not
!        do jjj=1,d_nkindof(4)           ! THIS IS IMPORTANT, IF YOU WANT TO MAKE A NOTCH YOU SHOULD MAKE ITS PREVIOUS FORM AND TRANSCRIBE IT
!          j=d_wkindof(4,jjj)     
!          a=d_gex(i,j)
!          if (d_gen(j)%npost>0) then 
!            do kki=1,d_gen(j)%npost   
!              kkk=d_gen(j)%post(kki); if(d_gen(kkk)%kindof==8) cycle !>>Miquel21-8-14
!              smd=0.0d0
!              do k=1,d_nwpost(j,kki)
!                kk=int(d_wpost(j,kki,k,1)) ! RZ 17-11-14 INT
!                smd=smd + d_wpost(j,kki,k,2)*d_gex(i,kk)   !loss of product by transforming into another form !>>> Miquel20-12-13
!              end do                           
!              if (smd<0.0) smd=0.0   ! Is 9-2-14
!              b=smd*a/(1+a)         ! Michaelis-Menten
!              gext(kkk) = gext(kkk) + b ! that is the gain by the post due to its transformation from j
!              gext(j)   = gext(j)   - b ! that is the lost j has due to its transformation to j
!            end do
!            !gext(j) = gext(j) - gen(j)%mu*a + did(j) !>>> Is 9-2-14
!          end if
!
!          ! NOTICE: growth factors do not return any thing to pre forms, there is not back reaction, because they get secreted as they get produced
!          gext(j) = gext(j)  - d_gen(j)%mu*a + did(j)  ! Is 9-2-14
!        end do 
!
!        do jj=5,7 !********PRODUCTION OF OTHER SECONDARY FORMS: they are never transcripts : so it is the same whether we are in the nucleus or not
!        !  if (jj==7) cycle   !>>> Is 1-3-14 ! outcommented RZ 4-3-14
!          if (jj==7) then ! >>> RZ 4-3-14
!            do jjj=1,d_nkindof(jj)
!              j=d_wkindof(jj,jjj)     
!              a=d_gex(i,j)
!              gext(j) = gext(j)  - d_gen(j)%mu*a + did(j)  ! Is 9-2-14
!            end do
!          else   ! <<< RZ 4-3-14
!            do jjj=1,d_nkindof(jj)           ! THIS IS IMPORTANT, IF YOU WANT TO MAKE A NOTCH YOU SHOULD MAKE ITS PREVIOUS FORM AND TRANSCRIBE IT
!              j=d_wkindof(jj,jjj)     
!              a=d_gex(i,j)
!              if (d_gen(j)%npost>0) then 
!                do kki=1,d_gen(j)%npost   
!                  kkk=d_gen(j)%post(kki); if(d_gen(kkk)%kindof==8) cycle !>>Miquel21-8-14
!                  smd=0.0d0
!                  do k=1,d_nwpost(j,kki)
!                    kk=int(d_wpost(j,kki,k,1)) ! RZ 17-11-14 INT
!                    smd=smd + d_wpost(j,kki,k,2)*d_gex(i,kk)   !loss of product by transforming into another form !>>> Miquel20-12-13
!                  end do                           
!                  if (smd<0.0) smd=0.0   ! Is 9-2-14
!                  b=smd*a/(1+a)         ! Michaelis-Menten
!                  gext(kkk) = gext(kkk) + b ! that is the gain by the post due to its transformation from j
!                  gext(j)   = gext(j)   - b ! that is the lost j has due to its transformation to j
!                end do
!                !gext(j) = gext(j) - gen(j)%mu*a + did(j) !>>> Is 9-2-14
!              end if
!              gext(j) = gext(j)  - d_gen(j)%mu*a + did(j)  ! Is 9-2-14
!            end do 
!          end if  ! RZ 4-3-14
!        end do
!
!        !bound receptors !>>Miquel18-8-14
!        do jjj=1,d_nkindof(8)
!          j=d_wkindof(8,jjj)
!          if (d_gen(j)%npost<1) cycle !they should though ! Is 1-10-14
!          a=d_gex(i,j)
!          !the dissociation reaction
!          krev=0d0
!          kkk=d_gen(j)%post(1) !the post form (either ligand or receptor)
!          iii=d_gen(j)%post(2) !the post form (either ligand or receptor)
!          krev= d_wpost(j,1,1,2) !dissociation constant
!          b=krev*a       !loss of bound form by dissociation
!          gext(kkk) = gext(kkk) + b !this is the gain of ligand or receptor by dissociation of bound form
!          gext(iii) = gext(iii) + b !this is the gain of ligand or receptor by dissociation of bound form
!
!          gext(j)   = gext(j)   - b !this is the loss of bound form by dissociation (it loses 1 for every couple of ligand-receptor released)
!
!          !the binding reaction
!          if (d_gen(j)%npre<1) cycle !they should though !>>Miquel1-10-14
!          kbound=0d0
!          kkk=d_gen(j)%pre(1) !the pre form (either ligand or receptor)
!          iii=d_gen(j)%pre(2) !the pre form (either ligand or receptor)
!          kbound=d_wpre(j,1,1,2) !binding constant
!          smd=d_gex(i,iii)*d_gex(i,kkk) !productory of the concentrations of ligand and receptor
!          b=kbound*smd         ! gain of bound form by binding of receptor and ligand
!
!          gext(kkk) = gext(kkk) - b !this is the loss of dissociated form due to binding
!          gext(iii) = gext(iii) - b !this is the loss of dissociated form due to binding
!
!          gext(j)   = gext(j)   + b !this is the gain of bound form due to binding
!          gext(j) = gext(j)  - d_gen(j)%mu*a + did(j)  ! Is 9-2-14
!        end do
!
!!print*,"debug: genestep_kernel tall4",i 
!        !do j=1,ng
!        !print*,"debug: genestep_kernel tall4.",j,gext(j)
!        !  d_dgex(i,j)=gext(j)
!        !end do
!        digex(i,1:ng)=gext(1:ng)  ! this is the increment in the gene this step
!!print*,"debug: genestep_kernel tall5",i 
!457 continue
!end subroutine genestep_reaction_kernel


!!!!!!!!!!!!!!!!!!

end module device
