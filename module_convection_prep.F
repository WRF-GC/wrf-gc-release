!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               __          _______  ______       _____  _____                 !
!               \ \        / /  __ \|  ____|     / ____|/ ____|                !
!                \ \  /\  / /| |__) | |__ ______| |  __| |                     !
!                 \ \/  \/ / |  _  /|  __|______| | |_ | |                     !
!                  \  /\  /  | | \ \| |         | |__| | |____                 !
!                   \/  \/   |_|  \_\_|          \_____|\_____|                !
!                                                                              !
!----------------------- ALPHA VERSION, v0.1 (20190101) -----------------------!
!
! WRF-GC: GEOS-Chem High Performance-powered Chemistry Add-On for WRF Model
! Developed by Xu Feng <fengx7@pku.edu.cn> (Convection module)
!    November 2018, Peking University, Dept of Atmospheric and Oceanic Sciences
!    Correspondence to: Tzung-May Fu <tmfu@pku.edu.cn>
!
! ALPHA INFORMATION:
!    WRF-GC Alpha (version 0.1) is experimental. Please notify the authors of
!    any bugs, suggestions and feature requests through email or the GitHub
!    repository.
!
! COPYRIGHT STATEMENT:
!    Permission is hereby granted, free of charge, to any person obtaining a copy
!   of this software and associated documentation files (the "Software"), to 
!   use, copy, modify the Software, and to permit persons to whom the Software is
!   furnished to do so, subject to the following conditions:
!
!   - The above copyright notice and this permission notice shall be included in all
!   copies or substantial portions of the Software.
! 
!   - The Software, modified in part or in full may not be redistributed without
!   express permission from the copyright holder(s).
! 
!   Except as contained in this notice or in attribution, the name of the WRF-GC model
!   shall not be used as an endorsement for distributing modified copies of the
!   Software without prior written permission from the copyright holder(s).
! 
!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!   SOFTWARE.
! 
!  WRF and the GEOS-Chem model, GCHP are (c) their original authors.
!
!------------------------------------------------------------------------------- 
module module_convection_prep
     use module_model_constants, only:rd=>r_d, rv=>r_v, &
   &       cpd=>cp, alv=>xlv, als=>xls, alf=>xlf, t13=>Prandtl, g              

     implicit none
     real,private :: rcpd,vtmpc1,tmelt,                &
             c1es,c2es,c3les,c3ies,c4les,c4ies,c5les,c5ies,zrg

     real,private :: r5alvcp,r5alscp,ralvdcp,ralsdcp,ralfdcp,rtwat,rtber,rtice 
     real,private :: entrdd,cmfcmax,cmfcmin,cmfdeps,zdnoprc,cprcon,pgcoef
     integer,private :: momtrans

     parameter(         &
      rcpd=1.0/cpd,     &
      tmelt=273.16,     &
      zrg=1.0/g,        &
      c1es=610.78,      &
      c2es=c1es*rd/rv,  &
      c3les=17.2693882, &
      c3ies=21.875,     &
      c4les=35.86,      &
      c4ies=7.66,       &
      c5les=c3les*(tmelt-c4les),  &
      c5ies=c3ies*(tmelt-c4ies),  &
      r5alvcp=c5les*alv*rcpd,     &
      r5alscp=c5ies*als*rcpd,     &
      ralvdcp=alv*rcpd,           &
      ralsdcp=als*rcpd,           &
      ralfdcp=alf*rcpd,           &
      rtwat=tmelt,                &
      rtber=tmelt-5.,             &
      rtice=tmelt-23.,            &
      vtmpc1=rv/rd-1.0 )
      parameter(entrdd = 2.0e-4)
      parameter(cmfcmax = 1.0)
      parameter(cmfcmin = 1.e-10)
      parameter(cmfdeps = 0.30)
      parameter(zdnoprc = 2.0e4)
      parameter(cprcon = 1.4e-3)
      parameter(momtrans = 2 )
      parameter(pgcoef=0.7)
      logical :: nonequil
      parameter(nonequil = .true. )
      logical :: lmfpen,lmfmid,lmfscv,lmfdd,lmfdudv
      parameter(lmfpen=.true.,lmfmid=.true.,lmfscv=.true.,lmfdd=.true.,lmfdudv=.true.)
contains
     subroutine cu_ntiedtke(                                    &
                 dt,itimestep,stepcu                            &
                ,raincv,pratec,qfx,hfx                          &
                ,u3d,v3d,w,t3d,qv3d,qc3d,qi3d,pi3d,rho3d        &
                ,qvften,thften                                  &
                ,dz8w,pcps,p8w,xland,cu_act_flag,dx             &
                ,ids,ide, jds,jde, kds,kde                      &
                ,ims,ime, jms,jme, kms,kme                      &
                ,its,ite, jts,jte, kts,kte                      &
                ,rthcuten,rqvcuten,rqccuten,rqicuten            &
                ,rucuten, rvcuten                               &

                ,pmflxrain, pmflxsnow                           &
                ,cmfmc, dqrcu, cloud_bot                        &
                ,dtrainu, dtraind, dtrain, reevapcn             &
                                                                )
      implicit none
      integer, intent(in) ::            ids,ide, jds,jde, kds,kde,      &
                                        ims,ime, jms,jme, kms,kme,      &
                                        its,ite, jts,jte, kts,kte,      &
                                        itimestep,                      &
                                        stepcu

      real,    intent(in) ::                                            &
                                        dt,                             &
                                        dx


      real,    dimension(ims:ime, jms:jme), intent(in) ::               &
                                        xland

      real,    dimension(ims:ime, jms:jme), intent(inout) ::            &
                                        raincv, pratec

      real,    dimension(ims:ime, kms:kme, jms:jme), intent(inout) ::   &
                                        pmflxrain, pmflxsnow, cmfmc

      real,    dimension(ims:ime, kms:kme, jms:jme), intent(inout) ::   &
                                        dqrcu
      integer, dimension(ims:ime, jms:jme), intent(inout) ::            &
                                        cloud_bot
      real,    dimension(ims:ime, kms:kme, jms:jme), intent(inout) ::   &
                                        dtrainu, dtraind, dtrain
      real,    dimension(ims:ime, kms:kme, jms:jme), intent(inout) ::   &
                                        reevapcn
      logical, dimension(ims:ime,jms:jme), intent(inout) ::             &
                                        cu_act_flag


      real,    dimension(ims:ime, kms:kme, jms:jme), intent(in) ::      &
                                        dz8w,                           &
                                        pcps,                           &
                                        p8w,                            &
                                        pi3d,                           &
                                        qc3d,                           &
                                        qvften,                         &
                                        thften,                         &
                                        qi3d,                           &
                                        qv3d,                           &
                                        rho3d,                          &
                                        t3d,                            &
                                        u3d,                            &
                                        v3d,                            &
                                        w
      real,    dimension(ims:ime, jms:jme) ::                           &
                                        qfx,                            &
                                        hfx                            

      real, dimension(ims:ime, kms:kme, jms:jme),                       &
               optional, intent(inout) ::                               &
                                        rqccuten,                       &
                                        rqicuten,                       &
                                        rqvcuten,                       &
                                        rthcuten,                       &
                                        rucuten,                        &
                                        rvcuten

      real      ::                                      &
                                        delt,                           &
                                        rdelt

      real     , dimension(its:ite) ::                  &
                                        rcs,                            &
                                        rn,                             &
                                        evap,                           &
                                        heatflux
      integer  , dimension(its:ite) ::  slimsk

      real     , dimension(its:ite, kts:kte+1) ::         &
                                        pmflxr1,          &
                                        pmflxs1
      real     , dimension(its:ite, kts:kte)   ::         &
                                        zmfu1,            &
                                        zmfd1
      integer  , dimension(its:ite) ::  icbot1
      
      real     , dimension(its:ite, kts:kte)  ::          &
                                        pmfude_rate1,      &
                                        pmfdde_rate1    
      real     , dimension(its:ite, kts:kte)  ::          &
                                        zdmfdp1
       
      real     , dimension(its:ite, kts:kte+1) ::         &
                                        prsi,           &
                                        ghti,           &
                                          zi 

      real     , dimension(its:ite, kts:kte) ::         &
                                        dot,                            &
                                        prsl,                           &
                                        q1,                             &
                                        q2,                             &
                                        q3,                             &
                                        q1b,                            &
                                        t1b,                            &
                                        q11,                            &
                                        q12,                            &
                                        t1,                             &
                                        u1,                             &
                                        v1,                             &
                                        zl,                             &
                                        omg,                            &
                                        ghtl                           

      integer, dimension(its:ite) ::                                    &
                                        kbot,                           &
                                        ktop

      integer ::                                                        &
                                        i,                              &
                                        im,                             &
                                        j,                              &
                                        k,                              &
                                        km,                             &
                                        kp,                             &
                                        kx,                             &
                                        kx1


      integer                      :: zz, pp
      do j=jts,jte
         do i=its,ite
            cu_act_flag(i,j)=.true.
         enddo
      enddo

      im=ite-its+1
      kx=kte-kts+1
      kx1=kx+1
      delt=dt*stepcu
      rdelt=1./delt

   do j=jts,jte

      do i=its,ite
        zi(i,kts)=0.0
      enddo
      do k=kts,kte
        do i=its,ite
          zi(i,k+1)=zi(i,k)+dz8w(i,k,j)
        enddo
      enddo
      do k=kts,kte
        do i=its,ite
          zl(i,k)=0.5*(zi(i,k)+zi(i,k+1))
        enddo
      enddo

      do i=its,ite
        slimsk(i)=int(abs(xland(i,j)-2.))
      enddo

      do k=kts,kte
        kp=k+1
        do i=its,ite
          dot(i,k)=-0.5*g*rho3d(i,k,j)*(w(i,k,j)+w(i,kp,j))
        enddo
      enddo

      pp = 0
      do k=kts,kte
        zz = kte-pp
        do i=its,ite
          u1(i,zz)=u3d(i,k,j)
          v1(i,zz)=v3d(i,k,j)
          t1(i,zz)=t3d(i,k,j)
          q1(i,zz)=qv3d(i,k,j)
          if(itimestep == 1) then
             q1b(i,zz)=0.
             t1b(i,zz)=0.
          else
             q1b(i,zz)=qvften(i,k,j)
             t1b(i,zz)=thften(i,k,j)
          endif
          q2(i,zz)=qc3d(i,k,j)
          q3(i,zz)=qi3d(i,k,j)
          omg(i,zz)=dot(i,k)
          ghtl(i,zz)=zl(i,k)
          prsl(i,zz) = pcps(i,k,j)
        enddo
        pp = pp + 1
      enddo

      pp = 0
      do k=kts,kte+1
        zz = kte+1-pp
        do i=its,ite
          ghti(i,zz) = zi(i,k)
          prsi(i,zz) = p8w(i,k,j) 
        enddo
        pp = pp + 1
      enddo
      do i=its,ite
        evap(i)    = qfx(i,j)
        heatflux(i)= hfx(i,j)
      enddo
      call tiecnvn(u1,v1,t1,q1,q2,q3,q1b,t1b,ghtl,ghti,omg,prsl,prsi,evap,heatflux,  &
                  rn,pmflxr1,pmflxs1,zmfu1,zmfd1,icbot1,pmfude_rate1,pmfdde_rate1,zdmfdp1,  &
                  slimsk,im,kx,kx1,delt,dx)

      do i=its,ite
         raincv(i,j)=rn(i)/stepcu
         pratec(i,j)=rn(i)/(stepcu * dt)
      enddo
    
      pp = 0 
      do k=kts,kte+1
      	zz = kte+1-pp
        do i=its,ite
          pmflxrain(i,k,j)=pmflxr1(i,zz)
          pmflxsnow(i,k,j)=pmflxs1(i,zz)
        enddo
        pp = pp + 1
      enddo
      do i=its,ite
         cloud_bot(i,j)=kte - icbot1(i) + 1
      enddo

      pp=0
      do k=kts,kte
        zz = kte-pp
        do i=its,ite
          cmfmc(i,k,j)=zmfu1(i,zz) + zmfd1(i,zz)

          if (k <= cloud_bot(i,j)) then
             cmfmc(i,k,j)=0.0
          endif

          cmfmc(i,k,j)=min(cmfmc(i,k,j),0.2)
          cmfmc(i,k,j)=max(cmfmc(i,k,j),0.0)
        enddo
        pp = pp + 1
      enddo
      pp=0
      do k=kts,kte
        zz = kte-pp
        do i=its,ite
          dtrainu(i,k,j)=pmfude_rate1(i,zz)
          dtraind(i,k,j)=pmfdde_rate1(i,zz)
          reevapcn(i,k,j)=zdmfdp1(i,zz)
        enddo
        pp = pp + 1
      enddo
      
      pp = 0
      do k=kts,kte
        zz = kte-pp
        do i=its,ite
          rthcuten(i,k,j)=(t1(i,zz)-t3d(i,k,j))/pi3d(i,k,j)*rdelt
          rqvcuten(i,k,j)=(q1(i,zz)-qv3d(i,k,j))*rdelt
          rucuten(i,k,j) =(u1(i,zz)-u3d(i,k,j))*rdelt
          rvcuten(i,k,j) =(v1(i,zz)-v3d(i,k,j))*rdelt
        enddo
        pp = pp + 1
      enddo

      if(present(rqccuten))then
          pp = 0
          do k=kts,kte
            zz = kte-pp
            do i=its,ite
              rqccuten(i,k,j)=(q2(i,zz)-qc3d(i,k,j))*rdelt
            enddo
            pp = pp + 1
          enddo
        endif

      if(present(rqicuten))then
          pp = 0
          do k=kts,kte
            zz = kte-pp
            do i=its,ite
              rqicuten(i,k,j)=(q3(i,zz)-qi3d(i,k,j))*rdelt
            enddo
            pp = pp + 1
          enddo
        endif


   enddo

      do j=jts,jte
         do k=kts,kte
            do i=its,ite
               dqrcu(i, k, j) = 0.
            enddo
         enddo
      enddo

      do j=jts,jte
         do k=kts,kte
            do i=its,ite
               if(k .eq. cloud_bot(i,j)) then
                 dqrcu(i, k, j) = 1.0
               endif
            enddo
          enddo
       enddo
      
      do j=jts,jte
         do k=kts,kte
            do i=its,ite
              dtrain(i,k,j) = dtrainu(i,k,j) - dtraind(i,k,j)

	      if (k <= cloud_bot(i,j)) then
                 dtrain(i,k,j)=0.0
              endif

              dtrain(i,k,j) = min(dtrain(i,k,j),0.2)
              dtrain(i,k,j) = max(dtrain(i,k,j),0.0)
            enddo
         enddo
      enddo
 
   end subroutine cu_ntiedtke
   subroutine ntiedtkeinit(rthcuten,rqvcuten,rqccuten,rqicuten,         &
                     rucuten,rvcuten,rthften,rqvften,                   &
                     restart,p_qc,p_qi,p_first_scalar,                  &
                     allowed_to_read,                                   &
                     ids, ide, jds, jde, kds, kde,                      &
                     ims, ime, jms, jme, kms, kme,                      &
                     its, ite, jts, jte, kts, kte)
   implicit none
   logical , intent(in)           ::  allowed_to_read,restart
   integer , intent(in)           ::  ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte
   integer , intent(in)           ::  p_first_scalar, p_qi, p_qc

   real,     dimension( ims:ime , kms:kme , jms:jme ) , intent(out) ::  &
                                                              rthcuten, &
                                                              rqvcuten, &
                                                              rqccuten, &
                                                              rqicuten, &
                                                              rucuten,rvcuten,&
                                                              rthften,rqvften

   integer :: i, j, k, itf, jtf, ktf

   jtf=min0(jte,jde-1)
   ktf=min0(kte,kde-1)
   itf=min0(ite,ide-1)

   if(.not.restart)then
     do j=jts,jtf
     do k=kts,ktf
     do i=its,itf
       rthcuten(i,k,j)=0.
       rqvcuten(i,k,j)=0.
       rucuten(i,k,j)=0.
       rvcuten(i,k,j)=0.
     enddo
     enddo
     enddo

     DO j=jts,jtf
     DO k=kts,ktf
     DO i=its,itf
        rthften(i,k,j)=0.
        rqvften(i,k,j)=0.
     ENDDO
     ENDDO
     ENDDO

     if (p_qc .ge. p_first_scalar) then
        do j=jts,jtf
        do k=kts,ktf
        do i=its,itf
           rqccuten(i,k,j)=0.
        enddo
        enddo
        enddo
     endif

     if (p_qi .ge. p_first_scalar) then
        do j=jts,jtf
        do k=kts,ktf
        do i=its,itf
           rqicuten(i,k,j)=0.
        enddo
        enddo
        enddo
     endif
   endif

   end subroutine ntiedtkeinit

      subroutine tiecnvn(pu,pv,pt,pqv,pqc,pqi,pqvf,ptf,poz,pzz,pomg, &
     &         pap,paph,evap,hfx,zprecc,pmflxr,pmflxs,zmfu,zmfd,icbot, &
     &         pmfude_rate,pmfdde_rate,zdmfdp,lndj,lq,km,km1,dt,dx)

      implicit none

      real pu(lq,km),   pv(lq,km),   pt(lq,km),   pqv(lq,km)
      real poz(lq,km),  pomg(lq,km), evap(lq),    zprecc(lq)
      real pzz(lq,km1)

      real pum1(lq,km),    pvm1(lq,km),  ztt(lq,km),                &
     &     ptte(lq,km),    pqte(lq,km),  pvom(lq,km),  pvol(lq,km), &
     &     pverv(lq,km),   pgeo(lq,km),  pap(lq,km),   paph(lq,km1)
      real pqhfl(lq),      zqq(lq,km),    &
     &     prsfc(lq),      pssfc(lq),    pcte(lq,km), &
     &     phhfl(lq),      hfx(lq),      pgeoh(lq,km1)
      real ztp1(lq,km),    zqp1(lq,km),  ztu(lq,km),   zqu(lq,km),  &
     &     zlu(lq,km),     zlude(lq,km), zmfu(lq,km),  zmfd(lq,km), &
     &     zqsat(lq,km),   pqc(lq,km),   pqi(lq,km),   zrain(lq)
      real pqvf(lq,km),    ptf(lq,km)
      real pmflxr(lq,km1), pmflxs(lq,km1)
      real pmfude_rate(lq,km), pmfdde_rate(lq,km)
      real zdmfdp(lq,km)


      integer icbot(lq),   ictop(lq),     ktype(lq),   lndj(lq)
      logical locum(lq)

      real ztmst,fliq,fice,ztc,zalf,tt
      integer i,j,k,lq,km,km1
      real dt,dx,ztpp1
      real zew,zqs,zcor

      ztmst=dt

      do j=1,lq
        zrain(j)=0.0
        locum(j)=.false.
        prsfc(j)=0.0
        pssfc(j)=0.0
        pqhfl(j)=evap(j)
        phhfl(j)=hfx(j)
        pgeoh(j,km1)=g*pzz(j,km1)
      end do

      do k=1,km
        do j=1,lq
          pcte(j,k)=0.0
          pvom(j,k)=0.0
          pvol(j,k)=0.0
          ztp1(j,k)=pt(j,k)
          zqp1(j,k)=pqv(j,k)/(1.0+pqv(j,k))
          pum1(j,k)=pu(j,k)
          pvm1(j,k)=pv(j,k)
          pverv(j,k)=pomg(j,k)
          pgeo(j,k)=g*poz(j,k)
          pgeoh(j,k)=g*pzz(j,k)
          tt=ztp1(j,k)
          zew  = foeewm(tt)                                      
          zqs  = zew/pap(j,k)                                      
          zqs  = min(0.5,zqs)                                         
          zcor = 1./(1.-vtmpc1*zqs)                                    
          zqsat(j,k)=zqs*zcor      
          pqte(j,k)=pqvf(j,k) 
          zqq(j,k) =pqte(j,k)
          ptte(j,k)=ptf(j,k)
          ztt(j,k) =ptte(j,k)
        end do
      end do

      call cumastrn        &
     &    (lq,       km,       km1,      km-1,    ztp1, &
     &     zqp1,     pum1,     pvm1,     pverv,   zqsat,&
     &     pqhfl,    ztmst,    pap,      paph,    pgeo, &
     &     ptte,     pqte,     pvom,     pvol,    prsfc,&
     &     pssfc,    locum,    pmflxr,   pmflxs,        &
     &     pmfude_rate, pmfdde_rate,     zdmfdp,        &
     &     ktype,    icbot,    ictop,    ztu,     zqu,  &
     &     zlu,      zlude,    zmfu,     zmfd,    zrain,&
     &     pcte,     phhfl,    lndj,     pgeoh,   dx)

      do k=1,km
      do j=1,lq
      if(pcte(j,k).gt.0.) then
        fliq=foealfa(ztp1(j,k))
        fice=1.0-fliq
        pqc(j,k)=pqc(j,k)+fliq*pcte(j,k)*ztmst
        pqi(j,k)=pqi(j,k)+fice*pcte(j,k)*ztmst
      endif
      end do
      end do

      do k=1,km
        do j=1,lq
          pt(j,k)=  ztp1(j,k)+(ptte(j,k)-ztt(j,k))*ztmst
          zqp1(j,k)=zqp1(j,k)+(pqte(j,k)-zqq(j,k))*ztmst
          pqv(j,k)=zqp1(j,k)/(1.0-zqp1(j,k))
        end do
      end do

      do j=1,lq
        zprecc(j)=amax1(0.0,(prsfc(j)+pssfc(j))*ztmst)
      end do

      if (lmfdudv) then
        do k=1,km
          do j=1,lq
            pu(j,k)=pu(j,k)+pvom(j,k)*ztmst
            pv(j,k)=pv(j,k)+pvol(j,k)*ztmst
          end do
        end do
      endif

      return
      end subroutine tiecnvn

      subroutine cumastrn  &
     &    (klon,     klev,     klevp1,   klevm1,   pten, &
     &     pqen,     puen,     pven,     pverv,    pqsen,&
     &     pqhfl,    ztmst,    pap,      paph,     pgeo, &
     &     ptte,     pqte,     pvom,     pvol,     prsfc,& 
     &     pssfc,    ldcum,    pmflxr,   pmflxs,         &
     &     pmfude_rate, pmfdde_rate,     zdmfdp,         &
     &     ktype,    kcbot,    kctop,    ptu,      pqu,  &
     &     plu,      plude,    pmfu,     pmfd,     prain,&
     &     pcte,     phhfl,    lndj,     zgeoh,   dx)
      implicit none

















































      integer  klev,klon,klevp1,klevm1
      real     pten(klon,klev),        pqen(klon,klev),& 
     &         puen(klon,klev),        pven(klon,klev),&
     &         ptte(klon,klev),        pqte(klon,klev),&
     &         pvom(klon,klev),        pvol(klon,klev),&
     &         pqsen(klon,klev),       pgeo(klon,klev),&
     &         pap(klon,klev),         paph(klon,klevp1),&
     &         pverv(klon,klev),       pqhfl(klon),&
     &         phhfl(klon)            
      real     ptu(klon,klev),         pqu(klon,klev),&
     &         plu(klon,klev),         plude(klon,klev),&
     &         pmfu(klon,klev),        pmfd(klon,klev),&
     &         prain(klon),&
     &         prsfc(klon),            pssfc(klon)
      real     ztenh(klon,klev),       zqenh(klon,klev),&
     &         zgeoh(klon,klevp1),     zqsenh(klon,klev),&
     &         ztd(klon,klev),         zqd(klon,klev),&
     &         zmfus(klon,klev),       zmfds(klon,klev),&
     &         zmfuq(klon,klev),       zmfdq(klon,klev),&
     &         zdmfup(klon,klev),      zdmfdp(klon,klev),&
     &         zmful(klon,klev),       zrfl(klon),&
     &         zuu(klon,klev),         zvu(klon,klev),&
     &         zud(klon,klev),         zvd(klon,klev),&
     &         zlglac(klon,klev)
      real     pmflxr(klon,klevp1),    pmflxs(klon,klevp1)
      real     zhcbase(klon),& 
     &         zmfub(klon),            zmfub1(klon),&
     &         zdhpbl(klon)          
      real     zsfl(klon),             zdpmel(klon,klev),&
     &         pcte(klon,klev),        zcape(klon),&
     &         zcape1(klon),           zcape2(klon),&
     &         ztauc(klon),            ztaubl(klon),&
     &         zheat(klon)            
      real     wup(klon),              zdqcv(klon)            
      real     wbase(klon),            zmfuub(klon)
      real     upbl(klon)
      real     dx
      real     pmfude_rate(klon,klev), pmfdde_rate(klon,klev)
      real     zmfuus(klon,klev),      zmfdus(klon,klev)
      real     zuv2(klon,klev),ztenu(klon,klev),ztenv(klon,klev)
      real     zmfuvb(klon),zsum12(klon),zsum22(klon)
      integer  ilab(klon,klev),        idtop(klon),&
     &         ictop0(klon),           ilwmin(klon)
      integer  kdpl(klon)
      integer  kcbot(klon),            kctop(klon),&
     &         ktype(klon),            lndj(klon)
      logical  ldcum(klon)
      logical  loddraf(klon),          llo1,   llo2(klon)

      real     zcons,zcons2,zqumqe,zdqmin,zdh,zmfmax
      real     zalfaw,zalv,zqalv,zc5ldcp,zc4les,zhsat,zgam,zzz,zhhat
      real     zpbmpt,zro,zdz,zdp,zeps,zfac,wspeed
      integer  jl,jk,ik
      integer  ikb,ikt,icum,itopm2
      real     ztmst,ztau,zerate,zderate,zmfa
      real     zmfs(klon),pmean(klev),zlon
      real     zduten,zdvten,ztdis,pgf_u,pgf_v



      zcons=1./(g*ztmst)
      zcons2=3./(g*ztmst)



      call cuinin &
     &    (klon,     klev,     klevp1,   klevm1,   pten, &
     &     pqen,     pqsen,    puen,     pven,     pverv,&
     &     pgeo,     paph,     zgeoh,    ztenh,    zqenh,&
     &     zqsenh,   ilwmin,   ptu,      pqu,      ztd,  &
     &     zqd,      zuu,      zvu,      zud,      zvd,  &
     &     pmfu,     pmfd,     zmfus,    zmfds,    zmfuq,&
     &     zmfdq,    zdmfup,   zdmfdp,   zdpmel,   plu,  &
     &     plude,    ilab)






       call cutypen &
     &     ( klon,     klev,     klevp1,   klevm1,     pqen,&
     &      ztenh,    zqenh,     zqsenh,    zgeoh,     paph,&
     &      phhfl,    pqhfl,       pgeo,    pqsen,      pap,&
     &       pten,     lndj,        ptu,      pqu,     ilab,&
     &      ldcum,    kcbot,     ictop0,    ktype,    wbase,    plu,   kdpl)


       do jl=1,klon
         zdhpbl(jl)=0.0
         upbl(jl) = 0.0
         idtop(jl)=0
       end do

       do jk=2,klev
       do jl=1,klon
         if(jk.ge.kcbot(jl) .and. ldcum(jl)) then
            zdhpbl(jl)=zdhpbl(jl)+(alv*pqte(jl,jk)+cpd*ptte(jl,jk))&
     &                 *(paph(jl,jk+1)-paph(jl,jk))
            if(lndj(jl) .eq. 0) then
              wspeed = sqrt(puen(jl,jk)**2 + pven(jl,jk)**2) 
              upbl(jl) = upbl(jl) + wspeed*(paph(jl,jk+1)-paph(jl,jk))
            end if
         end if
       end do
       end do

      do jl=1,klon
        if(ldcum(jl)) then
           ikb=kcbot(jl)
           zmfmax = (paph(jl,ikb)-paph(jl,ikb-1))*zcons2
           if(ktype(jl) == 1) then
             zmfub(jl)= 0.1*zmfmax
           else if ( ktype(jl) == 2 ) then
             zqumqe = pqu(jl,ikb) + plu(jl,ikb) - zqenh(jl,ikb)
             zdqmin = max(0.01*zqenh(jl,ikb),1.e-10)
             zdh = cpd*(ptu(jl,ikb)-ztenh(jl,ikb)) + alv*zqumqe
             zdh = g*max(zdh,1.e5*zdqmin)
             if ( zdhpbl(jl) > 0. ) then
               zmfub(jl) = zdhpbl(jl)/zdh
               zmfub(jl) = min(zmfub(jl),zmfmax)
             else
               zmfub(jl) = 0.1*zmfmax
               ldcum(jl) = .false.
             end if
            end if  
        else
           zmfub(jl) = 0.
        end if
      end do





      call cuascn &
     &    (klon,     klev,     klevp1,   klevm1,   ztenh,&
     &     zqenh,    puen,     pven,     pten,     pqen,&
     &     pqsen,    pgeo,     zgeoh,    pap,      paph,&
     &     pqte,     pverv,    ilwmin,   ldcum,    zhcbase,&
     &     ktype,    ilab,     ptu,      pqu,      plu,&
     &     zuu,      zvu,      pmfu,     zmfub,&
     &     zmfus,    zmfuq,    zmful,    plude,    zdmfup,&
     &     kcbot,    kctop,    ictop0,   icum,     ztmst,&
     &     zqsenh,   zlglac,   lndj,     wup,      wbase,   kdpl,  pmfude_rate )



      do jl=1,klon
        if ( ldcum(jl) ) then
          ikb = kcbot(jl)
          itopm2 = kctop(jl)
          zpbmpt = paph(jl,ikb) - paph(jl,itopm2)
          if ( ktype(jl) == 1 .and. zpbmpt <  zdnoprc ) ktype(jl) = 2
          if ( ktype(jl) == 2 .and. zpbmpt >= zdnoprc ) ktype(jl) = 1
          ictop0(jl) = kctop(jl)
        end if
        zrfl(jl)=zdmfup(jl,1)
      end do

      do jk=2,klev
        do jl=1,klon
          zrfl(jl)=zrfl(jl)+zdmfup(jl,jk)
        end do
      end do

      do jk = 1,klev
      do jl = 1,klon
        pmfd(jl,jk) = 0.
        zmfds(jl,jk) = 0.
        zmfdq(jl,jk) = 0.
        zdmfdp(jl,jk) = 0.
        zdpmel(jl,jk) = 0.
      end do
      end do



      if(lmfdd) then


        call cudlfsn    &
     &    (klon,     klev,&
     &     kcbot,    kctop,    lndj,   ldcum,   &
     &     ztenh,    zqenh,    puen,     pven,   &
     &     pten,     pqsen,    pgeo,              &
     &     zgeoh,    paph,     ptu,      pqu,      plu,   &
     &     zuu,      zvu,      zmfub,    zrfl,             &
     &     ztd,      zqd,      zud,      zvd,               &
     &     pmfd,     zmfds,    zmfdq,    zdmfdp,             &
     &     idtop,    loddraf)


        call cuddrafn                                         &
     &   ( klon,     klev,     loddraf,                      &
     &     ztenh,    zqenh,    puen,     pven,                 &
     &     pgeo,     zgeoh,    paph,     zrfl,                  &
     &     ztd,      zqd,      zud,      zvd,      pmfu,         &
     &     pmfd,     zmfds,    zmfdq,    zdmfdp,   pmfdde_rate      )

      end if







      do jl=1,klon
      if(ldcum(jl) .and. ktype(jl) .eq. 1) then
        ikb = kcbot(jl)
        ikt = kctop(jl)
        zheat(jl)=0.0
        zcape(jl)=0.0
        zcape1(jl)=0.0
        zcape2(jl)=0.0
        zmfub1(jl)=zmfub(jl)
    
        ztauc(jl)  = (zgeoh(jl,ikt)-zgeoh(jl,ikb)) / &
                   ((2.+ min(15.0,wup(jl)))*g)
        if(lndj(jl) .eq. 0) then 
          upbl(jl) = 2.+ upbl(jl)/(paph(jl,klev+1)-paph(jl,ikb))
          ztaubl(jl) = (zgeoh(jl,ikb)-zgeoh(jl,klev+1))/(g*upbl(jl))
          ztaubl(jl) = min(300., ztaubl(jl))
        else
          ztaubl(jl) = ztauc(jl)
        end if
      end if    
      end do

      do jk = 1 , klev
      do jl = 1 , klon
        llo1 = ldcum(jl) .and. ktype(jl) .eq. 1
        if ( llo1 .and. jk <= kcbot(jl) .and. jk > kctop(jl) ) then
          ikb = kcbot(jl)
          zdz = pgeo(jl,jk-1)-pgeo(jl,jk)
          zdp = pap(jl,jk)-pap(jl,jk-1)
          zheat(jl) = zheat(jl) + ((pten(jl,jk-1)-pten(jl,jk)+zdz*rcpd) / &
                      ztenh(jl,jk)+vtmpc1*(pqen(jl,jk-1)-pqen(jl,jk))) * &
                      (g*(pmfu(jl,jk)+pmfd(jl,jk)))
          zcape1(jl) = zcape1(jl) + ((ptu(jl,jk)-ztenh(jl,jk))/ztenh(jl,jk) + &
                      vtmpc1*(pqu(jl,jk)-zqenh(jl,jk))-plu(jl,jk))*zdp
        end if

        if ( llo1 .and. jk >= kcbot(jl) ) then
        if((paph(jl,klev+1)-paph(jl,kdpl(jl)))<50.e2) then
          zdp = paph(jl,jk+1)-paph(jl,jk)
          zcape2(jl) = zcape2(jl) + ztaubl(jl)* &
                     ((1.+vtmpc1*pqen(jl,jk))*ptte(jl,jk)+vtmpc1*pten(jl,jk)*pqte(jl,jk))*zdp 
        end if
        end if
      end do
      end do

      do jl=1,klon
       if(ldcum(jl).and.ktype(jl).eq.1) then
           ikb = kcbot(jl)
           ikt = kctop(jl)
           ztau = ztauc(jl) * (1.+1.33e-5*dx)
           ztau = max(ztmst,ztau)
           ztau = max(360.,ztau)
           ztau = min(10800.,ztau)
           if(nonequil) then
             zcape2(jl)= max(0.,zcape2(jl))
             zcape(jl) = max(0.,min(zcape1(jl)-zcape2(jl),5000.))
           else
             zcape(jl) = max(0.,min(zcape1(jl),5000.))
           end if
           zheat(jl) = max(1.e-4,zheat(jl))
           zmfub1(jl) = (zcape(jl)*zmfub(jl))/(zheat(jl)*ztau)
           zmfub1(jl) = max(zmfub1(jl),0.001)
           zmfmax=(paph(jl,ikb)-paph(jl,ikb-1))*zcons2
           zmfub1(jl)=min(zmfub1(jl),zmfmax)
       end if
      end do




       do jl=1,klon
         if(ldcum(jl) .and. ktype(jl) .eq. 2) then
           ikb=kcbot(jl)
           if(pmfd(jl,ikb).lt.0.0 .and. loddraf(jl)) then
              zeps=-pmfd(jl,ikb)/max(zmfub(jl),cmfcmin)
           else
              zeps=0.
           endif
           zqumqe=pqu(jl,ikb)+plu(jl,ikb)-  &
     &            zeps*zqd(jl,ikb)-(1.-zeps)*zqenh(jl,ikb)
           zdqmin=max(0.01*zqenh(jl,ikb),cmfcmin)
           zmfmax=(paph(jl,ikb)-paph(jl,ikb-1))*zcons2

           zdh=cpd*(ptu(jl,ikb)-zeps*ztd(jl,ikb)- &
     &       (1.-zeps)*ztenh(jl,ikb))+alv*zqumqe
           zdh=g*max(zdh,1.e5*zdqmin)
           if(zdhpbl(jl).gt.0.)then
             zmfub1(jl)=zdhpbl(jl)/zdh
           else
             zmfub1(jl) = zmfub(jl)
           end if
           zmfub1(jl) = min(zmfub1(jl),zmfmax)
         end if


         if(ldcum(jl) .and. ktype(jl) .eq. 3 ) then
            zmfub1(jl) = zmfub(jl)
         end if

       end do


       do jk=1,klev
       do jl=1,klon
        if( ldcum(jl) ) then
           zfac=zmfub1(jl)/max(zmfub(jl),cmfcmin)
           pmfd(jl,jk)=pmfd(jl,jk)*zfac
           zmfds(jl,jk)=zmfds(jl,jk)*zfac
           zmfdq(jl,jk)=zmfdq(jl,jk)*zfac
           zdmfdp(jl,jk)=zdmfdp(jl,jk)*zfac
           pmfdde_rate(jl,jk) = pmfdde_rate(jl,jk)*zfac
        end if
       end do
       end do


    do jl = 1,klon
      if ( ldcum(jl) ) zmfs(jl) = zmfub1(jl)/max(cmfcmin,zmfub(jl))
    end do
    do jk = 2 , klev
      do jl = 1,klon
        if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
          ikb = kcbot(jl)
          if ( jk>ikb ) then
            zdz = ((paph(jl,klev+1)-paph(jl,jk))/(paph(jl,klev+1)-paph(jl,ikb)))
            pmfu(jl,jk) = pmfu(jl,ikb)*zdz
          end if
          zmfmax = (paph(jl,jk)-paph(jl,jk-1))*zcons2
          if ( pmfu(jl,jk)*zmfs(jl) > zmfmax ) then
            zmfs(jl) = min(zmfs(jl),zmfmax/pmfu(jl,jk))
          end if
        end if
      end do
    end do
    do jk = 2 , klev
      do jl = 1,klon
        if ( ldcum(jl) .and. jk <= kcbot(jl) .and. jk >= kctop(jl)-1 ) then
          pmfu(jl,jk) = pmfu(jl,jk)*zmfs(jl)
          zmfus(jl,jk) = zmfus(jl,jk)*zmfs(jl)
          zmfuq(jl,jk) = zmfuq(jl,jk)*zmfs(jl)
          zmful(jl,jk) = zmful(jl,jk)*zmfs(jl)
          zdmfup(jl,jk) = zdmfup(jl,jk)*zmfs(jl)
          plude(jl,jk) = plude(jl,jk)*zmfs(jl)
          pmfude_rate(jl,jk) = pmfude_rate(jl,jk)*zmfs(jl)
        end if
      end do
    end do


    do jl = 1,klon
      if ( ktype(jl) == 2 .and. &
           kcbot(jl) == kctop(jl) .and. kcbot(jl) >= klev-1 ) then
        ldcum(jl) = .false.
        ktype(jl) = 0
      end if
    end do

    if ( .not. lmfscv .or. .not. lmfpen ) then
      do jl = 1,klon
        llo2(jl) = .false.
        if ( (.not. lmfscv .and. ktype(jl) == 2) .or. &
             (.not. lmfpen .and. ktype(jl) == 1) ) then
          llo2(jl) = .true.
          ldcum(jl) = .false.
        end if
      end do
    end if


    do jl = 1,klon
      if ( loddraf(jl) .and. idtop(jl) <= kctop(jl) ) then
        idtop(jl) = kctop(jl) + 1
      end if
    end do
    do jk = 2 , klev
      do jl = 1,klon
        if ( loddraf(jl) ) then
          if ( jk < idtop(jl) ) then
            pmfd(jl,jk) = 0.
            zmfds(jl,jk) = 0.
            zmfdq(jl,jk) = 0.
            pmfdde_rate(jl,jk) = 0.
            zdmfdp(jl,jk) = 0.
          else if ( jk == idtop(jl) ) then
            pmfdde_rate(jl,jk) = 0.
          end if
        end if
      end do
    end do



       call cuflxn                                        &                  
     &  (  klon,     klev,     ztmst                       &                                       
     &  ,  pten,     pqen,     pqsen,    ztenh,   zqenh     &                  
     &  ,  paph,     pap,      zgeoh,    lndj,    ldcum      &                 
     &  ,  kcbot,    kctop,    idtop,    itopm2               &                 
     &  ,  ktype,    loddraf                                   &                 
     &  ,  pmfu,     pmfd,     zmfus,    zmfds                  &               
     &  ,  zmfuq,    zmfdq,    zmful,    plude                   &              
     &  ,  zdmfup,   zdmfdp,   zdpmel,   zlglac                   &             
     &  ,  prain,    pmfdde_rate, pmflxr, pmflxs )     

    do jl=1,klon
      zmfs(jl) = 1.
      zmfuub(jl)=0.
    end do
    do jk = 2 , klev
      do jl = 1,klon
        if ( loddraf(jl) .and. jk >= idtop(jl)-1 ) then
          zmfmax = pmfu(jl,jk)*0.98
          if ( pmfd(jl,jk)+zmfmax+1.e-15 < 0. ) then
            zmfs(jl) = min(zmfs(jl),-zmfmax/pmfd(jl,jk))
          end if
        end if
      end do
    end do

    do jk = 2 , klev
      do jl = 1 , klon
        if ( zmfs(jl) < 1. .and. jk >= idtop(jl)-1 ) then
          pmfd(jl,jk) = pmfd(jl,jk)*zmfs(jl)
          zmfds(jl,jk) = zmfds(jl,jk)*zmfs(jl)
          zmfdq(jl,jk) = zmfdq(jl,jk)*zmfs(jl)
          pmfdde_rate(jl,jk) = pmfdde_rate(jl,jk)*zmfs(jl)
          zmfuub(jl) = zmfuub(jl) - (1.-zmfs(jl))*zdmfdp(jl,jk)
          pmflxr(jl,jk+1) = pmflxr(jl,jk+1) + zmfuub(jl)
          zdmfdp(jl,jk) = zdmfdp(jl,jk)*zmfs(jl)
        end if
      end do
    end do

    do jk = 2 , klev - 1
      do jl = 1, klon
        if ( loddraf(jl) .and. jk >= idtop(jl)-1 ) then
          zerate = -pmfd(jl,jk) + pmfd(jl,jk-1) + pmfdde_rate(jl,jk)
          if ( zerate < 0. ) then
            pmfdde_rate(jl,jk) = pmfdde_rate(jl,jk) - zerate
          end if
        end if
        if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
          zerate = pmfu(jl,jk) - pmfu(jl,jk+1) + pmfude_rate(jl,jk)
          if ( zerate < 0. ) then
            pmfude_rate(jl,jk) = pmfude_rate(jl,jk) - zerate
          end if
          zdmfup(jl,jk) = pmflxr(jl,jk+1) + pmflxs(jl,jk+1) - &
                          pmflxr(jl,jk) - pmflxs(jl,jk)
          zdmfdp(jl,jk) = 0.
        end if
      end do
    end do

    do jl = 1,klon
      if ( loddraf(jl) ) then
        jk = idtop(jl)
        ik = min(jk+1,klev)
        if ( zmfdq(jl,jk) < 0.3*zmfdq(jl,ik) ) then
            zmfdq(jl,jk) = 0.3*zmfdq(jl,ik)
        end if
      end if
    end do


    do jk = 2 , klev
      do jl = 1, klon
        if ( ldcum(jl) .and. jk >= kctop(jl)-1 .and. jk < kcbot(jl) ) then
          zdz = ztmst*g/(paph(jl,jk+1)-paph(jl,jk))
          zmfa = zmfuq(jl,jk+1) + zmfdq(jl,jk+1) - &
                 zmfuq(jl,jk) - zmfdq(jl,jk) + &
                 zmful(jl,jk+1) - zmful(jl,jk) + zdmfup(jl,jk)
          zmfa = (zmfa-plude(jl,jk))*zdz
          if ( pqen(jl,jk)+zmfa < 0. ) then
            plude(jl,jk) = plude(jl,jk) + 2.*(pqen(jl,jk)+zmfa)/zdz
          end if
          if ( plude(jl,jk) < 0. ) plude(jl,jk) = 0.
        end if
        if ( .not. ldcum(jl) ) pmfude_rate(jl,jk) = 0.
        if ( abs(pmfd(jl,jk-1)) < 1.0e-20 ) pmfdde_rate(jl,jk) = 0.
      end do
    end do

    do jl=1,klon
      prsfc(jl) = pmflxr(jl,klev+1)
      pssfc(jl) = pmflxs(jl,klev+1)
    end do



      call cudtdqn(klon,klev,itopm2,kctop,idtop,ldcum,loddraf, &
                 ztmst,paph,zgeoh,pgeo,pten,ztenh,pqen,zqenh,pqsen,     &
                 zlglac,plude,pmfu,pmfd,zmfus,zmfds,zmfuq,zmfdq,zmful,   &
                 zdmfup,zdmfdp,zdpmel,ptte,pqte,pcte)



      if(lmfdudv) then
      do jk = klev-1 , 2 , -1
        ik = jk + 1
        do jl = 1,klon
          if ( ldcum(jl) ) then
            if ( jk == kcbot(jl) .and. ktype(jl) < 3 ) then
              ikb = kdpl(jl)
              zuu(jl,jk) = puen(jl,ikb-1)
              zvu(jl,jk) = pven(jl,ikb-1)
            else if ( jk == kcbot(jl) .and. ktype(jl) == 3 ) then
              zuu(jl,jk) = puen(jl,jk-1)
              zvu(jl,jk) = pven(jl,jk-1)
            end if
            if ( jk < kcbot(jl) .and. jk >= kctop(jl) ) then
            if(momtrans .eq. 1)then
              zfac = 0.
              if ( ktype(jl) == 1 .or. ktype(jl) == 3 ) zfac = 2.
              if ( ktype(jl) == 1 .and. jk <= kctop(jl)+2 ) zfac = 3.
              zerate = pmfu(jl,jk) - pmfu(jl,ik) + &
                (1.+zfac)*pmfude_rate(jl,jk)
              zderate = (1.+zfac)*pmfude_rate(jl,jk)
              zmfa = 1./max(cmfcmin,pmfu(jl,jk))
              zuu(jl,jk) = (zuu(jl,ik)*pmfu(jl,ik) + &
                zerate*puen(jl,jk)-zderate*zuu(jl,ik))*zmfa
              zvu(jl,jk) = (zvu(jl,ik)*pmfu(jl,ik) + &
                zerate*pven(jl,jk)-zderate*zvu(jl,ik))*zmfa
            else
              pgf_u = -pgcoef*0.5*(pmfu(jl,ik)*(puen(jl,ik)-puen(jl,jk))+&
                                   pmfu(jl,jk)*(puen(jl,jk)-puen(jl,jk-1)))
              pgf_v = -pgcoef*0.5*(pmfu(jl,ik)*(pven(jl,ik)-pven(jl,jk))+&
                                   pmfu(jl,jk)*(pven(jl,jk)-pven(jl,jk-1)))
              zerate = pmfu(jl,jk) - pmfu(jl,ik) + pmfude_rate(jl,jk)
              zderate = pmfude_rate(jl,jk)
              zmfa = 1./max(cmfcmin,pmfu(jl,jk))
              zuu(jl,jk) = (zuu(jl,ik)*pmfu(jl,ik) + &
                zerate*puen(jl,jk)-zderate*zuu(jl,ik)+pgf_u)*zmfa
              zvu(jl,jk) = (zvu(jl,ik)*pmfu(jl,ik) + &
                zerate*pven(jl,jk)-zderate*zvu(jl,ik)+pgf_v)*zmfa
            end if
            end if
          end if
        end do
      end do

      if(lmfdd) then
      do jk = 3 , klev
        ik = jk - 1
        do jl = 1,klon
          if ( ldcum(jl) ) then
            if ( jk == idtop(jl) ) then
              zud(jl,jk) = 0.5*(zuu(jl,jk)+puen(jl,ik))
              zvd(jl,jk) = 0.5*(zvu(jl,jk)+pven(jl,ik))
            else if ( jk > idtop(jl) ) then
              zerate = -pmfd(jl,jk) + pmfd(jl,ik) + pmfdde_rate(jl,jk)
              zmfa = 1./min(-cmfcmin,pmfd(jl,jk))
              zud(jl,jk) = (zud(jl,ik)*pmfd(jl,ik) - &
                zerate*puen(jl,ik)+pmfdde_rate(jl,jk)*zud(jl,ik))*zmfa
              zvd(jl,jk) = (zvd(jl,ik)*pmfd(jl,ik) - &
                zerate*pven(jl,ik)+pmfdde_rate(jl,jk)*zvd(jl,ik))*zmfa
            end if
          end if
        end do
      end do
      end if



      zmfs(:) = 1.
      do jk = 2 , klev
        do jl = 1, klon
          if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
            zmfmax = (paph(jl,jk)-paph(jl,jk-1))*zcons
            if ( pmfu(jl,jk) > zmfmax .and. jk >= kctop(jl) ) then
              zmfs(jl) = min(zmfs(jl),zmfmax/pmfu(jl,jk))
            end if
          end if
        end do
      end do
      do jk = 1 , klev
        do jl = 1, klon
          zmfuus(jl,jk) = pmfu(jl,jk)
          zmfdus(jl,jk) = pmfd(jl,jk)
          if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
            zmfuus(jl,jk) = pmfu(jl,jk)*zmfs(jl)
            zmfdus(jl,jk) = pmfd(jl,jk)*zmfs(jl)
          end if
        end do
      end do


     do jk = 1 , klev
        do jl = 1, klon
          ztenu(jl,jk) = pvom(jl,jk)
          ztenv(jl,jk) = pvol(jl,jk)
        end do
      end do

      call cududvn(klon,klev,itopm2,ktype,kcbot,kctop, &
                  ldcum,ztmst,paph,puen,pven,zmfuus,zmfdus,zuu,  &
                  zud,zvu,zvd,pvom,pvol)

      do jl = 1, klon
        zsum12(jl) = 0.
        zsum22(jl) = 0.
      end do
        do jk = 1 , klev
          do jl = 1, klon
            zuv2(jl,jk) = 0.
            if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
              zdz = (paph(jl,jk+1)-paph(jl,jk))
              zduten = pvom(jl,jk) - ztenu(jl,jk)
              zdvten = pvol(jl,jk) - ztenv(jl,jk)
              zuv2(jl,jk) = sqrt(zduten**2+zdvten**2)
              zsum22(jl) = zsum22(jl) + zuv2(jl,jk)*zdz
              zsum12(jl) = zsum12(jl) - &
                (puen(jl,jk)*zduten+pven(jl,jk)*zdvten)*zdz
            end if
          end do
        end do
        do jk = 1 , klev
          do jl = 1, klon
            if ( ldcum(jl) .and. jk>=kctop(jl)-1 ) then
              ztdis = rcpd*zsum12(jl)*zuv2(jl,jk)/max(1.e-15,zsum22(jl))
              ptte(jl,jk) = ptte(jl,jk) + ztdis
            end if
          end do
        end do

      end if




    if ( .not. lmfscv .or. .not. lmfpen ) then
      do jk = 2 , klev
        do jl = 1, klon
          if ( llo2(jl) .and. jk >= kctop(jl)-1 ) then
            ptu(jl,jk) = pten(jl,jk)
            pqu(jl,jk) = pqen(jl,jk)
            plu(jl,jk) = 0.
            pmfude_rate(jl,jk) = 0.
            pmfdde_rate(jl,jk) = 0.
          end if
        end do
      end do
      do jl = 1, klon
        if ( llo2(jl) ) then
          kctop(jl) = klev - 1
          kcbot(jl) = klev - 1
        end if
      end do
    end if

      return
      end subroutine cumastrn



     subroutine cuinin &
     &    (klon,     klev,     klevp1,   klevm1,   pten,&
     &     pqen,     pqsen,    puen,     pven,     pverv,&
     &     pgeo,     paph,     pgeoh,    ptenh,    pqenh,&
     &     pqsenh,   klwmin,   ptu,      pqu,      ptd,&
     &     pqd,      puu,      pvu,      pud,      pvd,&
     &     pmfu,     pmfd,     pmfus,    pmfds,    pmfuq,&
     &     pmfdq,    pdmfup,   pdmfdp,   pdpmel,   plu,&
     &     plude,    klab)
      implicit none
















      integer  klon,klev,klevp1,klevm1
      real     pten(klon,klev),        pqen(klon,klev),&
     &         puen(klon,klev),        pven(klon,klev),&
     &         pqsen(klon,klev),       pverv(klon,klev),&
     &         pgeo(klon,klev),        pgeoh(klon,klevp1),&
     &         paph(klon,klevp1),      ptenh(klon,klev),&
     &         pqenh(klon,klev),       pqsenh(klon,klev)
      real     ptu(klon,klev),         pqu(klon,klev),&
     &         ptd(klon,klev),         pqd(klon,klev),&
     &         puu(klon,klev),         pud(klon,klev),&
     &         pvu(klon,klev),         pvd(klon,klev),&
     &         pmfu(klon,klev),        pmfd(klon,klev),&
     &         pmfus(klon,klev),       pmfds(klon,klev),&
     &         pmfuq(klon,klev),       pmfdq(klon,klev),&
     &         pdmfup(klon,klev),      pdmfdp(klon,klev),&
     &         plu(klon,klev),         plude(klon,klev)
      real     zwmax(klon),            zph(klon), &
     &         pdpmel(klon,klev)
      integer  klab(klon,klev),        klwmin(klon)
      logical  loflag(klon)

      integer  jl,jk
      integer  icall,ik
      real     zzs





      do jk=2,klev
      do jl=1,klon
        ptenh(jl,jk)=(max(cpd*pten(jl,jk-1)+pgeo(jl,jk-1), &
     &             cpd*pten(jl,jk)+pgeo(jl,jk))-pgeoh(jl,jk))*rcpd
        pqenh(jl,jk) = pqen(jl,jk-1)
        pqsenh(jl,jk)= pqsen(jl,jk-1)
        zph(jl)=paph(jl,jk)
        loflag(jl)=.true.
      end do

      if ( jk >= klev-1 .or. jk < 2 ) cycle
      ik=jk
      icall=0
      call cuadjtqn(klon,klev,ik,zph,ptenh,pqsenh,loflag,icall)
      do jl=1,klon
        pqenh(jl,jk)=min(pqen(jl,jk-1),pqsen(jl,jk-1)) &
     &            +(pqsenh(jl,jk)-pqsen(jl,jk-1))
        pqenh(jl,jk)=max(pqenh(jl,jk),0.)
      end do
      end do

      do jl=1,klon
        ptenh(jl,klev)=(cpd*pten(jl,klev)+pgeo(jl,klev)- &
     &                pgeoh(jl,klev))*rcpd
        pqenh(jl,klev)=pqen(jl,klev)
        ptenh(jl,1)=pten(jl,1)
        pqenh(jl,1)=pqen(jl,1)
        klwmin(jl)=klev
        zwmax(jl)=0.
      end do

      do jk=klevm1,2,-1
      do jl=1,klon
        zzs=max(cpd*ptenh(jl,jk)+pgeoh(jl,jk), &
     &        cpd*ptenh(jl,jk+1)+pgeoh(jl,jk+1))
        ptenh(jl,jk)=(zzs-pgeoh(jl,jk))*rcpd
      end do
      end do

      do jk=klev,3,-1
      do jl=1,klon
        if(pverv(jl,jk).lt.zwmax(jl)) then
           zwmax(jl)=pverv(jl,jk)
           klwmin(jl)=jk
        end if
      end do
      end do



      do jk=1,klev
      ik=jk-1
      if(jk.eq.1) ik=1
      do jl=1,klon
      ptu(jl,jk)=ptenh(jl,jk)
      ptd(jl,jk)=ptenh(jl,jk)
      pqu(jl,jk)=pqenh(jl,jk)
      pqd(jl,jk)=pqenh(jl,jk)
      plu(jl,jk)=0.
      puu(jl,jk)=puen(jl,ik)
      pud(jl,jk)=puen(jl,ik)
      pvu(jl,jk)=pven(jl,ik)
      pvd(jl,jk)=pven(jl,ik)
      klab(jl,jk)=0
      end do
      end do
      return
      end subroutine cuinin



      subroutine cutypen & 
     &   (  klon,    klev,     klevp1,   klevm1,     pqen,&
     &     ptenh,   pqenh,     pqsenh,    pgeoh,     paph,& 
     &       hfx,     qfx,       pgeo,    pqsen,      pap,&
     &      pten,    lndj,       cutu,     cuqu,    culab,&
     &     ldcum,   cubot,      cutop,    ktype,    wbase,    culu,   kdpl  )









































      implicit none

      integer  klon, klev, klevp1, klevm1
      real     ptenh(klon,klev),       pqenh(klon,klev),& 
     &         pqsen(klon,klev),       pqsenh(klon,klev),&
     &         pgeoh(klon,klevp1),     paph(klon,klevp1),&
     &         pap(klon,klev),         pqen(klon,klev)
      real     pten(klon,klev)
      real     ptu(klon,klev),pqu(klon,klev),plu(klon,klev)
      real     pgeo(klon,klev)
      integer  klab(klon,klev)
      integer  kctop(klon),kcbot(klon)

      real     qfx(klon),hfx(klon)
      real     zph(klon)
      integer  lndj(klon)
      logical  loflag(klon), deepflag(klon), resetflag(klon)

      real     cutu(klon,klev), cuqu(klon,klev), culu(klon,klev)
      integer  culab(klon,klev)
      real     wbase(klon)
      integer  ktype(klon),cubot(klon),cutop(klon),kdpl(klon)
      logical  ldcum(klon)

      real     zqold(klon)
      real     rho, part1, part2, root, conw, deltt, deltq
      real     eta(klon),dz(klon),coef(klon)
      real     dhen(klon,klev), dh(klon,klev)
      real     plude(klon,klev)
      real     kup(klon,klev)
      real     vptu(klon,klev),vten(klon,klev)
      real     zbuo(klon,klev),abuoy(klon,klev)
  
      real     zz,zdken,zdq
      real     fscale,crirh1,pp
      real     atop1,atop2,abot
      real     tmix,zmix,qmix,pmix
      real     zlglac,dp
      integer  nk,is,ikb,ikt

      real     zqsu,zcor,zdp,zesdp,zalfaw,zfacw,zfaci,zfac,zdsdp,zdqsdt,zdtdp
      real     zpdifftop, zpdiffbot
      integer  zcbase(klon), itoppacel(klon)
      integer  jl,jk,ik,icall,levels
      logical  needreset, lldcum(klon)

      do jl=1,klon
        kcbot(jl)=klev
        kctop(jl)=klev
        kdpl(jl) =klev
        ktype(jl)=0
        wbase(jl)=0.
        ldcum(jl)=.false.
      end do





      do jk=1,klev
      do jl=1,klon
          plu(jl,jk)=culu(jl,jk)  ! parcel liquid water
          ptu(jl,jk)=cutu(jl,jk)  ! parcel temperature
          pqu(jl,jk)=cuqu(jl,jk)  ! parcel specific humidity
          klab(jl,jk)=culab(jl,jk)
           dh(jl,jk)=0.0  ! parcel dry static energy
         dhen(jl,jk)=0.0  ! environment dry static energy
          kup(jl,jk)=0.0  ! updraught kinetic energy for parcel
         vptu(jl,jk)=0.0  ! parcel virtual temperature considering water-loading
         vten(jl,jk)=0.0  ! environment virtual temperature
         zbuo(jl,jk)=0.0  ! parcel buoyancy
         abuoy(jl,jk)=0.0
      end do
      end do

      do jl=1,klon
         zqold(jl) = 0.
         lldcum(jl) = .false.
         loflag(jl) = .true.
      end do

      do jk=klevm1,2,-1

      if(jk .eq. klevm1) then
      do jl=1,klon
        rho=pap(jl,klev)/ &
     &         (rd*(pten(jl,klev)*(1.+vtmpc1*pqen(jl,klev))))
        part1 = 1.5*0.4*pgeo(jl,klev)/ &
     &              (rho*pten(jl,klev))
        part2 = -hfx(jl)*rcpd-vtmpc1*pten(jl,klev)*qfx(jl)
        root  = 0.001-part1*part2
        if(part2 .lt. 0.) then
           conw  = 1.2*(root)**t13
           deltt = max(1.5*hfx(jl)/(rho*cpd*conw),0.)
           deltq = max(1.5*qfx(jl)/(rho*conw),0.)
           kup(jl,klev) = 0.5*(conw**2)
           pqu(jl,klev)= pqenh(jl,klev) + deltq
          dhen(jl,klev)= pgeoh(jl,klev) + ptenh(jl,klev)*cpd
           dh(jl,klev) = dhen(jl,klev)  + deltt*cpd
          ptu(jl,klev) = (dh(jl,klev)-pgeoh(jl,klev))*rcpd 
          vptu(jl,klev)=ptu(jl,klev)*(1.+vtmpc1*pqu(jl,klev))
          vten(jl,klev)=ptenh(jl,klev)*(1.+vtmpc1*pqenh(jl,klev))
          zbuo(jl,klev)=(vptu(jl,klev)-vten(jl,klev))/vten(jl,klev)
          klab(jl,klev) = 1
        else
          loflag(jl) = .false.
        end if
      end do
      end if
 
      is=0
      do jl=1,klon
         if(loflag(jl))then
            is=is+1
         endif
      enddo
      if(is.eq.0) exit

      do jl=1,klon
      if(loflag(jl)) then
        eta(jl) = 0.8/(pgeo(jl,jk)*zrg)+2.e-4
        dz(jl)  = (pgeoh(jl,jk)-pgeoh(jl,jk+1))*zrg
        coef(jl)= 0.5*eta(jl)*dz(jl)
        dhen(jl,jk) = pgeoh(jl,jk) + cpd*ptenh(jl,jk)
        dh(jl,jk) = (coef(jl)*(dhen(jl,jk+1)+dhen(jl,jk))&
     &              +(1.-coef(jl))*dh(jl,jk+1))/(1.+coef(jl))
        pqu(jl,jk) =(coef(jl)*(pqenh(jl,jk+1)+pqenh(jl,jk))&
     &              +(1.-coef(jl))*pqu(jl,jk+1))/(1.+coef(jl))
        ptu(jl,jk) = (dh(jl,jk)-pgeoh(jl,jk))*rcpd
        zqold(jl) = pqu(jl,jk)
        zph(jl)=paph(jl,jk)
      end if
      end do

      ik=jk
      icall=1
      call cuadjtqn(klon,klev,ik,zph,ptu,pqu,loflag,icall)
      do jl=1,klon
        if( loflag(jl) ) then
          zdq = max((zqold(jl) - pqu(jl,jk)),0.)
          plu(jl,jk) = plu(jl,jk+1) + zdq
          zlglac=zdq*((1.-foealfa(ptu(jl,jk))) - &
                      (1.-foealfa(ptu(jl,jk+1))))
          plu(jl,jk) = min(plu(jl,jk),5.e-3)
          dh(jl,jk) =  pgeoh(jl,jk) + cpd*(ptu(jl,jk)+ralfdcp*zlglac)

          vptu(jl,jk) = ptu(jl,jk)*(1.+vtmpc1*pqu(jl,jk)-plu(jl,jk))+&
                        ralfdcp*zlglac
          vten(jl,jk) = ptenh(jl,jk)*(1.+vtmpc1*pqenh(jl,jk))
          zbuo(jl,jk) = (vptu(jl,jk) - vten(jl,jk))/vten(jl,jk)
          abuoy(jl,jk)=(zbuo(jl,jk)+zbuo(jl,jk+1))*0.5*g
          atop1 = 1.0 - 2.*coef(jl)
          atop2 = 2.0*dz(jl)*abuoy(jl,jk)
          abot =  1.0 + 2.*coef(jl)
          kup(jl,jk)  = (atop1*kup(jl,jk+1) + atop2) / abot

         if ( plu(jl,jk) > 0. .and. klab(jl,jk+1) == 1 ) then
              ik = jk + 1
              zqsu = foeewm(ptu(jl,ik))/paph(jl,ik)
              zqsu = min(0.5,zqsu)
              zcor = 1./(1.-vtmpc1*zqsu)
              zqsu = zqsu*zcor
              zdq = min(0.,pqu(jl,ik)-zqsu)
              zalfaw = foealfa(ptu(jl,ik))
              zfacw = c5les/((ptu(jl,ik)-c4les)**2)
              zfaci = c5ies/((ptu(jl,ik)-c4ies)**2)
              zfac = zalfaw*zfacw + (1.-zalfaw)*zfaci
              zesdp = foeewm(ptu(jl,ik))/paph(jl,ik)
              zcor = 1./(1.-vtmpc1*zesdp)
              zdqsdt = zfac*zcor*zqsu
              zdtdp = rd*ptu(jl,ik)/(cpd*paph(jl,ik))
              zdp = zdq/(zdqsdt*zdtdp)
              zcbase(jl) = paph(jl,ik) + zdp

              zpdifftop = zcbase(jl) - paph(jl,jk)
              zpdiffbot = paph(jl,jk+1) - zcbase(jl)
              if ( zpdifftop > zpdiffbot .and. kup(jl,jk+1) > 0. ) then
                ikb = min(klev-1,jk+1)
                klab(jl,ikb) = 2
                klab(jl,jk) = 2
                kcbot(jl) = ikb
                plu(jl,jk+1) = 1.0e-8
              else if ( zpdifftop <= zpdiffbot .and.kup(jl,jk) > 0. ) then
                klab(jl,jk) = 2
                kcbot(jl) = jk
              end if
          end if

          if(kup(jl,jk) .lt. 0.)then
            loflag(jl) = .false.
            if(plu(jl,jk+1) .gt. 0.) then
              kctop(jl) = jk
              lldcum(jl) = .true.
            else
              lldcum(jl) = .false.
            end if
          else 
            if(plu(jl,jk) .gt. 0.)then
              klab(jl,jk)=2
            else
              klab(jl,jk)=1
            end if
          end if
        end if
      end do

      end do ! end all the levels

      do jl=1,klon
        ikb = kcbot(jl)
        ikt = kctop(jl)
        if(paph(jl,ikb) - paph(jl,ikt) > zdnoprc) lldcum(jl) = .false.
        if(lldcum(jl)) then
            ktype(jl) = 2
            ldcum(jl) = .true.
            wbase(jl) = sqrt(max(2.*kup(jl,ikb),0.))
            cubot(jl) = ikb
            cutop(jl) = ikt
            kdpl(jl)  = klev
        else
            cutop(jl) = -1
            cubot(jl) = -1
            kdpl(jl)  = klev - 1
            ldcum(jl) = .false.
            wbase(jl) = 0.
        end if
      end do

      do jk=klev,1,-1
       do jl=1,klon
           ikt = kctop(jl)
           if(jk .ge. ikt)then
             culab(jl,jk) = klab(jl,jk)
             cutu(jl,jk)  = ptu(jl,jk)
             cuqu(jl,jk)  = pqu(jl,jk)
             culu(jl,jk)  = plu(jl,jk)
           end if
       end do
      end do
      







      deltt = 0.2
      deltq = 1.0e-4
      do jl=1,klon
        deepflag(jl) = .false.
      end do

      do jk=klev,1,-1
       do jl=1,klon
         if((paph(jl,klev+1)-paph(jl,jk)) .lt. 350.e2) itoppacel(jl) = jk
       end do
      end do

      do levels=klevm1-1,klev/2+1,-1 ! loop starts
        do jk=1,klev
          do jl=1,klon
             plu(jl,jk)=0.0  ! parcel liquid water
             ptu(jl,jk)=0.0  ! parcel temperature
             pqu(jl,jk)=0.0  ! parcel specific humidity
             dh(jl,jk)=0.0   ! parcel dry static energy
             dhen(jl,jk)=0.0  ! environment dry static energy
             kup(jl,jk)=0.0   ! updraught kinetic energy for parcel
             vptu(jl,jk)=0.0  ! parcel virtual temperature consideringwater-loading
             vten(jl,jk)=0.0  ! environment virtual temperature
             abuoy(jl,jk)=0.0
             zbuo(jl,jk)=0.0
             klab(jl,jk)=0
          end do
        end do

        do jl=1,klon
           kcbot(jl)    =  levels
           kctop(jl)    =  levels
           zqold(jl)    = 0.
           lldcum(jl)   = .false.
           resetflag(jl)= .false.
           loflag(jl)   = (.not. deepflag(jl)) .and. (levels.ge.itoppacel(jl))
        end do

      do jk=levels,2,-1
        is=0
        do jl=1,klon
         if(loflag(jl))then
            is=is+1
         endif
        enddo
        if(is.eq.0) exit

        if(jk .eq. levels) then
          do jl=1,klon
          if(loflag(jl)) then
            if((paph(jl,klev+1)-paph(jl,jk)) < 60.e2) then
              tmix=0.
              qmix=0.
              zmix=0.
              pmix=0.
              do nk=jk+2,jk,-1
              if(pmix < 50.e2) then
                dp = paph(jl,nk) - paph(jl,nk-1)
                tmix=tmix+dp*ptenh(jl,nk)
                qmix=qmix+dp*pqenh(jl,nk)
                zmix=zmix+dp*pgeoh(jl,nk)
                pmix=pmix+dp
              end if
              end do
              tmix=tmix/pmix
              qmix=qmix/pmix
              zmix=zmix/pmix
            else
              tmix=ptenh(jl,jk+1)
              qmix=pqenh(jl,jk+1)
              zmix=pgeoh(jl,jk+1)
            end if

            pqu(jl,jk+1) = qmix + deltq
            dhen(jl,jk+1)= zmix + tmix*cpd
            dh(jl,jk+1)  = dhen(jl,jk+1) + deltt*cpd
            ptu(jl,jk+1) = (dh(jl,jk+1)-pgeoh(jl,jk+1))*rcpd
            kup(jl,jk+1) = 0.5
            klab(jl,jk+1)= 1
            vptu(jl,jk+1)=ptu(jl,jk+1)*(1.+vtmpc1*pqu(jl,jk+1))
            vten(jl,jk+1)=ptenh(jl,jk+1)*(1.+vtmpc1*pqenh(jl,jk+1))
            zbuo(jl,jk+1)=(vptu(jl,jk+1)-vten(jl,jk+1))/vten(jl,jk+1)
          end if
          end do
        end if

        do jl=1,klon
           if(loflag(jl)) then

             fscale = min(1.,(pqsen(jl,jk)/pqsen(jl,levels))**3)
             eta(jl) = 1.75e-3*fscale
             dz(jl)  = (pgeoh(jl,jk)-pgeoh(jl,jk+1))*zrg
             coef(jl)= 0.5*eta(jl)*dz(jl)
             dhen(jl,jk) = pgeoh(jl,jk) + cpd*ptenh(jl,jk)
             dh(jl,jk) = (coef(jl)*(dhen(jl,jk+1)+dhen(jl,jk))&
     &              +(1.-coef(jl))*dh(jl,jk+1))/(1.+coef(jl))
             pqu(jl,jk) =(coef(jl)*(pqenh(jl,jk+1)+pqenh(jl,jk))&
     &              +(1.-coef(jl))*pqu(jl,jk+1))/(1.+coef(jl))
             ptu(jl,jk) = (dh(jl,jk)-pgeoh(jl,jk))*rcpd
             zqold(jl) = pqu(jl,jk)
             zph(jl)=paph(jl,jk)
           end if
        end do

        ik=jk
        icall=1
        call cuadjtqn(klon,klev,ik,zph,ptu,pqu,loflag,icall)
 
      do jl=1,klon
        if( loflag(jl) ) then
          zdq = max((zqold(jl) - pqu(jl,jk)),0.)
          plu(jl,jk) = plu(jl,jk+1) + zdq
          zlglac=zdq*((1.-foealfa(ptu(jl,jk))) - &
                      (1.-foealfa(ptu(jl,jk+1))))
          plu(jl,jk) = 0.5*plu(jl,jk)
          dh(jl,jk) =  pgeoh(jl,jk) + cpd*(ptu(jl,jk)+ralfdcp*zlglac)

          vptu(jl,jk) = ptu(jl,jk)*(1.+vtmpc1*pqu(jl,jk)-plu(jl,jk))+&
                        ralfdcp*zlglac
          vten(jl,jk) = ptenh(jl,jk)*(1.+vtmpc1*pqenh(jl,jk))
          zbuo(jl,jk) = (vptu(jl,jk) - vten(jl,jk))/vten(jl,jk)
          abuoy(jl,jk)=(zbuo(jl,jk)+zbuo(jl,jk+1))*0.5*g
          atop1 = 1.0 - 2.*coef(jl)
          atop2 = 2.0*dz(jl)*abuoy(jl,jk)
          abot =  1.0 + 2.*coef(jl)
          kup(jl,jk)  = (atop1*kup(jl,jk+1) + atop2) / abot

          if ( plu(jl,jk) > 0. .and. klab(jl,jk+1) == 1 ) then
              ik = jk + 1
              zqsu = foeewm(ptu(jl,ik))/paph(jl,ik)
              zqsu = min(0.5,zqsu)
              zcor = 1./(1.-vtmpc1*zqsu)
              zqsu = zqsu*zcor
              zdq = min(0.,pqu(jl,ik)-zqsu)
              zalfaw = foealfa(ptu(jl,ik))
              zfacw = c5les/((ptu(jl,ik)-c4les)**2)
              zfaci = c5ies/((ptu(jl,ik)-c4ies)**2)
              zfac = zalfaw*zfacw + (1.-zalfaw)*zfaci
              zesdp = foeewm(ptu(jl,ik))/paph(jl,ik)
              zcor = 1./(1.-vtmpc1*zesdp)
              zdqsdt = zfac*zcor*zqsu
              zdtdp = rd*ptu(jl,ik)/(cpd*paph(jl,ik))
              zdp = zdq/(zdqsdt*zdtdp)
              zcbase(jl) = paph(jl,ik) + zdp

              zpdifftop = zcbase(jl) - paph(jl,jk)
              zpdiffbot = paph(jl,jk+1) - zcbase(jl)
              if ( zpdifftop > zpdiffbot .and. kup(jl,jk+1) > 0. ) then
                ikb = min(klev-1,jk+1)
                klab(jl,ikb) = 2
                klab(jl,jk) = 2
                kcbot(jl) = ikb
                plu(jl,jk+1) = 1.0e-8
              else if ( zpdifftop <= zpdiffbot .and.kup(jl,jk) > 0. ) then
                klab(jl,jk) = 2
                kcbot(jl) = jk
              end if
          end if

          if(kup(jl,jk) .lt. 0.)then
            loflag(jl) = .false.
            if(plu(jl,jk+1) .gt. 0.) then
              kctop(jl) = jk
              lldcum(jl) = .true.
            else
              lldcum(jl) = .false.
            end if
          else 
            if(plu(jl,jk) .gt. 0.)then
              klab(jl,jk)=2
            else
              klab(jl,jk)=1
            end if
          end if
        end if
      end do

      end do ! end all the levels

      needreset = .false.
      do jl=1,klon
        ikb = kcbot(jl)
        ikt = kctop(jl)
        if(paph(jl,ikb) - paph(jl,ikt) < zdnoprc) lldcum(jl) = .false.
        if(lldcum(jl)) then      
         ktype(jl)    = 1
         ldcum(jl)    = .true.
         deepflag(jl) = .true.
         wbase(jl)    = sqrt(max(2.*kup(jl,ikb),0.))
         cubot(jl)    = ikb
         cutop(jl)    = ikt
         kdpl(jl)     = levels+1
         needreset    = .true.
         resetflag(jl)= .true.
        end if
      end do

      if(needreset) then
      do jk=klev,1,-1
        do jl=1,klon
          if(resetflag(jl)) then
            ikt = kctop(jl)
            ikb = kdpl(jl)
            if(jk .le. ikb .and. jk .ge. ikt )then
             culab(jl,jk) = klab(jl,jk)
             cutu(jl,jk)  = ptu(jl,jk)
             cuqu(jl,jk)  = pqu(jl,jk)
             culu(jl,jk)  = plu(jl,jk)
            else
             culab(jl,jk) = 1
             cutu(jl,jk)  = ptenh(jl,jk)
             cuqu(jl,jk)  = pqenh(jl,jk)
             culu(jl,jk)  = 0.
            end if
            if ( jk .lt. ikt ) culab(jl,jk) = 0
          end if
        end do
      end do
      end if

      end do ! end all cycles

      return
      end subroutine cutypen



      subroutine cuascn &
     &    (klon,     klev,     klevp1,   klevm1,   ptenh,&
     &     pqenh,    puen,     pven,     pten,     pqen,&
     &     pqsen,    pgeo,     pgeoh,    pap,      paph,&
     &     pqte,     pverv,    klwmin,   ldcum,    phcbase,&
     &     ktype,    klab,     ptu,      pqu,      plu,&
     &     puu,      pvu,      pmfu,     pmfub,    &
     &     pmfus,    pmfuq,    pmful,    plude,    pdmfup,&
     &     kcbot,    kctop,    kctop0,   kcum,     ztmst,&
     &     pqsenh,   plglac,   lndj,     wup,      wbase,   kdpl, pmfude_rate)
      implicit none






































































      integer  klev,klon,klevp1,klevm1
      real     ptenh(klon,klev),       pqenh(klon,klev), &
     &         puen(klon,klev),        pven(klon,klev),&
     &         pten(klon,klev),        pqen(klon,klev),&
     &         pgeo(klon,klev),        pgeoh(klon,klevp1),&
     &         pap(klon,klev),         paph(klon,klevp1),&
     &         pqsen(klon,klev),       pqte(klon,klev),&
     &         pverv(klon,klev),       pqsenh(klon,klev)
      real     ptu(klon,klev),         pqu(klon,klev),&
     &         puu(klon,klev),         pvu(klon,klev),&
     &         pmfu(klon,klev),        zph(klon),&
     &         pmfub(klon),            &
     &         pmfus(klon,klev),       pmfuq(klon,klev),&
     &         plu(klon,klev),         plude(klon,klev),&
     &         pmful(klon,klev),       pdmfup(klon,klev)
      real     zdmfen(klon),           zdmfde(klon),&
     &         zmfuu(klon),            zmfuv(klon),&
     &         zpbase(klon),           zqold(klon)              
      real     phcbase(klon),          zluold(klon)
      real     zprecip(klon),          zlrain(klon,klev)
      real     zbuo(klon,klev),        kup(klon,klev)
      real     wup(klon)
      real     wbase(klon),            zodetr(klon,klev)
      real     plglac(klon,klev)

      real     eta(klon),dz(klon)

      integer  klwmin(klon),           ktype(klon),&
     &         klab(klon,klev),        kcbot(klon),&
     &         kctop(klon),            kctop0(klon)
      integer  lndj(klon)
      logical  ldcum(klon),            loflag(klon)
      logical  llo2,llo3,              llo1(klon)

      integer  kdpl(klon)
      real     zoentr(klon),           zdpmean(klon)
      real     pdmfen(klon,klev),      pmfude_rate(klon,klev)

      integer  jl,jk
      integer  ikb,icum,itopm2,ik,icall,is,kcum,jlm,jll
      integer  jlx(klon)
      real     ztmst,zcons2,zfacbuo,zprcdgw,z_cwdrag,z_cldmax,z_cwifrac,z_cprc2
      real     zmftest,zmfmax,zqeen,zseen,zscde,zqude
      real     zmfusk,zmfuqk,zmfulk
      real     zbc,zbe,zkedke,zmfun,zwu,zprcon,zdt,zcbf,zzco
      real     zlcrit,zdfi,zc,zd,zint,zlnew,zvw,zvi,zalfaw,zrold
      real     zrnew,zz,zdmfeu,zdmfdu,dp
      real     zfac,zbuoc,zdkbuo,zdken,zvv,zarg,zchange,zxe,zxs,zdshrd
      real     atop1,atop2,abot



      zcons2=3./(g*ztmst)
      zfacbuo = 0.5/(1.+0.5)
      zprcdgw = cprcon*zrg
      z_cldmax = 5.e-3
      z_cwifrac = 0.5
      z_cprc2 = 0.5
      z_cwdrag = (3.0/8.0)*0.506/0.2



      llo3 = .false.
      do jl=1,klon
        zluold(jl)=0.
        wup(jl)=0.
        zdpmean(jl)=0.
        zoentr(jl)=0.
        if(.not.ldcum(jl)) then
          ktype(jl)=0
          kcbot(jl) = -1
          pmfub(jl) = 0.
          pqu(jl,klev) = 0.
        end if
      end do


      do jk=1,klev
      do jl=1,klon
          if(jk.ne.kcbot(jl)) plu(jl,jk)=0.
          pmfu(jl,jk)=0.
          pmfus(jl,jk)=0.
          pmfuq(jl,jk)=0.
          pmful(jl,jk)=0.
          plude(jl,jk)=0.
          plglac(jl,jk)=0.
          pdmfup(jl,jk)=0.
          zlrain(jl,jk)=0.
          zbuo(jl,jk)=0.
          kup(jl,jk)=0.
          pdmfen(jl,jk) = 0.
          pmfude_rate(jl,jk) = 0.
          if(.not.ldcum(jl).or.ktype(jl).eq.3) klab(jl,jk)=0
          if(.not.ldcum(jl).and.paph(jl,jk).lt.4.e4) kctop0(jl)=jk
      end do
      end do

      do jl = 1,klon
        if ( ktype(jl) == 3 ) ldcum(jl) = .false.
      end do



      do jl=1,klon
        kctop(jl)=kcbot(jl)
        if(ldcum(jl)) then
          ikb = kcbot(jl)
          kup(jl,ikb) = 0.5*wbase(jl)**2
          pmfu(jl,ikb) = pmfub(jl)
          pmfus(jl,ikb) = pmfub(jl)*(cpd*ptu(jl,ikb)+pgeoh(jl,ikb))
          pmfuq(jl,ikb) = pmfub(jl)*pqu(jl,ikb)
          pmful(jl,ikb) = pmfub(jl)*plu(jl,ikb)
        end if
      end do








      do jk=klevm1,3,-1



      ik=jk
      call cubasmcn&
     &    (klon,     klev,     klevm1,   ik,      pten,&
     &     pqen,     pqsen,    puen,     pven,    pverv,&
     &     pgeo,     pgeoh,    ldcum,   ktype,   klab,  zlrain,&
     &     pmfu,     pmfub,    kcbot,   ptu,&
     &     pqu,      plu,      puu,     pvu,      pmfus,&
     &     pmfuq,    pmful,    pdmfup)
      is = 0
      jlm = 0
      do jl = 1,klon
        loflag(jl) = .false.
        zprecip(jl) = 0.
        llo1(jl) = .false.
        is = is + klab(jl,jk+1)
        if ( klab(jl,jk+1) == 0 ) klab(jl,jk) = 0
        if ( (ldcum(jl) .and. klab(jl,jk+1) == 2) .or.   &
             (ktype(jl) == 3 .and. klab(jl,jk+1) == 1) ) then
          loflag(jl) = .true.
          jlm = jlm + 1
          jlx(jlm) = jl
        end if
        zph(jl) = paph(jl,jk)
        if ( ktype(jl) == 3 .and. jk == kcbot(jl) ) then
          zmfmax = (paph(jl,jk)-paph(jl,jk-1))*zcons2
          if ( pmfub(jl) > zmfmax ) then
            zfac = zmfmax/pmfub(jl)
            pmfu(jl,jk+1) = pmfu(jl,jk+1)*zfac
            pmfus(jl,jk+1) = pmfus(jl,jk+1)*zfac
            pmfuq(jl,jk+1) = pmfuq(jl,jk+1)*zfac
            pmfub(jl) = zmfmax
          end if
          pmfub(jl)=min(pmfub(jl),zmfmax)
        end if
      end do

      if(is.gt.0) llo3 = .true.



      ik=jk
      call cuentrn(klon,klev,ik,kcbot,ldcum,llo3, &
                  pgeoh,pmfu,zdmfen,zdmfde)


      if(llo3) then


        do jl = 1,klon
          zqold(jl) = 0.
        end do
        do jll = 1 , jlm
          jl = jlx(jll)
          zdmfde(jl) = min(zdmfde(jl),0.75*pmfu(jl,jk+1))
          if ( jk == kcbot(jl) ) then
            zoentr(jl) = -1.75e-3*(min(1.,pqen(jl,jk)/pqsen(jl,jk)) - &
                         1.)*(pgeoh(jl,jk)-pgeoh(jl,jk+1))*zrg
            zoentr(jl) = min(0.4,zoentr(jl))*pmfu(jl,jk+1)
          end if
          if ( jk < kcbot(jl) ) then
            zmfmax = (paph(jl,jk)-paph(jl,jk-1))*zcons2
            zxs = max(pmfu(jl,jk+1)-zmfmax,0.)
            wup(jl) = wup(jl) + kup(jl,jk+1)*(pap(jl,jk+1)-pap(jl,jk))
            zdpmean(jl) = zdpmean(jl) + pap(jl,jk+1) - pap(jl,jk)
            zdmfen(jl) = zoentr(jl)
            if ( ktype(jl) >= 2 ) then
              zdmfen(jl) = 2.0*zdmfen(jl)
              zdmfde(jl) = zdmfen(jl)
            end if
            zdmfde(jl) = zdmfde(jl) * &
                         (1.6-min(1.,pqen(jl,jk)/pqsen(jl,jk)))
            zmftest = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
            zchange = max(zmftest-zmfmax,0.)
            zxe = max(zchange-zxs,0.)
            zdmfen(jl) = zdmfen(jl) - zxe
            zchange = zchange - zxe
            zdmfde(jl) = zdmfde(jl) + zchange
          end if
          pdmfen(jl,jk) = zdmfen(jl) - zdmfde(jl)
          pmfu(jl,jk) = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
          zqeen = pqenh(jl,jk+1)*zdmfen(jl)
          zseen = (cpd*ptenh(jl,jk+1)+pgeoh(jl,jk+1))*zdmfen(jl)
          zscde = (cpd*ptu(jl,jk+1)+pgeoh(jl,jk+1))*zdmfde(jl)
          zqude = pqu(jl,jk+1)*zdmfde(jl)
          plude(jl,jk) = plu(jl,jk+1)*zdmfde(jl)
          zmfusk = pmfus(jl,jk+1) + zseen - zscde
          zmfuqk = pmfuq(jl,jk+1) + zqeen - zqude
          zmfulk = pmful(jl,jk+1) - plude(jl,jk)
          plu(jl,jk) = zmfulk*(1./max(cmfcmin,pmfu(jl,jk)))
          pqu(jl,jk) = zmfuqk*(1./max(cmfcmin,pmfu(jl,jk)))
          ptu(jl,jk) = (zmfusk * &
            (1./max(cmfcmin,pmfu(jl,jk)))-pgeoh(jl,jk))*rcpd
          ptu(jl,jk) = max(100.,ptu(jl,jk))
          ptu(jl,jk) = min(400.,ptu(jl,jk))
          zqold(jl) = pqu(jl,jk)
          zlrain(jl,jk) = zlrain(jl,jk+1)*(pmfu(jl,jk+1)-zdmfde(jl)) * &
                          (1./max(cmfcmin,pmfu(jl,jk)))
          zluold(jl) = plu(jl,jk)
        end do

        do jl = 1,klon
          if ( jk > kdpl(jl) ) then
            ptu(jl,jk) = ptenh(jl,jk)
            pqu(jl,jk) = pqenh(jl,jk)
            plu(jl,jk) = 0.
            zluold(jl) = plu(jl,jk)
          end if
        end do



      ik=jk
      icall=1

      if ( jlm > 0 ) then
        call cuadjtqn(klon,klev,ik,zph,ptu,pqu,loflag,icall)
      end if

        do jll = 1 , jlm
          jl = jlx(jll)
          if ( pqu(jl,jk) /= zqold(jl) ) then
            plglac(jl,jk) = plu(jl,jk) * &
                           ((1.-foealfa(ptu(jl,jk)))- &
                            (1.-foealfa(ptu(jl,jk+1))))
            ptu(jl,jk) = ptu(jl,jk) + ralfdcp*plglac(jl,jk)
          end if
        end do
        do jll = 1 , jlm
          jl = jlx(jll)
          if ( pqu(jl,jk) /= zqold(jl) ) then
            klab(jl,jk) = 2
            plu(jl,jk) = plu(jl,jk) + zqold(jl) - pqu(jl,jk)
            zbc = ptu(jl,jk)*(1.+vtmpc1*pqu(jl,jk)-plu(jl,jk+1) - &
              zlrain(jl,jk+1))
            zbe = ptenh(jl,jk)*(1.+vtmpc1*pqenh(jl,jk))
            zbuo(jl,jk) = zbc - zbe

            if ( ktype(jl) == 3 .and. klab(jl,jk+1) == 1 ) then
              if ( zbuo(jl,jk) > -0.5 ) then
                ldcum(jl) = .true.
                kctop(jl) = jk
                kup(jl,jk) = 0.5
              else
                klab(jl,jk) = 0
                pmfu(jl,jk) = 0.
                plude(jl,jk) = 0.
                plu(jl,jk) = 0.
              end if
            end if
            if ( klab(jl,jk+1) == 2 ) then
              if ( zbuo(jl,jk) < 0. ) then
                ptenh(jl,jk) = 0.5*(pten(jl,jk)+pten(jl,jk-1))
                pqenh(jl,jk) = 0.5*(pqen(jl,jk)+pqen(jl,jk-1))
                zbuo(jl,jk) = zbc - ptenh(jl,jk)*(1.+vtmpc1*pqenh(jl,jk))
              end if
              zbuoc = (zbuo(jl,jk) / &
                (ptenh(jl,jk)*(1.+vtmpc1*pqenh(jl,jk)))+zbuo(jl,jk+1) / &
                (ptenh(jl,jk+1)*(1.+vtmpc1*pqenh(jl,jk+1))))*0.5
              zdkbuo = (pgeoh(jl,jk)-pgeoh(jl,jk+1))*zfacbuo*zbuoc

              if ( zdmfen(jl) > 0. ) then
                zdken = min(1.,(1.+z_cwdrag)*zdmfen(jl) / &
                        max(cmfcmin,pmfu(jl,jk+1)))
              else
                zdken = min(1.,(1.+z_cwdrag)*zdmfde(jl) / &
                        max(cmfcmin,pmfu(jl,jk+1)))
              end if
              kup(jl,jk) = (kup(jl,jk+1)*(1.-zdken)+zdkbuo) / &
                           (1.+zdken)
              if ( zbuo(jl,jk) < 0. ) then
                zkedke = kup(jl,jk)/max(1.e-10,kup(jl,jk+1))
                zkedke = max(0.,min(1.,zkedke))
                zmfun = sqrt(zkedke)*pmfu(jl,jk+1)
                zdmfde(jl) = max(zdmfde(jl),pmfu(jl,jk+1)-zmfun)
                plude(jl,jk) = plu(jl,jk+1)*zdmfde(jl)
                pmfu(jl,jk) = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
              end if
              if ( zbuo(jl,jk) > -0.2  ) then
                ikb = kcbot(jl)
                zoentr(jl) = 1.75e-3*(0.3-(min(1.,pqen(jl,jk-1) /    &
                  pqsen(jl,jk-1))-1.))*(pgeoh(jl,jk-1)-pgeoh(jl,jk)) * &
                  zrg*min(1.,pqsen(jl,jk)/pqsen(jl,ikb))**3
                zoentr(jl) = min(0.4,zoentr(jl))*pmfu(jl,jk)
              else
                zoentr(jl) = 0.
              end if

              if ( jk > kdpl(jl) ) then
                pmfu(jl,jk) = pmfu(jl,jk+1)
                kup(jl,jk) = 0.5
              end if
              if ( kup(jl,jk) > 0. .and. pmfu(jl,jk) > 0. ) then
                kctop(jl) = jk
                llo1(jl) = .true.
              else
                klab(jl,jk) = 0
                pmfu(jl,jk) = 0.
                kup(jl,jk) = 0.
                zdmfde(jl) = pmfu(jl,jk+1)
                plude(jl,jk) = plu(jl,jk+1)*zdmfde(jl)
              end if

              if ( pmfu(jl,jk+1) > 0. ) pmfude_rate(jl,jk) = zdmfde(jl)
            end if
          else if ( ktype(jl) == 2 .and. pqu(jl,jk) == zqold(jl) ) then
            klab(jl,jk) = 0
            pmfu(jl,jk) = 0.
            kup(jl,jk) = 0.
            zdmfde(jl) = pmfu(jl,jk+1)
            plude(jl,jk) = plu(jl,jk+1)*zdmfde(jl)
            pmfude_rate(jl,jk) = zdmfde(jl)
          end if
        end do

        do jl = 1,klon
          if ( llo1(jl) ) then



            if ( lndj(jl).eq.1 ) then
              zdshrd = 5.e-4
            else
              zdshrd = 3.e-4
            end if
            ikb=kcbot(jl)
            if ( plu(jl,jk) > zdshrd )then
              zwu = min(15.0,sqrt(2.*max(0.1,kup(jl,jk+1))))
              zprcon = zprcdgw/(0.75*zwu)

              zdt = min(rtber-rtice,max(rtber-ptu(jl,jk),0.))
              zcbf = 1. + z_cprc2*sqrt(zdt)
              zzco = zprcon*zcbf
              zlcrit = zdshrd/zcbf
              zdfi = pgeoh(jl,jk) - pgeoh(jl,jk+1)
              zc = (plu(jl,jk)-zluold(jl))
              zarg = (plu(jl,jk)/zlcrit)**2
              if ( zarg < 25.0 ) then
                zd = zzco*(1.-exp(-zarg))*zdfi
              else
                zd = zzco*zdfi
              end if
              zint = exp(-zd)
              zlnew = zluold(jl)*zint + zc/zd*(1.-zint)
              zlnew = max(0.,min(plu(jl,jk),zlnew))
              zlnew = min(z_cldmax,zlnew)
              zprecip(jl) = max(0.,zluold(jl)+zc-zlnew)
              pdmfup(jl,jk) = zprecip(jl)*pmfu(jl,jk)
              zlrain(jl,jk) = zlrain(jl,jk) + zprecip(jl)
              plu(jl,jk) = zlnew
            end if
          end if
        end do
        do jl = 1, klon
          if ( llo1(jl) ) then
            if ( zlrain(jl,jk) > 0. ) then
              zvw = 21.18*zlrain(jl,jk)**0.2
              zvi = z_cwifrac*zvw
              zalfaw = foealfa(ptu(jl,jk))
              zvv = zalfaw*zvw + (1.-zalfaw)*zvi
              zrold = zlrain(jl,jk) - zprecip(jl)
              zc = zprecip(jl)
              zwu = min(15.0,sqrt(2.*max(0.1,kup(jl,jk))))
              zd = zvv/zwu
              zint = exp(-zd)
              zrnew = zrold*zint + zc/zd*(1.-zint)
              zrnew = max(0.,min(zlrain(jl,jk),zrnew))
              zlrain(jl,jk) = zrnew
            end if
          end if
        end do
        do jll = 1 , jlm
          jl = jlx(jll)
          pmful(jl,jk) = plu(jl,jk)*pmfu(jl,jk)
          pmfus(jl,jk) = (cpd*ptu(jl,jk)+pgeoh(jl,jk))*pmfu(jl,jk)
          pmfuq(jl,jk) = pqu(jl,jk)*pmfu(jl,jk)
        end do
      end if
    end do



      do jl = 1,klon
       if ( kctop(jl) == -1 ) ldcum(jl) = .false.
        kcbot(jl) = max(kcbot(jl),kctop(jl))
        if ( ldcum(jl) ) then
          wup(jl) = max(1.e-2,wup(jl)/max(1.,zdpmean(jl)))
          wup(jl) = sqrt(2.*wup(jl))
        end if
      end do

      return
      end subroutine cuascn



      subroutine cudlfsn   &                                                      
     &    (klon,     klev,  &                       
     &     kcbot,    kctop,    lndj,   ldcum,   &                               
     &     ptenh,    pqenh,    puen,     pven,     &                              
     &     pten,     pqsen,    pgeo,                &                             
     &     pgeoh,    paph,     ptu,      pqu,      plu,&                          
     &     puu,      pvu,      pmfub,    prfl,          &                         
     &     ptd,      pqd,      pud,      pvd,            &                        
     &     pmfd,     pmfds,    pmfdq,    pdmfdp,          &                       
     &     kdtop,    lddraf)                                                     
                                                                                 


                                                                                 

                                                                                 




                                                                                 







                                                                                 

                                                                                 


                                                                                 



                                                                                 




                                                                                 

                                                                                 


                                                                                 

                                                                                 















                                                                                 

                                                                                 

                                                                                 

                                                                                 








                                                                                 

                                                                                 

                                                                                 

                                                                                 

                                                                                 




      implicit none                                                       
                   
      integer  klev,klon                                                            
      real     ptenh(klon,klev),       pqenh(klon,klev), &                        
     &         puen(klon,klev),        pven(klon,klev),  &                        
     &         pten(klon,klev),        pqsen(klon,klev),  &                       
     &         pgeo(klon,klev),                       &                           
     &         pgeoh(klon,klev+1),     paph(klon,klev+1),&                        
     &         ptu(klon,klev),         pqu(klon,klev),   &                        
     &         puu(klon,klev),         pvu(klon,klev),   &                        
     &         plu(klon,klev),                          &                         
     &         pmfub(klon),            prfl(klon)                                
                                                                                 
      real     ptd(klon,klev),         pqd(klon,klev),   &                        
     &         pud(klon,klev),         pvd(klon,klev),    &                       
     &         pmfd(klon,klev),        pmfds(klon,klev),  &                       
     &         pmfdq(klon,klev),       pdmfdp(klon,klev)                         
      integer  kcbot(klon),            kctop(klon),       &                       
     &         kdtop(klon),            ikhsmin(klon)                             
      logical  ldcum(klon),                              &
     &         lddraf(klon)   
      integer  lndj(klon)                                                   
                                                                                 
      real     ztenwb(klon,klev),      zqenwb(klon,klev), &                       
     &         zcond(klon),            zph(klon),         &                       
     &         zhsmin(klon)                                                      
      logical  llo2(klon)                                                        

      integer  jl,jk
      integer  is,ik,icall,ike
      real     zhsk,zttest,zqtest,zbuo,zmftop
                                                                                   

                                                                                 


      do jl=1,klon                                                          
        lddraf(jl)=.false.                                                       
        kdtop(jl)=klev+1      
        ikhsmin(jl)=klev+1
        zhsmin(jl)=1.e8                                                  
      enddo                                                                      

                                                                                 




                                                                                 

                                                                                 






                                                                                 





      do jk=3,klev-2
         do jl=1,klon
           zhsk=cpd*pten(jl,jk)+pgeo(jl,jk) +  &
     &         foelhm(pten(jl,jk))*pqsen(jl,jk)
           if(zhsk .lt. zhsmin(jl)) then
              zhsmin(jl) = zhsk
              ikhsmin(jl)= jk
           end if
         end do
      end do


      ike=klev-3                                                                 
      do jk=3,ike                                                                
                                                                                 



        is=0                                                                     
        do jl=1,klon                                                        
          ztenwb(jl,jk)=ptenh(jl,jk)                                             
          zqenwb(jl,jk)=pqenh(jl,jk)                                             
          zph(jl)=paph(jl,jk)                                                    
          llo2(jl)=ldcum(jl).and.prfl(jl).gt.0..and..not.lddraf(jl).and.   &      
     &     (jk.lt.kcbot(jl).and.jk.gt.kctop(jl)).and. jk.ge.ikhsmin(jl)                            
          if(llo2(jl))then                                                       
            is=is+1                                                              
          endif                                                                  
        enddo                                                                    
        if(is.eq.0) cycle                                                        
                                                                                 
        ik=jk                                                                    
        icall=2                                                                  
        call cuadjtqn                                           &                  
     &   ( klon, klev, ik, zph, ztenwb, zqenwb, llo2, icall)                        
                                                                                 




        do jl=1,klon                                                        
          if(llo2(jl)) then                                                      
            zttest=0.5*(ptu(jl,jk)+ztenwb(jl,jk))                                
            zqtest=0.5*(pqu(jl,jk)+zqenwb(jl,jk))                                
            zbuo=zttest*(1.+vtmpc1  *zqtest)-                    &                  
     &       ptenh(jl,jk)*(1.+vtmpc1  *pqenh(jl,jk))                               
            zcond(jl)=pqenh(jl,jk)-zqenwb(jl,jk)                                 
            zmftop=-cmfdeps*pmfub(jl)                                            
            if(zbuo.lt.0..and.prfl(jl).gt.10.*zmftop*zcond(jl)) then             
              kdtop(jl)=jk                                                       
              lddraf(jl)=.true.                                                  
              ptd(jl,jk)=zttest                                                  
              pqd(jl,jk)=zqtest                                                  
              pmfd(jl,jk)=zmftop                                                 
              pmfds(jl,jk)=pmfd(jl,jk)*(cpd*ptd(jl,jk)+pgeoh(jl,jk))            
              pmfdq(jl,jk)=pmfd(jl,jk)*pqd(jl,jk)                                
              pdmfdp(jl,jk-1)=-0.5*pmfd(jl,jk)*zcond(jl)                         
              prfl(jl)=prfl(jl)+pdmfdp(jl,jk-1)                                  
            endif                                                                
          endif                                                                  
        enddo                                                                    
                                                                                 
      enddo                                                                      
                                                                                 
      return                                                                     
      end subroutine cudlfsn  





       subroutine cuddrafn                                 &                
     &   ( klon,     klev,    lddraf                      &                                   
     &   , ptenh,    pqenh,    puen,     pven             &                
     &   , pgeo,     pgeoh,    paph,     prfl             &                
     &   , ptd,      pqd,      pud,      pvd,      pmfu   &                
     &   , pmfd,     pmfds,    pmfdq,    pdmfdp,   pmfdde_rate     )                          
                                                                          

                                                                          

                                                                          




                                                                          


                                                                          




                                                                          





                                                                          



                                                                          


                                                                          

                                                                          

                                                                          

                                                                          








                                                                          

                                                                          

                                                                          

                                                                          








                                                                          





      implicit none
      
      integer  klev,klon                                                             
      real     ptenh(klon,klev),       pqenh(klon,klev),   &               
     &         puen(klon,klev),        pven(klon,klev),    &               
     &         pgeoh(klon,klev+1),     paph(klon,klev+1),  &               
     &         pgeo(klon,klev),        pmfu(klon,klev)                    
                                                                          
      real     ptd(klon,klev),         pqd(klon,klev),     &               
     &         pud(klon,klev),         pvd(klon,klev),     &               
     &         pmfd(klon,klev),        pmfds(klon,klev),   &               
     &         pmfdq(klon,klev),       pdmfdp(klon,klev),  &               
     &         prfl(klon)     
      real     pmfdde_rate(klon,klev)                                            
      logical  lddraf(klon)                                               
                                                                          
      real     zdmfen(klon),           zdmfde(klon),       &               
     &         zcond(klon),            zoentr(klon),       &               
     &         zbuoy(klon)                                                
      real     zph(klon)                         
      logical  llo2(klon)                                   
      logical  llo1

      integer  jl,jk
      integer  is,ik,icall,ike, itopde(klon)
      real     zentr,zdz,zzentr,zseen,zqeen,zsdde,zqdde,zdmfdp
      real     zmfdsk,zmfdqk,zbuo,zrain,zbuoyz,zmfduk,zmfdvk
                                                                











      do jl=1,klon                                                   
        zoentr(jl)=0.                                                     
        zbuoy(jl)=0.                                                      
        zdmfen(jl)=0.                                                     
        zdmfde(jl)=0.
      enddo 

      do jk=klev,1,-1
       do jl=1,klon
         pmfdde_rate(jl,jk) = 0.
         if((paph(jl,klev+1)-paph(jl,jk)).lt. 60.e2) itopde(jl) = jk
       end do
      end do                                                              
                                                                 
      do jk=3,klev                                                        
        is=0                                                              
        do jl=1,klon                                                 
          zph(jl)=paph(jl,jk)                                             
          llo2(jl)=lddraf(jl).and.pmfd(jl,jk-1).lt.0.                     
          if(llo2(jl)) then                                               
            is=is+1                                                       
          endif                                                           
        enddo                                                             
        if(is.eq.0) cycle                                                 
                                                                          
        do jl=1,klon                                                 
          if(llo2(jl)) then    
            zentr = entrdd*pmfd(jl,jk-1)*(pgeoh(jl,jk-1)-pgeoh(jl,jk))*zrg
            zdmfen(jl)=zentr                                              
            zdmfde(jl)=zentr                                              
          endif                                                           
        enddo 
                                                            
        do jl=1,klon
          if(llo2(jl)) then
          if(jk.gt.itopde(jl)) then
            zdmfen(jl)=0.
            zdmfde(jl)=pmfd(jl,itopde(jl))*                    &
     &       (paph(jl,jk)-paph(jl,jk-1))/                  &
     &       (paph(jl,klev+1)-paph(jl,itopde(jl)))
          endif
          endif
        enddo

        do jl=1,klon
          if(llo2(jl)) then
          if(jk.le.itopde(jl)) then
            zdz=-(pgeoh(jl,jk-1)-pgeoh(jl,jk))*zrg
            zzentr=zoentr(jl)*zdz*pmfd(jl,jk-1)
            zdmfen(jl)=zdmfen(jl)+zzentr
            zdmfen(jl)=max(zdmfen(jl),0.3*pmfd(jl,jk-1))
            zdmfen(jl)=max(zdmfen(jl),-0.75*pmfu(jl,jk)-   &
   &         (pmfd(jl,jk-1)-zdmfde(jl)))
            zdmfen(jl)=min(zdmfen(jl),0.)
          endif
          endif
        enddo

        do jl=1,klon                                                 
          if(llo2(jl)) then                                               
            pmfd(jl,jk)=pmfd(jl,jk-1)+zdmfen(jl)-zdmfde(jl)
            zseen=(cpd*ptenh(jl,jk-1)+pgeoh(jl,jk-1))*zdmfen(jl)         
            zqeen=pqenh(jl,jk-1)*zdmfen(jl)                               
            zsdde=(cpd*ptd(jl,jk-1)+pgeoh(jl,jk-1))*zdmfde(jl)           
            zqdde=pqd(jl,jk-1)*zdmfde(jl)                                 
            zmfdsk=pmfds(jl,jk-1)+zseen-zsdde                             
            zmfdqk=pmfdq(jl,jk-1)+zqeen-zqdde                             
            pqd(jl,jk)=zmfdqk*(1./min(-cmfcmin,pmfd(jl,jk)))              
            ptd(jl,jk)=(zmfdsk*(1./min(-cmfcmin,pmfd(jl,jk)))-&            
     &                  pgeoh(jl,jk))*rcpd                                
            ptd(jl,jk)=min(400.,ptd(jl,jk))                               
            ptd(jl,jk)=max(100.,ptd(jl,jk))                               
            zcond(jl)=pqd(jl,jk)                                          
          endif                                                           
        enddo                                                             
                                                                          
        ik=jk                                                             
        icall=2                                                           
        call cuadjtqn(klon, klev, ik, zph, ptd, pqd, llo2, icall )                
                                                                          
        do jl=1,klon                                                 
          if(llo2(jl)) then                                               
            zcond(jl)=zcond(jl)-pqd(jl,jk)
            zbuo=ptd(jl,jk)*(1.+vtmpc1  *pqd(jl,jk))-          &             
     &      ptenh(jl,jk)*(1.+vtmpc1  *pqenh(jl,jk))                         
            if(prfl(jl).gt.0..and.pmfu(jl,jk).gt.0.) then                 
              zrain=prfl(jl)/pmfu(jl,jk)                                  
              zbuo=zbuo-ptd(jl,jk)*zrain
            endif                                                         
            if(zbuo.ge.0 .or. prfl(jl).le.(pmfd(jl,jk)*zcond(jl))) then    
              pmfd(jl,jk)=0.
              zbuo=0.                                              
            endif
            pmfds(jl,jk)=(cpd*ptd(jl,jk)+pgeoh(jl,jk))*pmfd(jl,jk)       
            pmfdq(jl,jk)=pqd(jl,jk)*pmfd(jl,jk)                           
            zdmfdp=-pmfd(jl,jk)*zcond(jl)                                 
            pdmfdp(jl,jk-1)=zdmfdp                                        
            prfl(jl)=prfl(jl)+zdmfdp                                      
                                                                          

            zbuoyz=zbuo/ptenh(jl,jk)                                      
            zbuoyz=min(zbuoyz,0.0)                                        
            zdz=-(pgeo(jl,jk-1)-pgeo(jl,jk))
            zbuoy(jl)=zbuoy(jl)+zbuoyz*zdz
            zoentr(jl)=g*zbuoyz*0.5/(1.+zbuoy(jl)) 
            pmfdde_rate(jl,jk) = -zdmfde(jl)
          endif                                                           
        enddo
                                                             
      enddo                                                               
                                                                          
      return                                                              
      end subroutine cuddrafn



       subroutine cuflxn         &                                                 
     &  (  klon,     klev,     ztmst &                                                             
     &  ,  pten,     pqen,     pqsen,    ptenh,    pqenh   &                    
     &  ,  paph,     pap,      pgeoh,    lndj,   ldcum   &                    
     &  ,  kcbot,    kctop,    kdtop,    ktopm2            &                    
     &  ,  ktype,    lddraf                                &                    
     &  ,  pmfu,     pmfd,     pmfus,    pmfds             &                    
     &  ,  pmfuq,    pmfdq,    pmful,    plude             &                    
     &  ,  pdmfup,   pdmfdp,   pdpmel,   plglac            &                    
     &  ,  prain,    pmfdde_rate, pmflxr, pmflxs )                                         
                                                                               

                                                                               


                                                                               


                                                                               



                                                                               
                                                                               



                                                                               





                                                                               

                                                                               


                                                                               

                                                                               









                                                                               

                                                                               

                                                                               

                                                                               

                                                                               

                                                                               










                                                                               

                                                                               






                                                                               




      implicit none                                                     

      integer  klev,klon,ktopm2                                                                         
      real     pten(klon,klev),        ptenh(klon,klev),          &             
     &         pqen(klon,klev),        pqsen(klon,klev),          &             
     &         pqenh(klon,klev),       pap(klon,klev),            &             
     &         paph(klon,klev+1),      pgeoh(klon,klev+1),        &             
     &         plglac(klon,klev)                                               
                                                                               
      real     pmfu(klon,klev),        pmfd(klon,klev),           &             
     &         pmfus(klon,klev),       pmfds(klon,klev),          &             
     &         pmfuq(klon,klev),       pmfdq(klon,klev),          &             
     &         pdmfup(klon,klev),      pdmfdp(klon,klev),         &             
     &         pdpmel(klon,klev),      prain(klon),               &             
     &         pmful(klon,klev),       plude(klon,klev),          &             
     &         pmflxr(klon,klev+1),    pmflxs(klon,klev+1)  
      real     pmfdde_rate(klon,klev)               
      integer  kcbot(klon),            kctop(klon),               &             
     &         kdtop(klon),            ktype(klon)                             
      logical  lddraf(klon),                                &
     &         ldcum(klon)  
      integer  lndj(klon)

      integer  jl,jk
      integer  is,ik,icall,ike,ikb
      real     ztmst,ztaumel,zcons1a,zcons1,zcons2,zcucov,zcpecons
      real     zalfaw,zrfl,zdrfl1,zrnew,zrmin,zrfln,zdrfl,zdenom
      real     zpdr,zpds,zzp,zfac,zsnmlt
      real     rhevap(klon)
      integer  idbas(klon)
      logical  llddraf



      ztaumel=18000.
      zcons1a=cpd/(alf*g*ztaumel)
      zcons2=3./(g*ztmst)
      zcucov=0.05
      zcpecons=5.44e-4/g


      do jl=1,klon
        prain(jl)=0.
        if(.not.ldcum(jl).or.kdtop(jl).lt.kctop(jl)) lddraf(jl)=.false.
        if(.not.ldcum(jl)) ktype(jl)=0
        idbas(jl) = klev
        if(lndj(jl) .eq. 1) then
          rhevap(jl) = 0.7
        else
          rhevap(jl) = 0.9
        end if
      enddo

      ktopm2= 2
      do jk=ktopm2,klev
        ikb = min(jk+1,klev)
        do jl=1,klon
          pmflxr(jl,jk) = 0.
          pmflxs(jl,jk) = 0.
          pdpmel(jl,jk) = 0.
          if(ldcum(jl).and.jk.ge.kctop(jl)) then
            pmfus(jl,jk)=pmfus(jl,jk)-pmfu(jl,jk)*           &
     &       (cpd*ptenh(jl,jk)+pgeoh(jl,jk))
            pmfuq(jl,jk)=pmfuq(jl,jk)-pmfu(jl,jk)*pqenh(jl,jk)
            plglac(jl,jk)=pmfu(jl,jk)*plglac(jl,jk)
            llddraf = lddraf(jl) .and. jk >= kdtop(jl)
            if ( llddraf .and.jk.ge.kdtop(jl)) then
              pmfds(jl,jk) = pmfds(jl,jk)-pmfd(jl,jk) * &
                           (cpd*ptenh(jl,jk)+pgeoh(jl,jk))
              pmfdq(jl,jk) = pmfdq(jl,jk)-pmfd(jl,jk)*pqenh(jl,jk)
            else
              pmfd(jl,jk) = 0.
              pmfds(jl,jk) = 0.
              pmfdq(jl,jk) = 0.
              pdmfdp(jl,jk-1) = 0.
            end if
            if ( llddraf .and. pmfd(jl,jk) < 0. .and. &
               abs(pmfd(jl,ikb)) < 1.e-20 ) then
               idbas(jl) = jk
            end if
          else
            pmfu(jl,jk)=0.
            pmfd(jl,jk)=0.
            pmfus(jl,jk)=0.
            pmfds(jl,jk)=0.
            pmfuq(jl,jk)=0.
            pmfdq(jl,jk)=0.
            pmful(jl,jk)=0.
            plglac(jl,jk)=0.
            pdmfup(jl,jk-1)=0.
            pdmfdp(jl,jk-1)=0.
            plude(jl,jk-1)=0.
          endif
        enddo
      enddo

      do jl=1,klon
        pmflxr(jl,klev+1) = 0.
        pmflxs(jl,klev+1) = 0.
      end do
      do jl=1,klon
        if(ldcum(jl)) then
          ikb=kcbot(jl)
          ik=ikb+1
          zzp=((paph(jl,klev+1)-paph(jl,ik))/                  &
     &     (paph(jl,klev+1)-paph(jl,ikb)))
          if(ktype(jl).eq.3) then
            zzp=zzp**2
          endif
          pmfu(jl,ik)=pmfu(jl,ikb)*zzp
          pmfus(jl,ik)=(pmfus(jl,ikb)-                         &
     &         foelhm(ptenh(jl,ikb))*pmful(jl,ikb))*zzp
          pmfuq(jl,ik)=(pmfuq(jl,ikb)+pmful(jl,ikb))*zzp
          pmful(jl,ik)=0.
        endif
      enddo

      do jk=ktopm2,klev
        do jl=1,klon
          if(ldcum(jl).and.jk.gt.kcbot(jl)+1) then
            ikb=kcbot(jl)+1
            zzp=((paph(jl,klev+1)-paph(jl,jk))/                &
     &       (paph(jl,klev+1)-paph(jl,ikb)))
            if(ktype(jl).eq.3) then
              zzp=zzp**2
            endif
            pmfu(jl,jk)=pmfu(jl,ikb)*zzp
            pmfus(jl,jk)=pmfus(jl,ikb)*zzp
            pmfuq(jl,jk)=pmfuq(jl,ikb)*zzp
            pmful(jl,jk)=0.
          endif
          ik = idbas(jl)
          llddraf = lddraf(jl) .and. jk > ik .and. ik < klev
          if ( llddraf .and. ik == kcbot(jl)+1 ) then
           zzp = ((paph(jl,klev+1)-paph(jl,jk))/(paph(jl,klev+1)-paph(jl,ik)))
           if ( ktype(jl) == 3 ) zzp = zzp*zzp
            pmfd(jl,jk) = pmfd(jl,ik)*zzp
            pmfds(jl,jk) = pmfds(jl,ik)*zzp
            pmfdq(jl,jk) = pmfdq(jl,ik)*zzp
            pmfdde_rate(jl,jk) = -(pmfd(jl,jk-1)-pmfd(jl,jk))
           else if ( llddraf .and. ik /= kcbot(jl)+1 .and. jk == ik+1 ) then
            pmfdde_rate(jl,jk) = -(pmfd(jl,jk-1)-pmfd(jl,jk))
           end if
        enddo
      enddo





        do jk=ktopm2,klev
        do jl=1,klon
          if(ldcum(jl) .and. jk >=kctop(jl)-1 ) then
            prain(jl)=prain(jl)+pdmfup(jl,jk)
            if(pmflxs(jl,jk).gt.0..and.pten(jl,jk).gt.tmelt) then
              zcons1=zcons1a*(1.+0.5*(pten(jl,jk)-tmelt))
              zfac=zcons1*(paph(jl,jk+1)-paph(jl,jk))
              zsnmlt=min(pmflxs(jl,jk),zfac*(pten(jl,jk)-tmelt))
              pdpmel(jl,jk)=zsnmlt
              pqsen(jl,jk)=foeewm(pten(jl,jk)-zsnmlt/zfac)/pap(jl,jk)
            endif
            zalfaw=foealfa(pten(jl,jk))



            if ( pten(jl,jk) < tmelt .and. zalfaw > 0. ) then
              plglac(jl,jk) = plglac(jl,jk)+zalfaw*(pdmfup(jl,jk)+pdmfdp(jl,jk))
              zalfaw = 0.
            end if
            pmflxr(jl,jk+1)=pmflxr(jl,jk)+zalfaw*               &
     &       (pdmfup(jl,jk)+pdmfdp(jl,jk))+pdpmel(jl,jk)
            pmflxs(jl,jk+1)=pmflxs(jl,jk)+(1.-zalfaw)*          &
     &       (pdmfup(jl,jk)+pdmfdp(jl,jk))-pdpmel(jl,jk)
            if(pmflxr(jl,jk+1)+pmflxs(jl,jk+1).lt.0.0) then
              pdmfdp(jl,jk)=-(pmflxr(jl,jk)+pmflxs(jl,jk)+pdmfup(jl,jk))
              pmflxr(jl,jk+1)=0.0
              pmflxs(jl,jk+1)=0.0
              pdpmel(jl,jk)  =0.0
            else if ( pmflxr(jl,jk+1) < 0. ) then
              pmflxs(jl,jk+1) = pmflxs(jl,jk+1)+pmflxr(jl,jk+1)
              pmflxr(jl,jk+1) = 0.
            else if ( pmflxs(jl,jk+1) < 0. ) then
              pmflxr(jl,jk+1) = pmflxr(jl,jk+1)+pmflxs(jl,jk+1)
              pmflxs(jl,jk+1) = 0.
            end if
          endif
        enddo
      enddo
      do jk=ktopm2,klev
        do jl=1,klon
          if(ldcum(jl).and.jk.ge.kcbot(jl)) then
            zrfl=pmflxr(jl,jk)+pmflxs(jl,jk)
            if(zrfl.gt.1.e-20) then
              zdrfl1=zcpecons*max(0.,pqsen(jl,jk)-pqen(jl,jk))*zcucov* &
     &         (sqrt(paph(jl,jk)/paph(jl,klev+1))/5.09e-3*             &
     &         zrfl/zcucov)**0.5777*                                   &
     &         (paph(jl,jk+1)-paph(jl,jk))
              zrnew=zrfl-zdrfl1
              zrmin=zrfl-zcucov*max(0.,rhevap(jl)*pqsen(jl,jk) &
     &              -pqen(jl,jk)) *zcons2*(paph(jl,jk+1)-paph(jl,jk))
              zrnew=max(zrnew,zrmin)
              zrfln=max(zrnew,0.)
              zdrfl=min(0.,zrfln-zrfl)
              zdenom=1./max(1.e-20,pmflxr(jl,jk)+pmflxs(jl,jk))
              zalfaw=foealfa(pten(jl,jk))
              if ( pten(jl,jk) < tmelt ) zalfaw = 0.
              zpdr=zalfaw*pdmfdp(jl,jk)
              zpds=(1.-zalfaw)*pdmfdp(jl,jk)
              pmflxr(jl,jk+1)=pmflxr(jl,jk)+zpdr         &
     &         +pdpmel(jl,jk)+zdrfl*pmflxr(jl,jk)*zdenom
              pmflxs(jl,jk+1)=pmflxs(jl,jk)+zpds         &
     &         -pdpmel(jl,jk)+zdrfl*pmflxs(jl,jk)*zdenom
              pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
              if ( pmflxr(jl,jk+1)+pmflxs(jl,jk+1) < 0. ) then
                pdmfup(jl,jk) = pdmfup(jl,jk)-(pmflxr(jl,jk+1)+pmflxs(jl,jk+1))
                pmflxr(jl,jk+1) = 0.
                pmflxs(jl,jk+1) = 0.
                pdpmel(jl,jk)   = 0.
              else if ( pmflxr(jl,jk+1) < 0. ) then
                pmflxs(jl,jk+1) = pmflxs(jl,jk+1)+pmflxr(jl,jk+1)
                pmflxr(jl,jk+1) = 0.
              else if ( pmflxs(jl,jk+1) < 0. ) then
                pmflxr(jl,jk+1) = pmflxr(jl,jk+1)+pmflxs(jl,jk+1)
                pmflxs(jl,jk+1) = 0.
              end if
            else
              pmflxr(jl,jk+1)=0.0
              pmflxs(jl,jk+1)=0.0
              pdmfdp(jl,jk)=0.0
              pdpmel(jl,jk)=0.0
            endif
          endif
        enddo
      enddo
                                      
      return
      end subroutine cuflxn 



     subroutine cudtdqn(klon,klev,ktopm2,kctop,kdtop,ldcum, &
                     lddraf,ztmst,paph,pgeoh,pgeo,pten,ptenh,pqen,  &
                     pqenh,pqsen,plglac,plude,pmfu,pmfd,pmfus,pmfds, &
                     pmfuq,pmfdq,pmful,pdmfup,pdmfdp,pdpmel,ptent,ptenq,pcte)
    implicit none
    integer  klon,klev,ktopm2
    integer  kctop(klon),  kdtop(klon)
    logical  ldcum(klon),  lddraf(klon)
    real     ztmst
    real     paph(klon,klev+1), pgeoh(klon,klev+1)
    real     pgeo(klon,klev),   pten(klon,klev), &
             pqen(klon,klev),   ptenh(klon,klev),&
             pqenh(klon,klev),  pqsen(klon,klev),&
             plglac(klon,klev), plude(klon,klev)
    real     pmfu(klon,klev),   pmfd(klon,klev),&
             pmfus(klon,klev),  pmfds(klon,klev),&
             pmfuq(klon,klev),  pmfdq(klon,klev),&
             pmful(klon,klev),  pdmfup(klon,klev),&
             pdpmel(klon,klev), pdmfdp(klon,klev)
    real     ptent(klon,klev),  ptenq(klon,klev)
    real     pcte(klon,klev)

    integer  jk , ik , jl
    real     zalv , zzp
    real     zdtdt(klon,klev) , zdqdt(klon,klev) , zdp(klon,klev)


    do jk = 1 , klev
      do jl = 1, klon
        if ( ldcum(jl) ) then
          zdp(jl,jk) = g/(paph(jl,jk+1)-paph(jl,jk))
        end if
      end do
    end do



    do jk = ktopm2 , klev
      if ( jk < klev ) then
        do jl = 1,klon
          if ( ldcum(jl) ) then
            zalv = foelhm(pten(jl,jk))
            zdtdt(jl,jk) = zdp(jl,jk)*rcpd * &
              (pmfus(jl,jk+1)-pmfus(jl,jk)+pmfds(jl,jk+1) - &
               pmfds(jl,jk)+alf*plglac(jl,jk)-alf*pdpmel(jl,jk) - &
               zalv*(pmful(jl,jk+1)-pmful(jl,jk)-plude(jl,jk)-pdmfup(jl,jk)-pdmfdp(jl,jk)))
            zdqdt(jl,jk) = zdp(jl,jk)*(pmfuq(jl,jk+1) - &
              pmfuq(jl,jk)+pmfdq(jl,jk+1)-pmfdq(jl,jk)+pmful(jl,jk+1) - &
              pmful(jl,jk)-plude(jl,jk)-pdmfup(jl,jk)-pdmfdp(jl,jk))
          end if
        end do
      else
        do jl = 1,klon
          if ( ldcum(jl) ) then
            zalv = foelhm(pten(jl,jk))
            zdtdt(jl,jk) = -zdp(jl,jk)*rcpd * &
              (pmfus(jl,jk)+pmfds(jl,jk)+alf*pdpmel(jl,jk) - &
               zalv*(pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)+plude(jl,jk)))
            zdqdt(jl,jk) = -zdp(jl,jk)*(pmfuq(jl,jk) + plude(jl,jk) + &
              pmfdq(jl,jk)+(pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
          end if
        end do
      end if
    end do



    do jk = ktopm2 , klev
     do jl = 1, klon
       if ( ldcum(jl) ) then
         ptent(jl,jk) = ptent(jl,jk) + zdtdt(jl,jk)
         ptenq(jl,jk) = ptenq(jl,jk) + zdqdt(jl,jk)
         pcte(jl,jk)  = zdp(jl,jk)*plude(jl,jk)
       end if
     end do
    end do

    return
  end subroutine cudtdqn



    subroutine cududvn(klon,klev,ktopm2,ktype,kcbot,kctop,ldcum,  &
                    ztmst,paph,puen,pven,pmfu,pmfd,puu,pud,pvu,pvd,ptenu, &
                    ptenv)
    implicit none
    integer  klon,klev,ktopm2
    integer  ktype(klon), kcbot(klon), kctop(klon)
    logical  ldcum(klon)
    real     ztmst
    real     paph(klon,klev+1)
    real     puen(klon,klev),  pven(klon,klev),&
             pmfu(klon,klev),  pmfd(klon,klev),&
             puu(klon,klev),   pud(klon,klev),&
             pvu(klon,klev),   pvd(klon,klev)
    real     ptenu(klon,klev), ptenv(klon,klev)

    real     zuen(klon,klev) , zven(klon,klev) , zmfuu(klon,klev), &
             zmfdu(klon,klev), zmfuv(klon,klev), zmfdv(klon,klev)

    integer  ik , ikb , jk , jl
    real     zzp, zdtdt
   
    real     zdudt(klon,klev), zdvdt(klon,klev), zdp(klon,klev)

    do jk = 1 , klev
      do jl = 1, klon
        if ( ldcum(jl) ) then
          zuen(jl,jk) = puen(jl,jk)
          zven(jl,jk) = pven(jl,jk)
          zdp(jl,jk) = g/(paph(jl,jk+1)-paph(jl,jk))
        end if
      end do
    end do



    do jk = ktopm2 , klev
      ik = jk - 1
      do jl = 1,klon
        if ( ldcum(jl) ) then
          zmfuu(jl,jk) = pmfu(jl,jk)*(puu(jl,jk)-zuen(jl,ik))
          zmfuv(jl,jk) = pmfu(jl,jk)*(pvu(jl,jk)-zven(jl,ik))
          zmfdu(jl,jk) = pmfd(jl,jk)*(pud(jl,jk)-zuen(jl,ik))
          zmfdv(jl,jk) = pmfd(jl,jk)*(pvd(jl,jk)-zven(jl,ik))
        end if
      end do
    end do

      do jk = ktopm2 , klev
        do jl = 1, klon
          if ( ldcum(jl) .and. jk > kcbot(jl) ) then
            ikb = kcbot(jl)
            zzp = ((paph(jl,klev+1)-paph(jl,jk))/(paph(jl,klev+1)-paph(jl,ikb)))
            if ( ktype(jl) == 3 ) zzp = zzp*zzp
            zmfuu(jl,jk) = zmfuu(jl,ikb)*zzp
            zmfuv(jl,jk) = zmfuv(jl,ikb)*zzp
            zmfdu(jl,jk) = zmfdu(jl,ikb)*zzp
            zmfdv(jl,jk) = zmfdv(jl,ikb)*zzp
          end if
        end do
      end do



    do jk = ktopm2 , klev
      if ( jk < klev ) then
        ik = jk + 1
        do jl = 1,klon
          if ( ldcum(jl) ) then
            zdudt(jl,jk) = zdp(jl,jk) * &
                          (zmfuu(jl,ik)-zmfuu(jl,jk)+zmfdu(jl,ik)-zmfdu(jl,jk))
            zdvdt(jl,jk) = zdp(jl,jk) * &
                          (zmfuv(jl,ik)-zmfuv(jl,jk)+zmfdv(jl,ik)-zmfdv(jl,jk))
          end if
        end do
      else
        do jl = 1,klon
          if ( ldcum(jl) ) then
            zdudt(jl,jk) = -zdp(jl,jk)*(zmfuu(jl,jk)+zmfdu(jl,jk))
            zdvdt(jl,jk) = -zdp(jl,jk)*(zmfuv(jl,jk)+zmfdv(jl,jk))
          end if
        end do
      end if
    end do



    do jk = ktopm2 , klev
      do jl = 1, klon
        if ( ldcum(jl) ) then
          ptenu(jl,jk) = ptenu(jl,jk) + zdudt(jl,jk)
          ptenv(jl,jk) = ptenv(jl,jk) + zdvdt(jl,jk)
        end if
      end do
    end do

    return
    end subroutine cududvn



      subroutine cuadjtqn                 &                                           
     &    (klon, klev, kk, psp, pt, pq, ldflag,  kcall)                           




                                                                                    









                                                                                    



                                                                                    








                                                                                    

                                                                                    

                                                                                    






                                                                                    

                                                                                    
      implicit none   
                 
      integer  klev,klon                                                                   
      real     pt(klon,klev),          pq(klon,klev),  &                             
     &         psp(klon)                                                            
      logical  ldflag(klon)                                                         

      integer  jl,jk
      integer  isum,kcall,kk
      real     zqmax,zqsat,zcor,zqp,zcond,zcond1,zl,zi,zf



      zqmax=0.5  
                                                                                    


                                                                                    
      if ( kcall == 1 ) then
      do jl = 1,klon
        if ( ldflag(jl) ) then
          zqp = 1./psp(jl)
          zl = 1./(pt(jl,kk)-c4les)
          zi = 1./(pt(jl,kk)-c4ies)
          zqsat = c2es*(foealfa(pt(jl,kk))*exp(c3les*(pt(jl,kk)-tmelt)*zl) + &
                (1.-foealfa(pt(jl,kk)))*exp(c3ies*(pt(jl,kk)-tmelt)*zi))
          zqsat = zqsat*zqp
          zqsat = min(0.5,zqsat)
          zcor = 1. - vtmpc1*zqsat
          zf = foealfa(pt(jl,kk))*r5alvcp*zl**2 + &
               (1.-foealfa(pt(jl,kk)))*r5alscp*zi**2
          zcond = (pq(jl,kk)*zcor**2-zqsat*zcor)/(zcor**2+zqsat*zf)
          if ( zcond > 0. ) then
            pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond
            pq(jl,kk) = pq(jl,kk) - zcond
            zl = 1./(pt(jl,kk)-c4les)
            zi = 1./(pt(jl,kk)-c4ies)
            zqsat = c2es*(foealfa(pt(jl,kk)) * &
              exp(c3les*(pt(jl,kk)-tmelt)*zl)+(1.-foealfa(pt(jl,kk))) * &
              exp(c3ies*(pt(jl,kk)-tmelt)*zi))
            zqsat = zqsat*zqp
            zqsat = min(0.5,zqsat)
            zcor = 1. - vtmpc1*zqsat
            zf = foealfa(pt(jl,kk))*r5alvcp*zl**2 + &
                 (1.-foealfa(pt(jl,kk)))*r5alscp*zi**2
            zcond1 = (pq(jl,kk)*zcor**2-zqsat*zcor)/(zcor**2+zqsat*zf)
            if ( abs(zcond) < 1.e-20 ) zcond1 = 0.
            pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk) - zcond1
          end if
        end if
      end do
      elseif ( kcall == 2 ) then
      do jl = 1,klon
        if ( ldflag(jl) ) then
          zqp = 1./psp(jl)
          zqsat = foeewm(pt(jl,kk))*zqp
          zqsat = min(0.5,zqsat)
          zcor = 1./(1.-vtmpc1*zqsat)
          zqsat = zqsat*zcor
          zcond = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*foedem(pt(jl,kk)))
          zcond = min(zcond,0.)
          pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond
          pq(jl,kk) = pq(jl,kk) - zcond
          zqsat = foeewm(pt(jl,kk))*zqp
          zqsat = min(0.5,zqsat)
          zcor = 1./(1.-vtmpc1*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*foedem(pt(jl,kk)))
          if ( abs(zcond) < 1.e-20 ) zcond1 = min(zcond1,0.)
          pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond1
          pq(jl,kk) = pq(jl,kk) - zcond1
        end if
      end do
      else if ( kcall == 0 ) then
      do jl = 1,klon
        zqp = 1./psp(jl)
        zqsat = foeewm(pt(jl,kk))*zqp
        zqsat = min(0.5,zqsat)
        zcor = 1./(1.-vtmpc1*zqsat)
        zqsat = zqsat*zcor
        zcond1 = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*foedem(pt(jl,kk)))
        pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond1
        pq(jl,kk) = pq(jl,kk) - zcond1
        zqsat = foeewm(pt(jl,kk))*zqp
        zqsat = min(0.5,zqsat)
        zcor = 1./(1.-vtmpc1*zqsat)
        zqsat = zqsat*zcor
        zcond1 = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*foedem(pt(jl,kk)))
        pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond1
        pq(jl,kk) = pq(jl,kk) - zcond1
      end do
      end if

      return
      end subroutine cuadjtqn



      subroutine cubasmcn &
     &    (klon,     klev,     klevm1,  kk,     pten,&
     &     pqen,     pqsen,    puen,    pven,   pverv,&
     &     pgeo,     pgeoh,    ldcum,   ktype,  klab,  plrain,&
     &     pmfu,     pmfub,    kcbot,   ptu,&
     &     pqu,      plu,      puu,     pvu,    pmfus,&
     &     pmfuq,    pmful,    pdmfup)
      implicit none


















      real     pten(klon,klev),        pqen(klon,klev),&
     &         puen(klon,klev),        pven(klon,klev),&
     &         pqsen(klon,klev),       pverv(klon,klev),&
     &         pgeo(klon,klev),        pgeoh(klon,klev+1)
      real     ptu(klon,klev),         pqu(klon,klev),&
     &         puu(klon,klev),         pvu(klon,klev),&
     &         plu(klon,klev),         pmfu(klon,klev),&
     &         pmfub(klon),   &         
     &         pmfus(klon,klev),       pmfuq(klon,klev),&
     &         pmful(klon,klev),       pdmfup(klon,klev),&
     &         plrain(klon,klev)
      integer  ktype(klon),            kcbot(klon),&
     &         klab(klon,klev)
      logical  ldcum(klon)

      integer  jl,kk,klev,klon,klevp1,klevm1
      real     zzzmb



         do jl=1,klon
          if(.not.ldcum(jl) .and. klab(jl,kk+1).eq.0) then
            if(lmfmid .and. pqen(jl,kk) .gt. 0.80*pqsen(jl,kk).and. &
              pgeo(jl,kk)*zrg .gt. 5.0e2  .and.  &
     &        pgeo(jl,kk)*zrg .lt. 1.0e4 )  then
            ptu(jl,kk+1)=(cpd*pten(jl,kk)+pgeo(jl,kk)-pgeoh(jl,kk+1))&
     &                          *rcpd
            pqu(jl,kk+1)=pqen(jl,kk)
            plu(jl,kk+1)=0.
            zzzmb=max(cmfcmin,-pverv(jl,kk)*zrg)
            zzzmb=min(zzzmb,cmfcmax)
            pmfub(jl)=zzzmb
            pmfu(jl,kk+1)=pmfub(jl)
            pmfus(jl,kk+1)=pmfub(jl)*(cpd*ptu(jl,kk+1)+pgeoh(jl,kk+1))
            pmfuq(jl,kk+1)=pmfub(jl)*pqu(jl,kk+1)
            pmful(jl,kk+1)=0.
            pdmfup(jl,kk+1)=0.
            kcbot(jl)=kk
            klab(jl,kk+1)=1
            plrain(jl,kk+1)=0.0
            ktype(jl)=3
          end if
         end if
        end do
      return
      end subroutine cubasmcn



      subroutine cuentrn(klon,klev,kk,kcbot,ldcum,ldwork, &
                    pgeoh,pmfu,pdmfen,pdmfde)
       implicit none
       integer  klon,klev,kk
       integer  kcbot(klon)
       logical  ldcum(klon)
       logical  ldwork
       real  pgeoh(klon,klev+1)
       real  pmfu(klon,klev) 
       real  pdmfen(klon) 
       real  pdmfde(klon) 
       logical  llo1
       integer  jl
       real  zdz , zmf
       real  zentr(klon)



    if ( ldwork ) then
      do jl = 1,klon
        pdmfen(jl) = 0.
        pdmfde(jl) = 0.
        zentr(jl) = 0.
      end do



      do jl = 1, klon
        if ( ldcum(jl) ) then
          zdz = (pgeoh(jl,kk)-pgeoh(jl,kk+1))*zrg
          zmf = pmfu(jl,kk+1)*zdz
          llo1 = kk < kcbot(jl)
          if ( llo1 ) then
            pdmfen(jl) = zentr(jl)*zmf
            pdmfde(jl) = 0.75e-4*zmf
          end if
        end if
      end do
    end if
    end subroutine cuentrn



      real function foealfa(tt)








        implicit none
        real tt
         foealfa = min(1.,((max(rtice,min(rtwat,tt))-rtice) &
     &  /(rtwat-rtice))**2)

        return
      end function foealfa

      real function foelhm(tt)
        implicit none
        real tt
        foelhm = foealfa(tt)*alv + (1.-foealfa(tt))*als
      return
      end function foelhm

      real function foeewm(tt)
        implicit none
        real tt
        foeewm  = c2es * &
     &     (foealfa(tt)*exp(c3les*(tt-tmelt)/(tt-c4les))+ &
     &     (1.-foealfa(tt))*exp(c3ies*(tt-tmelt)/(tt-c4ies)))
      return
      end function foeewm

      real function foedem(tt)
        implicit none
        real tt
        foedem  = foealfa(tt)*r5alvcp*(1./(tt-c4les)**2)+ &
     &              (1.-foealfa(tt))*r5alscp*(1./(tt-c4ies)**2)
      return
      end function foedem

      real function foeldcpm(tt)
        implicit none
        real tt
        foeldcpm = foealfa(tt)*ralvdcp+ &
     &        (1.-foealfa(tt))*ralsdcp
      return
      end function  foeldcpm


end module module_convection_prep
