       subroutine em (Nm11t11, Nm11t12, Nm11t22, Nm12t11, Nm12t12,
     & Nm12t22,Nm22t11, Nm22t12, Nm22t22,del,
     &f11,f12,f21,f22)
       double precision f11,f12,f21,f22,fo,rnind,emdel,del,r2a3
       double precision Nm11t11,Nm11t12,Nm11t22,Nm12t11,Nm12t12,Nm12t22
     &,Nm22t11, Nm22t12, Nm22t22
       double precision ft1,ft2,fm1,fm2,ft1E,fm1E

         rnind=Nm11t11+ Nm11t12+ Nm11t22+Nm12t11+Nm12t12+Nm12t22
     &+Nm22t11+ Nm22t12+ Nm22t22
c            print*,rnind
        ft1=2.*Nm11t11+2.*Nm11t12+2.*Nm11t22+Nm12t11
     &+Nm12t12+ Nm12t22
        ft1=ft1/(2.*rnind)
        ft2=1.-ft1

         fm1=2.*Nm11t11+ 2.*Nm12t11+2.*Nm22t11+Nm11t12+Nm12t12
     &+Nm22t12
        fm1=fm1/(2.*rnind)
        fm2=1.-fm1
c        print*, 'ft1,fm1 EM',ft1,fm1                

       f11=ft1*fm1
       f12=ft1*fm2
       f21=ft2*fm1
       f22=ft2*fm2

        nj=0

        rnind=rnind*2.

 121   continue

        fo=(f11*f22)/((f11*f22)+(f12*f21))
        f11=((2.*Nm11t11)+Nm11t12+Nm12t11+fo*Nm12t12)/rnind
        f12=((2.*Nm11t22)+Nm11t12+Nm12t22+(1.-fo)*Nm12t12)/rnind
        f21=((2.*Nm22t11)+Nm12t11+Nm22t12+(1.-fo)*Nm12t12)/rnind
        f22=((2.*Nm22t22)+Nm12t22+Nm22t12+fo*Nm12t12)/rnind
c         print*, 'fo, f11,f12,f21,f22 ', fo, f11,f12,f21,f22 
         nj=nj+1

         emdel=f11*f22-f12*f21

          del=emdel
         if (nj.eq.1) d1=del
         if (nj.eq.1) go to 121

          if (abs(del-d1).gt.0.0000001.and.nj.le.10000) then
c          write(56,*) d1                                                                                                                 
          d1=del
          go to 121
          else
          ft1E=f11+f12
          fm1E=f11+f21

          return
          endif
         stop
         end
