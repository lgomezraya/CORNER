       program corner
       integer ni,sizew,pr
       character*3 encas,encasa
c     pr option to print individual marker pairs pr=1 means print
       parameter (ni=435)
c   ni=total number of individuals
       integer nmar ,nsnp(18),non1(1000),non2(1000)
       integer nac1(100),nac2(100),mma
       double precision rz,rb,rb1(18),rb2(18),rema
 
c     nmar=maximum total number of markers
c     marmax= maximum number markers within any given chromosome       
       parameter (nmar= 20339,maxmar= 5000)
       double precision sumi,sumi1,nsumi,dsumi,sumi2,fi(nmar),uuA,uvA
c     real, allocatable :: uu(:,:)
       double precision Nm11t11x, Nm11t12x, Nm11t22x,
     &  Nm12t11x, Nm12t12x, Nm12t22x
     &,Nm22t11x, Nm22t12x, Nm22t22x,rtot
       double precision uu(nmar),r1,r11,suma(ni,nmar),tsuma(nmar,ni)
       double precision ma(ni,ni),imahw(ni,ni),tsu(nmar,ni),ttt
       double precision lan(ni),land(ni),dr1,su(ni,nmar),gs(ni,nmar)
       double precision mx(ni,ni),mx1(ni,ni),mx2(ni,ni),tgs(nmar,ni)
       double precision ff,fff,angle, GHW(ni,ni),GGS(ni,ni),ort
       double precision mx3(ni,ni),r(ni),r2(ni),abserr, maxr2(3,3)
       double precision mahw(ni,ni), zhw(ni,nmar), tzhw(nmar,ni),traceD
c     zhw1 auxiliar  matrix to help removing high LD markers
       double precision f11x,f22x,f12x,f21x,a1,a2,b1,b2,ri
       double precision f11,f12,f21,f22,del
       double precision cor1,cor2,cor3,cor4,cor5,cor7
       integer po1(18,nmar)
       double precision zhw1(ni,nmar), tzhw1(nmar,ni)
       double precision manohw(ni,ni),f,varhw,varhwNO,traceG
       double precision orto(ni,nmar) ,torto(nmar,ni)
       integer pox(nmar),cro(nmar),croG(nmar)
       double precision f11a(nmar),f12a(nmar),f22a(nmar)
       character*1 alle1, alle2,al1(18,ni,maxmar),al2(18,ni,maxmar)
       character*14 mar
       character*3 prin
       character*2 x1,x2,x3,ival,iva(18,ni,maxmar)
       integer k,posa,in(ni),nni,positi(nmar),ac(nmar),cr(nmar)
       integer  ind,pos,nchro,gr,po(18,maxmar),lpc,lp,nnnG,ppa(nmar)
c       corne is a matrix with corner situation for any pair of markers: 1 is corner 0 is not corner     
       double precision fr(3,18,4000),orA,orD,rx
       integer  corne(maxmar,maxmar),cori(maxmar,maxmar)
       double precision  rem(maxmar,maxmar),rcor(maxmar,maxmar),ndif,dif       
c        character*1  a1,a2,b1,b2,c1,c2
       character*2 aa,ab,bb, ge(3,18,4000),sex,aax,bbx,ge1(3,18,4000)
c

C      unit 18 must be SORTED BY SNP and then INDIVIDUAL
       open (unit=88,file='corner.output')
       open (unit=25,file='corner.detail')
       open (unit=81,file='corner.all')       
       open (unit=89,file='consecutive')       
       open (unit=18,file='snpv11.all')
       open (unit=10,file='freq.all')
       open (unit=77,file='nonrec')
 222   print*, "ENTER CHROMOSOME TO ANALIZE"
       print*, ""
       print*, "Chromosome 1 : enter 1"
       print*, "Chromosome 2 : enter 2"       
       print*, "......."
       print*, "......."        
       print*, "Chromosome 18 : enter 18"
       print*, "Chromosomes 1-18 : enter 99"                    
       read*,nc
        nfd=0
       if(nc.ge.1.and.nc.le.18) then
          nc1=nc
          nc2=nc
          nd=1 

      print*,"Do you want partial results in output file:corner.detail?
     %"
       print*, " "       
       print*, " Enter yes or not"
       print*, ""        
       print*, ""        
       print*, ""        
       print*, ""        
       read*,prin       
       if(prin.eq."yes") pr=1
       if(prin.eq."not") pr=0
       if(prin.eq."Yes") pr=1
       if(prin.eq."Not") pr=0       
       if(prin.eq."y") pr=1
       if(prin.eq."n") pr=0                
       endif

       if(nc.eq.99) then
          nc1=1
          nc2=18
          nd=1
          pr=0
        endif
        if(nd.eq.0) then
           print*, " "
           print*, " "
           print*, " "
           print*, " "                      
           print*, " ENTER A VALID CHROMOSOME "
           print*, " "
           print*, " "
           print*, " "
           print*, " "           
           go to 222
           
         endif 
       
        ncount=0
        cor1=0
        cor2=0
        cor3=0
        cor4=0
        cor5=0
        cor7=0               
        corne=0
        cori=0
       po=0   
       non1=0
       non2=0
       nac1=0
       nac2=0       
        rb1=0
        rb2=0       
 80    nsnp=0
    
        kk=0
        mx=0
        mx1=0
        mx2=0
        mx3=0
        ma=0        
        uu=0 
        rdif=0
        dif=0
c     READ GENOTYPE FREQUENCIES

 1      continue        
       
         read (10,111,end=112) nchro,pos,k,aa,ab,bb,f11,f12,f22
 111     format (i3,2x,i10, 2x,i1,x, a2,x,a2,x, a2,x,3(f6.4,2x))

         nsnp(nchro)=nsnp(nchro)+1
         

          f=f11+0.5*f12
          mahw=0
          manohw=0
          if (nchroa.ne.nchro.and.nchro.ne.1) rb2(nchroa)=po(nchroa,kk)         
          if (nchroa.ne.nchro) kk=0
          
           kk=kk+1
           nchroa=nchro
           if(kk.eq.1) rb1(nchro)=pos   
           po(nchro,kk)=pos
           fr(1,nchro,kk)=f11
           fr(2,nchro,kk)=f12
           fr(3,nchro,kk)=f22

           ge(1,nchro,kk)=aa
           ge(2,nchro,kk)=ab
           ge(3,nchro,kk)=bb


           go to 1

 112       continue
       rb2(nchro)=po(nchro,kk)                    
       f=0  
       j=0
 12    j=j+1
       read (18,*,end=11) 

       go to 12

 11    continue
       nrec=j-1
       print*, " SNPs # PER CHROMOSOMES "
      print*, "  "       
       print*,' number of records in file ', nrec
       ntotal=0
       print*,' number of individuals ', ni
       do i=1,18       
       print*,' number of SNPs per chromosome ', i,': ', nsnp(i)
       ntotal=ntotal+nsnp(i)
      enddo

      print*, " "
      print*, " Total number of SNPs",ntotal
      print*, " "
      print*, " "
      print*, " LENGTH OF CHROMOSOMES "
      print*, " "            
c      print*,rb1
c      print*,rb2      
c      rz=0      
      do i=1,18       
         rb1(i)=rb1(i)/1000000.
         rb2(i)=rb2(i)/1000000.                     
         rb=rb2(i)-rb1(i)
            rz=rz+rb

      PRINT'(a21,I3,a3,f8.3,a3)','length for chromosome',i,': ',rb,' Mb'
      enddo      
       print*, " "
       rx=rz
       print*,'length of genome ',rz, " Mb"
       rewind 18

       do k=1,18 
          njc=0
          kiu=0
          do i=1,nsnp(k)      


         do 2222 j=1,ni
       
           read(18,1289) ind,pos,gr,mar,nchro,alle1,alle2
 1289      format(i5,i11,i2,a12,i5,2x,a1,a1) 

         al1(nchro,j,i)=alle1
         al2(nchro,j,i)=alle2
         iva(nchro,j,i)=trim(alle1//alle2)
         po1(nchro,i)=pos
       
 2222  continue
          enddo
          njc=0
       enddo

       ri=ni
c     k is chromosome number
       nd=0

       
          do k=nc1,nc2
            sto=0
            corne=0
            cori=0
            rem=0
            rcor=0
            do 27 i1=1,nsnp(k)-1    
      do 87 j1=i1+1,nsnp(k)
            Nm11t11=0
            Nm11t12=0
            Nm11t22=0
            Nm12t11=0
            Nm12t12=0
            Nm12t22=0
            Nm22t11=0
            Nm22t12=0
            Nm22t22=0            
         do 130 jx=1,ni

           if(iva(k,jx,i1).eq.ge(1,k,i1)) then

        if(iva(k,jx,j1).eq.ge(1,k,j1)) then
           Nm11t11=Nm11t11+1
          endif

        if(iva(k,jx,j1).eq.ge(2,k,j1)) then
            Nm11t12=Nm11t12+1

         endif

        if(iva(k,jx,j1).eq.ge(3,k,j1)) then
            Nm11t22=Nm11t22+1
         endif         
         
           endif

c****************
        if(iva(k,jx,i1).eq.ge(2,k,i1)) then

        if(iva(k,jx,j1).eq.ge(1,k,j1)) then
            Nm12t11=Nm12t11+1
         endif

        if(iva(k,jx,j1).eq.ge(2,k,j1)) then
            Nm12t12=Nm12t12+1
         endif

        if(iva(k,jx,j1).eq.ge(3,k,j1)) then
            Nm12t22=Nm12t22+1
         endif         
       
           endif

c****************
        if(iva(k,jx,i1).eq.ge(3,k,i1)) then

        if(iva(k,jx,j1).eq.ge(1,k,j1)) then
            Nm22t11=Nm22t11+1
         endif

        if(iva(k,jx,j1).eq.ge(2,k,j1)) then
            Nm22t12=Nm22t12+1
         endif

        if(iva(k,jx,j1).eq.ge(3,k,j1)) then
            Nm22t22=Nm22t22+1
         endif         
         
           endif
 130    continue
       if(Nm11t11+Nm11t12+Nm12t11.eq.0) then        
          if(Nm12t22+Nm22t12+Nm22t22.eq.0) then
          go to 27           
          endif
          endif

      if(Nm11t12+Nm11t22+Nm12t22.eq.0) then        
         if(Nm12t11+Nm22t11+Nm22t12.eq.0) then
           go to 27

         endif
      endif

c      print*, " "
      if (pr.eq.1) then
      write(25,*) " Chromosome", k
     %, "Counts markers ", i1, j1, "  in positions "
     @ ,po1(k,i1), "and ",po1(k,j1)

      write(25,*) " "
      write(25,*) " "
      write(25,*) "             ", ge(1,k,j1),"          ",ge(2,k,j1) ,
     $"            ",ge(3,k,j1) 
      write(25,*) "_________________________________________"  
      write(25,*) ge(1,k,i1),"|", Nm11t11, Nm11t12, Nm11t22
      write(25,*) ge(2,k,i1), "|", Nm12t11 ,Nm12t12, Nm12t22
      write(25,*) ge(3,k,i1),"|", Nm22t11, Nm22t12,  Nm22t22
      write(25,*)"_________________________________________"  
      write(25,*) " "
      write(25,*) " total counts"
      write(25,*) " "
      
      
      write(25,*) Nm11t11+ Nm11t12+  Nm11t22
     % +     Nm12t11 +Nm12t12+  Nm12t22
     % +     Nm22t11+ Nm22t12+  Nm22t22,ni
      write(25,*) " "
      endif
      
      if (Nm12t11+Nm11t11+Nm11t12.eq.0) then
         cori(i1,j1)=9                  
         corne(i1,j1)=1
      cor1=cor1+1
        go to 6566
      endif
      

      if (Nm11t12+Nm11t22+Nm12t22.eq.0) then 
      if (Nm12t11+Nm22t11+Nm22t12.eq.0) then
         corne(i1,j1)=1                  
       if(j1.eq.i1+1)  cor5=cor5+1
      else
         cor2=cor2+1
         corne(i1,j1)=1         
         
         cori(i1,j1)=9                           
       endif
        go to 6566
      endif
      
      if (Nm12t11+Nm22t11+Nm22t12.eq.0) then
         cor3=cor3+1
         corne(i1,j1)=1
         cori(i1,j1)=9                           
        
        go to 6566
      endif         
      if (Nm22t12+Nm22t22+Nm12t22.eq.0) then
         cor4=cor4+1
         corne(i1,j1)=1
         cori(i1,j1)=9                  
      if (Nm12t11+Nm11t11+Nm11t12.eq.0) then
         corne(i1,j1)=1
         cor7=cor7+1
      endif
        endif

 6566   continue
      ncount=ncount+1
       a11=0
       a12=0
       a22=0
       b11=0
       b12=0
       b22=0
        ort=0
        Nm11t11x=Nm11t11
        Nm11t12x=Nm11t12
        Nm11t22x=Nm11t22
        Nm12t11x=Nm12t11
        Nm12t12x=Nm12t12
        Nm12t22x=Nm12t22
        Nm22t11x=Nm22t11
        Nm22t12x=Nm22t12
        Nm22t22x=Nm22t22

        r1=0
         
        a1=0
        a2=0
        b1=0
        b2=0
        Nm11t11x=Nm11t11
        Nm11t12x=Nm11t12
        Nm11t22x=Nm11t22
        Nm12t11x=Nm12t11
        Nm12t12x=Nm12t12
        Nm12t22x=Nm12t22
        Nm22t11x=Nm22t11
        Nm22t12x=Nm22t12
        Nm22t22x=Nm22t22        
      a1=2.*(Nm11t11x+Nm11t12x+Nm11t22x)+(Nm12t11x+Nm12t12x+Nm12t22x)
      a2=2.*(Nm22t11x+Nm22t12x+Nm22t22x)+(Nm12t11x+Nm12t12x+Nm12t22x)
      b1=2.*(Nm11t11x+Nm12t11x+Nm22t11x)+(Nm11t12x+Nm12t12x+Nm22t12x)
      b2=2.*(Nm11t22x+Nm12t22x+Nm22t22x)+(Nm11t12x+Nm12t12x+Nm22t12x)
      a1=a1/(2.*ri)
      a2=a2/(2.*ri)
      b1=b1/(2.*ri)
      b2=b2/(2.*ri)
        r11=0   
      
      
       if(Nm11t11+Nm11t12+Nm12t11.eq.0) then        

      
      f11=0
      f12=(Nm12t12x+Nm11t12x+2.*Nm11t22x+Nm12t22x)/(2.*ri)
      f21=(Nm12t12x+Nm12t11x+2.*Nm22t11x+Nm22t12x)/(2.*ri)
      f22=(2.*Nm22t22x+Nm22t12x+Nm12t22x)/(2.*ri)

      r11=((-(f12*f21))**2)/((f11+f12)*(f21+f22)*(f11+f21)*(f12+f22))


      endif


      if(Nm11t12+Nm11t22+Nm12t22.eq.0) then
      f11=(Nm12t12x+Nm11t12x+2.*Nm11t11x+Nm12t11x)/(2.*ri)
      f12=0
      f21=(Nm12t11x+2.*Nm22t11x+Nm22t12x)/(2.*ri)
      f22=(2.*Nm22t22x+Nm22t12x+Nm12t22x+Nm12t12x)/(2.*ri)

      r11=((f11*f22)**2)/((f11+f12)*(f21+f22)*(f11+f21)*(f12+f22))
      endif

      if (Nm12t11+Nm22t11+Nm22t12.eq.0) then        
      f11=(Nm12t12x+Nm11t12x+2.*Nm11t11x+Nm12t11x)/(2.*ri)
      f12=(Nm12t22x+2.*Nm11t22x+Nm11t12x)/(2.*ri)
      f21=0
      f22=(Nm12t12x+2.*Nm22t22x+Nm22t12x+Nm12t22x)/(2.*ri)

      r11=(f11*f22)**2/((f11+f12)*(f21+f22)*(f11+f21)*(f12+f22))


      endif
      
      if(Nm12t22+Nm22t12+Nm22t22.eq.0) then        


      f11=(2.*Nm11t11x+Nm11t12x+Nm12t11x)/(2.*ri)
      f12=(Nm12t12x+Nm11t12x+2.*Nm11t22x+Nm12t22x)/(2.*ri)
      f21=(Nm12t12x+Nm12t11x+2.*Nm22t11x+Nm22t12x)/(2.*ri)
      f22=0

       r11=((-(f12*f21))**2)/((f11+f12)*(f21+f22)*(f11+f21)*(f12+f22))

      endif                  
      f11x=f11
      f12x=f12
      f21x=f21
      f22x=f22      
      
              

      f11=0
      f12=0
      f21=0
      f22=0
      rema=0
      del=0
      
       call em (Nm11t11x, Nm11t12x, Nm11t22x, Nm12t11x, Nm12t12x,
     & Nm12t22x,Nm22t11x, Nm22t12x, Nm22t22x,del,rema,
     &     f11,f12,f21,f22)
       
c       r1=rema
      if(cori(i1,j1).eq.9) rcor(i1,j1)=r11
      rem(i1,j1)=rema
      
       if(pr.eq.1) then
          write(25,*)  "LD estimates:"
      write(25,*)" EM algorithm r2= ", rema
      write(25,*) ""
      if(cori(i1,j1).eq.9) then
      write(25,*)" Corner algorithm r2= ", r11
      write(25,*)""
      write(25,*) "Haplotype Frequencies Corner Algortithm"
      write(25,*) " f11 f12 f21 f22"
      write(25,*) f11x,f12x,f21x,f22x," total sum ", f11x+f12x+f21x+f22x
      write(25,*) ""      
      endif
      write(25,*)""
      write(25,*) "Allele Frequencies"
      write(25,*) ""
      write(25,*) " f1 f2 f1 f2"      
      write(25,*) a1,a2,b1,b2
      write (25,*) ""
      write (25,*) ""      
      write(25,*) "Haplotype Frequencies EM"
      write(25,*) " f11 f12 f21 f22"
      write(25,*) f11,f12,f21,f22," total sum ", f11+f12+f21+f22
      write(25,*) ""

      write(25,*) "****************************************************"
      endif 
cccccccccccccccccccccccccccccccccccccccccccc
 
ccccccccccccccccccccccccc          

 87   continue
 27   continue
      kiu=0
      n1=0

      do i1=1,nsnp(k)-1
         if(corne(i1,i1+1).eq.1.and.n1.eq.0) then
            n1=i1
            go to 987
            endif
         if(corne(i1,i1+1).eq.0.and.n1.gt.0) then
            n2=i1
          kiu=kiu+1
          non1(kiu)=n1
          non2(kiu)=n2          
          
       write (85,190) k, n1,n2, n2-n1,po1(k,n1),po1(k,n2)                     
       n1=0
       n2=0
      endif
 987  continue


      if(cori(i1,i1+1).eq.9) then
      write (89,190) k,i1,i1+1, corne(i1,i1+1),po1(k,i1),po1(k,i1+1)
     %     ,cori(i1,i1+1) ,rem(i1,i1+1),rcor(i1,i1+1)
      dif=dif+abs(rem(i1,i1+1)-rcor(i1,i1+1))
      ndif=ndif+1  
      endif

      
      do j1=i1+1,nsnp(k)
       write (88,190)k, i1,j1, corne(i1,j1),po1(k,i1),po1(k,j1)
 190  format(i3,1x,2I9,2x, i3,I15,2x,2i15,2(2x(f10.8)))
      enddo
      enddo


c      Print*,non1
c      print*," "
c     Print*,non2

      nn=0
      do  nj1=1,nsnp(k)-1
         if(corne(nj1,nj1+1).eq.1) then
            np1=nj1
            np2=np1+1
            go to 41
           endif
        enddo

 41      nflag=0

c          print*,np1,np2         
        do 45 i1=np1,np2
           do 46 j1=i1+1,np2
              if(corne(i1,j1).eq.0) nflag=1
 46        continue
 45     continue
       
       if (nflag.eq.1) then   
         np2=np2-1
         if(np2-np1.ge.1) then         
       write (81,190) k,np1,np2,1-np1+np2
     $      ,po1(k,np1),po1(k,np2)                    
      endif
       np1=np2+1
       np2=np2+1      
       if(np2.ge.nsnp(k)) go to 49        
        go to 41
      else
         
       np2=np2+1
       if(np2.gt.nsnp(k)) then
       write (81,190) k,np1,np2-1,-np1+np2
     $      ,po1(k,np1),po1(k,np2-1)                              
          go to 49
          endif
       go to 41
           endif

 49     continue

       enddo             
       nndif=ndif
       print*, " "
       print*, " "
       Print*, " There were ",nndif,  " consecutive SNP pairs where the 
     % Corners Algorithm was applied"
       print*, " "       
       print*, " "
       Print*,"The absolute value of the difference EM vs Corners was"
     $ , dif/ndif
       print*, " "
       print*, " "
           end
      
