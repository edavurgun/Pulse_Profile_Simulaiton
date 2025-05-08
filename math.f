      SUBROUTINE qtrap(func,a,b,s,eps,par1,par2)
      Implicit None
      INTEGER JMAX
      REAL*8 a,b,func,s,EPS,par1,par2
      EXTERNAL func
      PARAMETER (JMAX=20)
CU    USES trapzd
      INTEGER j
      REAL*8 olds
      include 'common.f'

      olds=-1.e30
      do 11 j=1,JMAX
        call trapzd(func,a,b,s,j,par1,par2)
        if (abs(s-olds).lt.EPS*abs(olds)) return
        if (s.eq.0..and.olds.eq.0..and.j.gt.6) return
        olds=s
11    continue
c      pause 'too many steps in qtrap'
      END

      SUBROUTINE trapzd(func,a,b,s,n,par1,par2)
      Implicit None
      INTEGER n
      REAL*8 a,b,s,func,par1,par2
      EXTERNAL func
      INTEGER it,j
      REAL*8 del,sum,tnm,x
      include 'common.f'

      if (n.eq.1) then
        s=0.5*(b-a)*(func(a,par1,par2)+func(b,par1,par2))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x,par1,par2)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END













