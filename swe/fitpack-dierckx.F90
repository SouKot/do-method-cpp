!* ============================================================================ *
!*
!* Multi-Layer Shallow Water Model
!*
!* Version 1.4 - 20/05/2008 - ATvS
!*
!* ============================================================================ *
      subroutine spgrid(iopt,ider,mu,u,mv,v,r,r0,r1,s,nuest,nvest,nu,tu,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
!c  given the function values r(i,j) on the latitude-longitude grid
!c  (u(i),v(j)), i=1,...,mu ; j=1,...,mv , spgrid determines a smooth
!c  bicubic spline approximation on the rectangular domain 0<=u<=pi,
!c  vb<=v<=ve (vb = v(1), ve=vb+2*pi).
!c  this approximation s(u,v) will satisfy the properties
!c
!c    (1) s(0,v) = s(0,0) = dr(1)
!c
!c        d s(0,v)           d s(0,0)           d s(0,pi/2)
!c    (2) -------- = cos(v)* -------- + sin(v)* -----------
!c        d u                d u                d u
!c
!c                 = cos(v)*dr(2)+sin(v)*dr(3)
!c                                                     vb <= v <= ve
!c    (3) s(pi,v) = s(pi,0) = dr(4)
!c
!c        d s(pi,v)           d s(pi,0)           d s(pi,pi/2)
!c    (4) -------- = cos(v)*  --------- + sin(v)* ------------
!c        d u                 d u                 d u
!c
!c                 = cos(v)*dr(5)+sin(v)*dr(6)
!c
!c  and will be periodic in the variable v, i.e.
!c
!c         j           j
!c        d s(u,vb)   d s(u,ve)
!c    (5) --------- = ---------   0 <=u<= pi , j=0,1,2
!c           j           j
!c        d v         d v
!c
!c  the number of knots of s(u,v) and their position tu(i),i=1,2,...,nu;
!c  tv(j),j=1,2,...,nv, is chosen automatically by the routine. the
!c  smoothness of s(u,v) is achieved by minimalizing the discontinuity
!c  jumps of the derivatives of the spline at the knots. the amount of
!c  smoothness of s(u,v) is determined by the condition that
!c  fp=sumi=1,mu(sumj=1,mv((r(i,j)-s(u(i),v(j)))**2))+(r0-s(0,v))**2
!c  + (r1-s(pi,v))**2 <= s, with s a given non-negative constant.
!c  the fit s(u,v) is given in its b-spline representation and can be
!c  evaluated by means of routine bispev
!c
!c calling sequence:
!c     call spgrid(iopt,ider,mu,u,mv,v,r,r0,r1,s,nuest,nvest,nu,tu,
!c    *  ,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
!c
!c parameters:
!c  iopt  : integer array of dimension 3, specifying different options.
!c          unchanged on exit.
!c  iopt(1):on entry iopt(1) must specify whether a least-squares spline
!c          (iopt(1)=-1) or a smoothing spline (iopt(1)=0 or 1) must be
!c          determined.
!c          if iopt(1)=0 the routine will start with an initial set of
!c          knots tu(i)=0,tu(i+4)=pi,i=1,...,4;tv(i)=v(1)+(i-4)*2*pi,
!c          i=1,...,8.
!c          if iopt(1)=1 the routine will continue with the set of knots
!c          found at the last call of the routine.
!c          attention: a call with iopt(1)=1 must always be immediately
!c          preceded by another call with iopt(1) = 1 or iopt(1) = 0.
!c  iopt(2):on entry iopt(2) must specify the requested order of conti-
!c          nuity at the pole u=0.
!c          if iopt(2)=0 only condition (1) must be fulfilled and
!c          if iopt(2)=1 conditions (1)+(2) must be fulfilled.
!c  iopt(3):on entry iopt(3) must specify the requested order of conti-
!c          nuity at the pole u=pi.
!c          if iopt(3)=0 only condition (3) must be fulfilled and
!c          if iopt(3)=1 conditions (3)+(4) must be fulfilled.
!c  ider  : integer array of dimension 4, specifying different options.
!c          unchanged on exit.
!c  ider(1):on entry ider(1) must specify whether (ider(1)=0 or 1) or not
!c          (ider(1)=-1) there is a data value r0 at the pole u=0.
!c          if ider(1)=1, r0 will be considered to be the right function
!c          value, and it will be fitted exactly (s(0,v)=r0).
!c          if ider(1)=0, r0 will be considered to be a data value just
!c          like the other data values r(i,j).
!c  ider(2):on entry ider(2) must specify whether (ider(2)=1) or not
!c          (ider(2)=0) the approximation has vanishing derivatives
!c          dr(2) and dr(3) at the pole u=0  (in case iopt(2)=1)
!c  ider(3):on entry ider(3) must specify whether (ider(3)=0 or 1) or not
!c          (ider(3)=-1) there is a data value r1 at the pole u=pi.
!c          if ider(3)=1, r1 will be considered to be the right function
!c          value, and it will be fitted exactly (s(pi,v)=r1).
!c          if ider(3)=0, r1 will be considered to be a data value just
!c          like the other data values r(i,j).
!c  ider(4):on entry ider(4) must specify whether (ider(4)=1) or not
!c          (ider(4)=0) the approximation has vanishing derivatives
!c          dr(5) and dr(6) at the pole u=pi (in case iopt(3)=1)
!c  mu    : integer. on entry mu must specify the number of grid points
!c          along the u-axis. unchanged on exit.
!c          mu >= 1, mu >=mumin=4-i0-i1-ider(2)-ider(4) with
!c            i0=min(1,ider(1)+1), i1=min(1,ider(3)+1)
!c  u     : real array of dimension at least (mu). before entry, u(i)
!c          must be set to the u-co-ordinate of the i-th grid point
!c          along the u-axis, for i=1,2,...,mu. these values must be
!c          supplied in strictly ascending order. unchanged on exit.
!c          0 < u(i) < pi.
!c  mv    : integer. on entry mv must specify the number of grid points
!c          along the v-axis. mv > 3 . unchanged on exit.
!c  v     : real array of dimension at least (mv). before entry, v(j)
!c          must be set to the v-co-ordinate of the j-th grid point
!c          along the v-axis, for j=1,2,...,mv. these values must be
!c          supplied in strictly ascending order. unchanged on exit.
!c          -pi <= v(1) < pi , v(mv) < v(1)+2*pi.
!c  r     : real array of dimension at least (mu*mv).
!c          before entry, r(mv*(i-1)+j) must be set to the data value at
!c          the grid point (u(i),v(j)) for i=1,...,mu and j=1,...,mv.
!c          unchanged on exit.
!c  r0    : real value. on entry (if ider(1) >=0 ) r0 must specify the
!c          data value at the pole u=0. unchanged on exit.
!c  r1    : real value. on entry (if ider(1) >=0 ) r1 must specify the
!c          data value at the pole u=pi. unchanged on exit.
!c  s     : real. on entry (if iopt(1)>=0) s must specify the smoothing
!c          factor. s >=0. unchanged on exit.
!c          for advice on the choice of s see further comments
!c  nuest : integer. unchanged on exit.
!c  nvest : integer. unchanged on exit.
!c          on entry, nuest and nvest must specify an upper bound for the
!c          number of knots required in the u- and v-directions respect.
!c          these numbers will also determine the storage space needed by
!c          the routine. nuest >= 8, nvest >= 8.
!c          in most practical situation nuest = mu/2, nvest=mv/2, will
!c          be sufficient. always large enough are nuest=mu+6+iopt(2)+
!c          iopt(3), nvest = mv+7, the number of knots needed for
!c          interpolation (s=0). see also further comments.
!c  nu    : integer.
!c          unless ier=10 (in case iopt(1)>=0), nu will contain the total
!c          number of knots with respect to the u-variable, of the spline
!c          approximation returned. if the computation mode iopt(1)=1 is
!c          used, the value of nu should be left unchanged between sub-
!c          sequent calls. in case iopt(1)=-1, the value of nu should be
!c          specified on entry.
!c  tu    : real array of dimension at least (nuest).
!c          on succesful exit, this array will contain the knots of the
!c          spline with respect to the u-variable, i.e. the position of
!c          the interior knots tu(5),...,tu(nu-4) as well as the position
!c          of the additional knots tu(1)=...=tu(4)=0 and tu(nu-3)=...=
!c          tu(nu)=pi needed for the b-spline representation.
!c          if the computation mode iopt(1)=1 is used,the values of tu(1)
!c          ...,tu(nu) should be left unchanged between subsequent calls.
!c          if the computation mode iopt(1)=-1 is used, the values tu(5),
!c          ...tu(nu-4) must be supplied by the user, before entry.
!c          see also the restrictions (ier=10).
!c  nv    : integer.
!c          unless ier=10 (in case iopt(1)>=0), nv will contain the total
!c          number of knots with respect to the v-variable, of the spline
!c          approximation returned. if the computation mode iopt(1)=1 is
!c          used, the value of nv should be left unchanged between sub-
!c          sequent calls. in case iopt(1) = -1, the value of nv should
!c          be specified on entry.
!c  tv    : real array of dimension at least (nvest).
!c          on succesful exit, this array will contain the knots of the
!c          spline with respect to the v-variable, i.e. the position of
!c          the interior knots tv(5),...,tv(nv-4) as well as the position
!c          of the additional knots tv(1),...,tv(4) and tv(nv-3),...,
!c          tv(nv) needed for the b-spline representation.
!c          if the computation mode iopt(1)=1 is used,the values of tv(1)
!c          ...,tv(nv) should be left unchanged between subsequent calls.
!c          if the computation mode iopt(1)=-1 is used, the values tv(5),
!c          ...tv(nv-4) must be supplied by the user, before entry.
!c          see also the restrictions (ier=10).
!c  c     : real array of dimension at least (nuest-4)*(nvest-4).
!c          on succesful exit, c contains the coefficients of the spline
!c          approximation s(u,v)
!c  fp    : real. unless ier=10, fp contains the sum of squared
!c          residuals of the spline approximation returned.
!c  wrk   : real array of dimension (lwrk). used as workspace.
!c          if the computation mode iopt(1)=1 is used the values of
!c          wrk(1),..,wrk(12) should be left unchanged between subsequent
!c          calls.
!c  lwrk  : integer. on entry lwrk must specify the actual dimension of
!c          the array wrk as declared in the calling (sub)program.
!c          lwrk must not be too small.
!c           lwrk >= 12+nuest*(mv+nvest+3)+nvest*24+4*mu+8*mv+q
!c           where q is the larger of (mv+nvest) and nuest.
!c  iwrk  : integer array of dimension (kwrk). used as workspace.
!c          if the computation mode iopt(1)=1 is used the values of
!c          iwrk(1),.,iwrk(5) should be left unchanged between subsequent
!c          calls.
!c  kwrk  : integer. on entry kwrk must specify the actual dimension of
!c          the array iwrk as declared in the calling (sub)program.
!c          kwrk >= 5+mu+mv+nuest+nvest.
!c  ier   : integer. unless the routine detects an error, ier contains a
!c          non-positive value on exit, i.e.
!c   ier=0  : normal return. the spline returned has a residual sum of
!c            squares fp such that abs(fp-s)/s <= tol with tol a relat-
!c            ive tolerance set to 0.001 by the program.
!c   ier=-1 : normal return. the spline returned is an interpolating
!c            spline (fp=0).
!c   ier=-2 : normal return. the spline returned is the least-squares
!c            constrained polynomial. in this extreme case fp gives the
!c            upper bound for the smoothing factor s.
!c   ier=1  : error. the required storage space exceeds the available
!c            storage space, as specified by the parameters nuest and
!c            nvest.
!c            probably causes : nuest or nvest too small. if these param-
!c            eters are already large, it may also indicate that s is
!c            too small
!c            the approximation returned is the least-squares spline
!c            according to the current set of knots. the parameter fp
!c            gives the corresponding sum of squared residuals (fp>s).
!c   ier=2  : error. a theoretically impossible result was found during
!c            the iteration proces for finding a smoothing spline with
!c            fp = s. probably causes : s too small.
!c            there is an approximation returned but the corresponding
!c            sum of squared residuals does not satisfy the condition
!c            abs(fp-s)/s < tol.
!c   ier=3  : error. the maximal number of iterations maxit (set to 20
!c            by the program) allowed for finding a smoothing spline
!c            with fp=s has been reached. probably causes : s too small
!c            there is an approximation returned but the corresponding
!c            sum of squared residuals does not satisfy the condition
!c            abs(fp-s)/s < tol.
!c   ier=10 : error. on entry, the input data are controlled on validity
!c            the following restrictions must be satisfied.
!c            -1<=iopt(1)<=1, 0<=iopt(2)<=1, 0<=iopt(3)<=1,
!c            -1<=ider(1)<=1, 0<=ider(2)<=1, ider(2)=0 if iopt(2)=0.
!c            -1<=ider(3)<=1, 0<=ider(4)<=1, ider(4)=0 if iopt(3)=0.
!c            mu >= mumin (see above), mv >= 4, nuest >=8, nvest >= 8,
!c            kwrk>=5+mu+mv+nuest+nvest,
!c            lwrk >= 12+nuest*(mv+nvest+3)+nvest*24+4*mu+8*mv+
!c             max(nuest,mv+nvest)
!c            0< u(i-1)<u(i)< pi,i=2,..,mu,
!c            -pi<=v(1)< pi, v(1)<v(i-1)<v(i)<v(1)+2*pi, i=3,...,mv
!c            if iopt(1)=-1: 8<=nu<=min(nuest,mu+6+iopt(2)+iopt(3))
!c                           0<tu(5)<tu(6)<...<tu(nu-4)< pi
!c                           8<=nv<=min(nvest,mv+7)
!c                           v(1)<tv(5)<tv(6)<...<tv(nv-4)<v(1)+2*pi
!c                    the schoenberg-whitney conditions, i.e. there must
!c                    be subset of grid co-ordinates uu(p) and vv(q) such
!c                    that   tu(p) < uu(p) < tu(p+4) ,p=1,...,nu-4
!c                     (iopt(2)=1 and iopt(3)=1 also count for a uu-value
!c                           tv(q) < vv(q) < tv(q+4) ,q=1,...,nv-4
!c                     (vv(q) is either a value v(j) or v(j)+2*pi)
!c            if iopt(1)>=0: s>=0
!c                       if s=0: nuest>=mu+6+iopt(2)+iopt(3), nvest>=mv+7
!c            if one of these conditions is found to be violated,control
!c            is immediately repassed to the calling program. in that
!c            case there is no approximation returned.
!c
!c further comments:
!c   spgrid does not allow individual weighting of the data-values.
!c   so, if these were determined to widely different accuracies, then
!c   perhaps the general data set routine sphere should rather be used
!c   in spite of efficiency.
!c   by means of the parameter s, the user can control the tradeoff
!c   between closeness of fit and smoothness of fit of the approximation.
!c   if s is too large, the spline will be too smooth and signal will be
!c   lost ; if s is too small the spline will pick up too much noise. in
!c   the extreme cases the program will return an interpolating spline if
!c   s=0 and the constrained least-squares polynomial(degrees 3,0)if s is
!c   very large. between these extremes, a properly chosen s will result
!c   in a good compromise between closeness of fit and smoothness of fit.
!c   to decide whether an approximation, corresponding to a certain s is
!c   satisfactory the user is highly recommended to inspect the fits
!c   graphically.
!c   recommended values for s depend on the accuracy of the data values.
!c   if the user has an idea of the statistical errors on the data, he
!c   can also find a proper estimate for s. for, by assuming that, if he
!c   specifies the right s, spgrid will return a spline s(u,v) which
!c   exactly reproduces the function underlying the data he can evaluate
!c   the sum((r(i,j)-s(u(i),v(j)))**2) to find a good estimate for this s
!c   for example, if he knows that the statistical errors on his r(i,j)-
!c   values is not greater than 0.1, he may expect that a good s should
!c   have a value not larger than mu*mv*(0.1)**2.
!c   if nothing is known about the statistical error in r(i,j), s must
!c   be determined by trial and error, taking account of the comments
!c   above. the best is then to start with a very large value of s (to
!c   determine the least-squares polynomial and the corresponding upper
!c   bound fp0 for s) and then to progressively decrease the value of s
!c   ( say by a factor 10 in the beginning, i.e. s=fp0/10,fp0/100,...
!c   and more carefully as the approximation shows more detail) to
!c   obtain closer fits.
!c   to economize the search for a good s-value the program provides with
!c   different modes of computation. at the first call of the routine, or
!c   whenever he wants to restart with the initial set of knots the user
!c   must set iopt(1)=0.
!c   if iopt(1) = 1 the program will continue with the knots found at
!c   the last call of the routine. this will save a lot of computation
!c   time if spgrid is called repeatedly for different values of s.
!c   the number of knots of the spline returned and their location will
!c   depend on the value of s and on the complexity of the shape of the
!c   function underlying the data. if the computation mode iopt(1) = 1
!c   is used, the knots returned may also depend on the s-values at
!c   previous calls (if these were smaller). therefore, if after a number
!c   of trials with different s-values and iopt(1)=1,the user can finally
!c   accept a fit as satisfactory, it may be worthwhile for him to call
!c   spgrid once more with the chosen value for s but now with iopt(1)=0.
!c   indeed, spgrid may then return an approximation of the same quality
!c   of fit but with fewer knots and therefore better if data reduction
!c   is also an important objective for the user.
!c   the number of knots may also depend on the upper bounds nuest and
!c   nvest. indeed, if at a certain stage in spgrid the number of knots
!c   in one direction (say nu) has reached the value of its upper bound
!c   (nuest), then from that moment on all subsequent knots are added
!c   in the other (v) direction. this may indicate that the value of
!c   nuest is too small. on the other hand, it gives the user the option
!c   of limiting the number of knots the routine locates in any direction
!c   for example, by setting nuest=8 (the lowest allowable value for
!c   nuest), the user can indicate that he wants an approximation which
!c   is a simple cubic polynomial in the variable u.
!c
!c  other subroutines required:
!c    fpspgr,fpchec,fpchep,fpknot,fpopsp,fprati,fpgrsp,fpsysy,fpback,
!c    fpbacp,fpbspl,fpcyt1,fpcyt2,fpdisc,fpgivs,fprota
!c
!c  references:
!c   dierckx p. : fast algorithms for smoothing data over a disc or a
!c                sphere using tensor product splines, in "algorithms
!c                for approximation", ed. j.c.mason and m.g.cox,
!c                clarendon press oxford, 1987, pp. 51-65
!c   dierckx p. : fast algorithms for smoothing data over a disc or a
!c                sphere using tensor product splines, report tw73, dept.
!c                computer science,k.u.leuven, 1985.
!c   dierckx p. : curve and surface fitting with splines, monographs on
!c                numerical analysis, oxford university press, 1993.
!c
!c  author:
!c    p.dierckx
!c    dept. computer science, k.u. leuven
!c    celestijnenlaan 200a, b-3001 heverlee, belgium.
!c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!c
!c  creation date : july 1985
!c  latest update : march 1989
!c
!c  ..
!c  ..scalar arguments..
      real r0,r1,s,fp
      integer mu,mv,nuest,nvest,nu,nv,lwrk,kwrk,ier
!c  ..array arguments..
      integer iopt(3),ider(4),iwrk(kwrk)
      real u(mu),v(mv),r(mu*mv),c((nuest-4)*(nvest-4)),tu(nuest), tv(nvest),wrk(lwrk)
!c  ..local scalars..
      real per,pi,tol,uu,ve,rmax,rmin,one,half,rn,rb,re
      integer i,i1,i2,j,jwrk,j1,j2,kndu,kndv,knru,knrv,kwest,l,ldr,lfpu,lfpv,lwest,lww,m,maxit,mumin,muu,nc
!c  ..function references..
      real atan2
      integer max0
!c  ..subroutine references..
!c    fpchec,fpchep,fpspgr
!c  ..
!c  set constants
      one = 1
      half = 0.5e0
      pi = atan2(0.,-one)
      per = pi+pi
      ve = v(1)+per
!c  we set up the parameters tol and maxit.
      maxit = 20
      tol = 0.1e-02
!c  before starting computations, a data check is made. if the input data
!c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(iopt(1).lt.(-1) .or. iopt(1).gt.1) go to 200
      if(iopt(2).lt.0 .or. iopt(2).gt.1) go to 200
      if(iopt(3).lt.0 .or. iopt(3).gt.1) go to 200
      if(ider(1).lt.(-1) .or. ider(1).gt.1) go to 200
      if(ider(2).lt.0 .or. ider(2).gt.1) go to 200
      if(ider(2).eq.1 .and. iopt(2).eq.0) go to 200
      if(ider(3).lt.(-1) .or. ider(3).gt.1) go to 200
      if(ider(4).lt.0 .or. ider(4).gt.1) go to 200
      if(ider(4).eq.1 .and. iopt(3).eq.0) go to 200
      mumin = 4
      if(ider(1).ge.0) mumin = mumin-1
      if(iopt(2).eq.1 .and. ider(2).eq.1) mumin = mumin-1
      if(ider(3).ge.0) mumin = mumin-1
      if(iopt(3).eq.1 .and. ider(4).eq.1) mumin = mumin-1
      if(mumin.eq.0) mumin = 1
      if(mu.lt.mumin .or. mv.lt.4) go to 200
      if(nuest.lt.8 .or. nvest.lt.8) go to 200
      m = mu*mv
      nc = (nuest-4)*(nvest-4)
      lwest = 12+nuest*(mv+nvest+3)+24*nvest+4*mu+8*mv+max0(nuest,mv+nvest)
      kwest = 5+mu+mv+nuest+nvest
      if(lwrk.lt.lwest .or. kwrk.lt.kwest) go to 200
      if(u(1).le.0. .or. u(mu).ge.pi) go to 200
      if(mu.eq.1) go to 30
      do 20 i=2,mu
        if(u(i-1).ge.u(i)) go to 200
  20  continue
  30  if(v(1).lt. (-pi) .or. v(1).ge.pi ) go to 200
      if(v(mv).ge.v(1)+per) go to 200
      do 40 i=2,mv
        if(v(i-1).ge.v(i)) go to 200
  40  continue
      if(iopt(1).gt.0) go to 140
!c  if not given, we compute an estimate for r0.
      rn = mv
      if(ider(1).lt.0) go to 45
      rb = r0
      go to 55
  45  rb = 0.
      do 50 i=1,mv
         rb = rb+r(i)
  50  continue
      rb = rb/rn
!c  if not given, we compute an estimate for r1.
  55  if(ider(3).lt.0) go to 60
      re = r1
      go to 70
  60  re = 0.
      j = m
      do 65 i=1,mv
         re = re+r(j)
         j = j-1
  65  continue
      re = re/rn
!c  we determine the range of r-values.
  70  rmin = rb
      rmax = re
      do 80 i=1,m
         if(r(i).lt.rmin) rmin = r(i)
         if(r(i).gt.rmax) rmax = r(i)
  80  continue
      wrk(5) = rb
      wrk(6) = 0.
      wrk(7) = 0.
      wrk(8) = re
      wrk(9) = 0.
      wrk(10) = 0.
      wrk(11) = rmax -rmin
      wrk(12) = wrk(11)
      iwrk(4) = mu
      iwrk(5) = mu
      if(iopt(1).eq.0) go to 140
      if(nu.lt.8 .or. nu.gt.nuest) go to 200
      if(nv.lt.11 .or. nv.gt.nvest) go to 200
      j = nu
      do 90 i=1,4
        tu(i) = 0.
        tu(j) = pi
        j = j-1
  90  continue
      l = 13
      wrk(l) = 0.
      if(iopt(2).eq.0) go to 100
      l = l+1
      uu = u(1)
      if(uu.gt.tu(5)) uu = tu(5)
      wrk(l) = uu*half
 100  do 110 i=1,mu
        l = l+1
        wrk(l) = u(i)
 110  continue
      if(iopt(3).eq.0) go to 120
      l = l+1
      uu = u(mu)
      if(uu.lt.tu(nu-4)) uu = tu(nu-4)
      wrk(l) = uu+(pi-uu)*half
 120  l = l+1
      wrk(l) = pi
      muu = l-12
      call fpchec(wrk(13),muu,tu,nu,3,ier)
      if(ier.ne.0) go to 200
      j1 = 4
      tv(j1) = v(1)
      i1 = nv-3
      tv(i1) = ve
      j2 = j1
      i2 = i1
      do 130 i=1,3
        i1 = i1+1
        i2 = i2-1
        j1 = j1+1
        j2 = j2-1
        tv(j2) = tv(i2)-per
        tv(i1) = tv(j1)+per
 130  continue
      l = 13
      do 135 i=1,mv
        wrk(l) = v(i)
        l = l+1
 135  continue
      wrk(l) = ve
      call fpchep(wrk(13),mv+1,tv,nv,3,ier)
      if(ier) 200,150,200
 140  if(s.lt.0.) go to 200
      if(s.eq.0. .and. (nuest.lt.(mu+6+iopt(2)+iopt(3)) .or. nvest.lt.(mv+7)) ) go to 200
!c  we partition the working space and determine the spline approximation
 150  ldr = 5
      lfpu = 13
      lfpv = lfpu+nuest
      lww = lfpv+nvest
      jwrk = lwrk-12-nuest-nvest
      knru = 6
      knrv = knru+mu
      kndu = knrv+mv
      kndv = kndu+nuest
      call fpspgr(iopt,ider,u,mu,v,mv,r,m,rb,re,s,nuest,nvest,tol,maxit,&
      nc,nu,tu,nv,tv,c,fp,wrk(1),wrk(2),wrk(3),wrk(4),wrk(lfpu),&
      wrk(lfpv),wrk(ldr),wrk(11),iwrk(1),iwrk(2),iwrk(3),iwrk(4),&
      iwrk(5),iwrk(knru),iwrk(knrv),iwrk(kndu),iwrk(kndv),wrk(lww),&
      jwrk,ier)
 200  return
      end
      subroutine bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,iwrk,kwrk,ier)
!c  subroutine bispev evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,...
!c  ,my a bivariate spline s(x,y) of degrees kx and ky, given in the
!c  b-spline representation.
!c
!c  calling sequence:
!c     call bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,
!c    * iwrk,kwrk,ier)
!c
!c  input parameters:
!c   tx    : real array, length nx, which contains the position of the
!c           knots in the x-direction.
!c   nx    : integer, giving the total number of knots in the x-direction
!c   ty    : real array, length ny, which contains the position of the
!c           knots in the y-direction.
!c   ny    : integer, giving the total number of knots in the y-direction
!c   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
!c           b-spline coefficients.
!c   kx,ky : integer values, giving the degrees of the spline.
!c   x     : real array of dimension (mx).
!c           before entry x(i) must be set to the x co-ordinate of the
!c           i-th grid point along the x-axis.
!c           tx(kx+1)<=x(i-1)<=x(i)<=tx(nx-kx), i=2,...,mx.
!c   mx    : on entry mx must specify the number of grid points along
!c           the x-axis. mx >=1.
!c   y     : real array of dimension (my).
!c           before entry y(j) must be set to the y co-ordinate of the
!c           j-th grid point along the y-axis.
!c           ty(ky+1)<=y(j-1)<=y(j)<=ty(ny-ky), j=2,...,my.
!c   my    : on entry my must specify the number of grid points along
!c           the y-axis. my >=1.
!c   wrk   : real array of dimension lwrk. used as workspace.
!c   lwrk  : integer, specifying the dimension of wrk.
!c           lwrk >= mx*(kx+1)+my*(ky+1)
!c   iwrk  : integer array of dimension kwrk. used as workspace.
!c   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my.
!c
!c  output parameters:
!c   z     : real array of dimension (mx*my).
!c           on succesful exit z(my*(i-1)+j) contains the value of s(x,y)
!c           at the point (x(i),y(j)),i=1,...,mx;j=1,...,my.
!c   ier   : integer error flag
!c    ier=0 : normal return
!c    ier=10: invalid input data (see restrictions)
!c
!c  restrictions:
!c   mx >=1, my >=1, lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my
!c   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
!c   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
!c
!c  other subroutines required:
!c    fpbisp,fpbspl
!c
!c  references :
!c    de boor c : on calculating with b-splines, j. approximation theory
!c                6 (1972) 50-62.
!c    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
!c                applics 10 (1972) 134-149.
!c    dierckx p. : curve and surface fitting with splines, monographs on
!c                 numerical analysis, oxford university press, 1993.
!c
!c  author :
!c    p.dierckx
!c    dept. computer science, k.u.leuven
!c    celestijnenlaan 200a, b-3001 heverlee, belgium.
!c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!c
!c  latest update : march 1987
!c
!c  ..scalar arguments..
      integer nx,ny,kx,ky,mx,my,lwrk,kwrk,ier
!c  ..array arguments..
      integer iwrk(kwrk)
      real tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my), wrk(lwrk)
!c  ..local scalars..
      integer i,iw,lwest
!c  ..
!c  before starting computations a data check is made. if the input data
!c  are invalid control is immediately repassed to the calling program.
      ier = 10
      lwest = (kx+1)*mx+(ky+1)*my
      if(lwrk.lt.lwest) go to 100
      if(kwrk.lt.(mx+my)) go to 100
      if(mx-1) 100,30,10
  10  do 20 i=2,mx
        if(x(i).lt.x(i-1)) go to 100
  20  continue
  30  if(my-1) 100,60,40
  40  do 50 i=2,my
        if(y(i).lt.y(i-1)) go to 100
  50  continue
  60  ier = 0
      iw = mx*(kx+1)+1
      call fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk(1),wrk(iw), iwrk(1),iwrk(mx+1))
 100  return
      end
      subroutine fpback(a,z,n,k,c,nest)
!c  subroutine fpback calculates the solution of the system of
!c  equations a*c = z with a a n x n upper triangular matrix
!c  of bandwidth k.
!c  ..
!c  ..scalar arguments..
      integer n,k,nest
!c  ..array arguments..
      real a(nest,k),z(n),c(n)
!c  ..local scalars..
      real store
      integer i,i1,j,k1,l,m
!c  ..
      k1 = k-1
      c(n) = z(n)/a(n,1)
      i = n-1
      if(i.eq.0) go to 30
      do 20 j=2,n
        store = z(i)
        i1 = k1
        if(j.le.k1) i1 = j-1
        m = i
        do 10 l=1,i1
          m = m+1
          store = store-c(m)*a(i,l+1)
  10    continue
        c(i) = store/a(i,1)
        i = i-1
  20  continue
  30  return
      end
      subroutine fpbacp(a,b,z,n,k,c,k1,nest)
!c  subroutine fpbacp calculates the solution of the system of equations
!c  g * c = z  with g  a n x n upper triangular matrix of the form
!c            ! a '   !
!c        g = !   ' b !
!c            ! 0 '   !
!c  with b a n x k matrix and a a (n-k) x (n-k) upper triangular
!c  matrix of bandwidth k1.
!c  ..
!c  ..scalar arguments..
      integer n,k,k1,nest
!c  ..array arguments..
      real a(nest,k1),b(nest,k),z(n),c(n)
!c  ..local scalars..
      integer i,i1,j,l,l0,l1,n2
      real store
!c  ..
      n2 = n-k
      l = n
      do 30 i=1,k
        store = z(l)
        j = k+2-i
        if(i.eq.1) go to 20
        l0 = l
        do 10 l1=j,k
          l0 = l0+1
          store = store-c(l0)*b(l,l1)
  10    continue
  20    c(l) = store/b(l,j-1)
        l = l-1
        if(l.eq.0) go to 80
  30  continue
      do 50 i=1,n2
        store = z(i)
        l = n2
        do 40 j=1,k
          l = l+1
          store = store-c(l)*b(i,j)
  40    continue
        c(i) = store
  50  continue
      i = n2
      c(i) = c(i)/a(i,1)
      if(i.eq.1) go to 80
      do 70 j=2,n2
        i = i-1
        store = c(i)
        i1 = k
        if(j.le.k) i1=j-1
        l = i
        do 60 l0=1,i1
          l = l+1
          store = store-c(l)*a(i,l0+1)
  60    continue
        c(i) = store/a(i,1)
  70  continue
  80  return
      end
      subroutine fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wx,wy,lx,ly)
!c  ..scalar arguments..
      integer nx,ny,kx,ky,mx,my
!c  ..array arguments..
      integer lx(mx),ly(my)
      real tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my), wx(mx,kx+1),wy(my,ky+1)
!c  ..local scalars..
      integer kx1,ky1,l,l1,l2,m,nkx1,nky1
      real arg,sp,tb,te
!c  ..local arrays..
      real h(6)
!c  ..subroutine references..
!c    fpbspl
!c  ..
      kx1 = kx+1
      nkx1 = nx-kx1
      tb = tx(kx1)
      te = tx(nkx1+1)
      l = kx1
      l1 = l+1
      do 40 i=1,mx
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
  10    if(arg.lt.tx(l1) .or. l.eq.nkx1) go to 20
        l = l1
        l1 = l+1
        go to 10
  20    call fpbspl(tx,nx,kx,arg,l,h)
        lx(i) = l-kx1
        do 30 j=1,kx1
          wx(i,j) = h(j)
  30    continue
  40  continue
      ky1 = ky+1
      nky1 = ny-ky1
      tb = ty(ky1)
      te = ty(nky1+1)
      l = ky1
      l1 = l+1
      do 80 i=1,my
        arg = y(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
  50    if(arg.lt.ty(l1) .or. l.eq.nky1) go to 60
        l = l1
        l1 = l+1
        go to 50
  60    call fpbspl(ty,ny,ky,arg,l,h)
        ly(i) = l-ky1
        do 70 j=1,ky1
          wy(i,j) = h(j)
  70    continue
  80  continue
      m = 0
      do 130 i=1,mx
        l = lx(i)*nky1
        do 90 i1=1,kx1
          h(i1) = wx(i,i1)
  90    continue
        do 120 j=1,my
          l1 = l+ly(j)
          sp = 0.
          do 110 i1=1,kx1
            l2 = l1
            do 100 j1=1,ky1
              l2 = l2+1
              sp = sp+c(l2)*h(i1)*wy(j,j1)
 100        continue
            l1 = l1+nky1
 110      continue
          m = m+1
          z(m) = sp
 120    continue
 130  continue
      return
      end

      subroutine fpbspl(t,n,k,x,l,h)
!c  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
!c  degree k at t(l) <= x < t(l+1) using the stable recurrence
!c  relation of de boor and cox.
!c  ..
!c  ..scalar arguments..
      real x
      integer n,k,l
!c  ..array arguments..
      real t(n),h(6)
!c  ..local scalars..
      real f,one
      integer i,j,li,lj
!c  ..local arrays..
      real hh(5)
!c  ..
      one = 0.1e+01
      h(1) = one
      do 20 j=1,k
        do 10 i=1,j
          hh(i) = h(i)
  10    continue
        h(1) = 0.
        do 20 i=1,j
          li = l+i
          lj = li-j
          f = hh(i)/(t(li)-t(lj))
          h(i) = h(i)+f*(t(li)-x)
          h(i+1) = f*(x-t(lj))
  20  continue
      return
      end
      subroutine fpchec(x,m,t,n,k,ier)
!c  subroutine fpchec verifies the number and the position of the knots
!c  t(j),j=1,2,...,n of a spline of degree k, in relation to the number
!c  and the position of the data points x(i),i=1,2,...,m. if all of the
!c  following conditions are fulfilled, the error parameter ier is set
!c  to zero. if one of the conditions is violated ier is set to ten.
!c      1) k+1 <= n-k-1 <= m
!c      2) t(1) <= t(2) <= ... <= t(k+1)
!c         t(n-k) <= t(n-k+1) <= ... <= t(n)
!c      3) t(k+1) < t(k+2) < ... < t(n-k)
!c      4) t(k+1) <= x(i) <= t(n-k)
!c      5) the conditions specified by schoenberg and whitney must hold
!c         for at least one subset of data points, i.e. there must be a
!c         subset of data points y(j) such that
!c             t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1
!c  ..
!c  ..scalar arguments..
      integer m,n,k,ier
!c  ..array arguments..
      real x(m),t(n)
!c  ..local scalars..
      integer i,j,k1,k2,l,nk1,nk2,nk3
      real tj,tl
!c  ..
      k1 = k+1
      k2 = k1+1
      nk1 = n-k1
      nk2 = nk1+1
      ier = 10
!c  check condition no 1
      if(nk1.lt.k1 .or. nk1.gt.m) go to 80
!c  check condition no 2
      j = n
      do 20 i=1,k
        if(t(i).gt.t(i+1)) go to 80
        if(t(j).lt.t(j-1)) go to 80
        j = j-1
  20  continue
!c  check condition no 3
      do 30 i=k2,nk2
        if(t(i).le.t(i-1)) go to 80
  30  continue
!c  check condition no 4
      if(x(1).lt.t(k1) .or. x(m).gt.t(nk2)) go to 80
!c  check condition no 5
      if(x(1).ge.t(k2) .or. x(m).le.t(nk1)) go to 80
      i = 1
      l = k2
      nk3 = nk1-1
      if(nk3.lt.2) go to 70
      do 60 j=2,nk3
        tj = t(j)
        l = l+1
        tl = t(l)
  40    i = i+1
        if(i.ge.m) go to 80
        if(x(i).le.tj) go to 40
        if(x(i).ge.tl) go to 80
  60  continue
  70  ier = 0
  80  return
      end
      subroutine fpchep(x,m,t,n,k,ier)
!c  subroutine fpchep verifies the number and the position of the knots
!c  t(j),j=1,2,...,n of a periodic spline of degree k, in relation to
!c  the number and the position of the data points x(i),i=1,2,...,m.
!c  if all of the following conditions are fulfilled, ier is set
!c  to zero. if one of the conditions is violated ier is set to ten.
!c      1) k+1 <= n-k-1 <= m+k-1
!c      2) t(1) <= t(2) <= ... <= t(k+1)
!c         t(n-k) <= t(n-k+1) <= ... <= t(n)
!c      3) t(k+1) < t(k+2) < ... < t(n-k)
!c      4) t(k+1) <= x(i) <= t(n-k)
!c      5) the conditions specified by schoenberg and whitney must hold
!c         for at least one subset of data points, i.e. there must be a
!c         subset of data points y(j) such that
!c             t(j) < y(j) < t(j+k+1), j=k+1,...,n-k-1
!c  ..
!c  ..scalar arguments..
      integer m,n,k,ier
!c  ..array arguments..
      real x(m),t(n)
!c  ..local scalars..
      integer i,i1,i2,j,j1,k1,k2,l,l1,l2,mm,m1,nk1,nk2
      real per,tj,tl,xi
!c  ..
      k1 = k+1
      k2 = k1+1
      nk1 = n-k1
      nk2 = nk1+1
      m1 = m-1
      ier = 10
!c  check condition no 1
      if(nk1.lt.k1 .or. n.gt.m+2*k) go to 130
!c  check condition no 2
      j = n
      do 20 i=1,k
        if(t(i).gt.t(i+1)) go to 130
        if(t(j).lt.t(j-1)) go to 130
        j = j-1
  20  continue
!c  check condition no 3
      do 30 i=k2,nk2
        if(t(i).le.t(i-1)) go to 130
  30  continue
!c  check condition no 4
      if(x(1).lt.t(k1) .or. x(m).gt.t(nk2)) go to 130
!c  check condition no 5
      l1 = k1
      l2 = 1
      do 50 l=1,m
         xi = x(l)
  40     if(xi.lt.t(l1+1) .or. l.eq.nk1) go to 50
         l1 = l1+1
         l2 = l2+1
         if(l2.gt.k1) go to 60
         go to 40
  50  continue
      l = m
  60  per = t(nk2)-t(k1)
      do 120 i1=2,l
         i = i1-1
         mm = i+m1
         do 110 j=k1,nk1
            tj = t(j)
            j1 = j+k1
            tl = t(j1)
  70        i = i+1
            if(i.gt.mm) go to 120
            i2 = i-m1
            if(i2) 80,80,90
  80        xi = x(i)
            go to 100
  90        xi = x(i2)+per
 100        if(xi.le.tj) go to 70
            if(xi.ge.tl) go to 120
 110     continue
         ier = 0
         go to 130
 120  continue
 130  return
      end
      subroutine fpcyt1(a,n,nn)
!c (l u)-decomposition of a cyclic tridiagonal matrix with the non-zero
!c elements stored as follows
!c
!c    | a(1,2) a(1,3)                                    a(1,1)  |
!c    | a(2,1) a(2,2) a(2,3)                                     |
!c    |        a(3,1) a(3,2) a(3,3)                              |
!c    |               ...............                            |
!c    |                               a(n-1,1) a(n-1,2) a(n-1,3) |
!c    | a(n,3)                                  a(n,1)   a(n,2)  |
!c
!c  ..
!c  ..scalar arguments..
      integer n,nn
!c  ..array arguments..
      real a(nn,6)
!c  ..local scalars..
      real aa,beta,gamma,sum,teta,v,one
      integer i,n1,n2
!c  ..
!c  set constant
      one = 1
      n2 = n-2
      beta = one/a(1,2)
      gamma = a(n,3)
      teta = a(1,1)*beta
      a(1,4) = beta
      a(1,5) = gamma
      a(1,6) = teta
      sum = gamma*teta
      do 10 i=2,n2
         v = a(i-1,3)*beta
         aa = a(i,1)
         beta = one/(a(i,2)-aa*v)
         gamma = -gamma*v
         teta = -teta*aa*beta
         a(i,4) = beta
         a(i,5) = gamma
         a(i,6) = teta
         sum = sum+gamma*teta
  10  continue
      n1 = n-1
      v = a(n2,3)*beta
      aa = a(n1,1)
      beta = one/(a(n1,2)-aa*v)
      gamma = a(n,1)-gamma*v
      teta = (a(n1,3)-teta*aa)*beta
      a(n1,4) = beta
      a(n1,5) = gamma
      a(n1,6) = teta
      a(n,4) = one/(a(n,2)-(sum+gamma*teta))
      return
      end
      subroutine fpcyt2(a,n,b,c,nn)
!c subroutine fpcyt2 solves a linear n x n system
!c         a * c = b
!c where matrix a is a cyclic tridiagonal matrix, decomposed
!c using subroutine fpsyt1.
!c  ..
!c  ..scalar arguments..
      integer n,nn
!c  ..array arguments..
      real a(nn,6),b(n),c(n)
!c  ..local scalars..
      real cc,sum
      integer i,j,j1,n1
!c  ..
      c(1) = b(1)*a(1,4)
      sum = c(1)*a(1,5)
      n1 = n-1
      do 10 i=2,n1
         c(i) = (b(i)-a(i,1)*c(i-1))*a(i,4)
         sum = sum+c(i)*a(i,5)
  10  continue
      cc = (b(n)-sum)*a(n,4)
      c(n) = cc
      c(n1) = c(n1)-cc*a(n1,6)
      j = n1
      do 20 i=3,n
         j1 = j-1
         c(j1) = c(j1)-c(j)*a(j1,3)*a(j1,4)-cc*a(j1,6)
         j = j1
  20  continue
      return
      end
      subroutine fpdisc(t,n,k2,b,nest)
!c  subroutine fpdisc calculates the discontinuity jumps of the kth
!c  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1)
!c  ..scalar arguments..
      integer n,k2,nest
!c  ..array arguments..
      real t(n),b(nest,k2)
!c  ..local scalars..
      real an,fac,prod
      integer i,ik,j,jk,k,k1,l,lj,lk,lmk,lp,nk1,nrint
!c  ..local array..
      real h(12)
!c  ..
      k1 = k2-1
      k = k1-1
      nk1 = n-k1
      nrint = nk1-k
      an = nrint
      fac = an/(t(nk1+1)-t(k1))
      do 40 l=k2,nk1
        lmk = l-k1
        do 10 j=1,k1
          ik = j+k1
          lj = l+j
          lk = lj-k2
          h(j) = t(l)-t(lk)
          h(ik) = t(l)-t(lj)
  10    continue
        lp = lmk
        do 30 j=1,k2
          jk = j
          prod = h(j)
          do 20 i=1,k
            jk = jk+1
            prod = prod*h(jk)*fac
  20      continue
          lk = lp+k1
          b(lmk,j) = (t(lk)-t(lp))/prod
          lp = lp+1
  30    continue
  40  continue
      return
      end
      subroutine fpgivs(piv,ww,cos,sin)
!c  subroutine fpgivs calculates the parameters of a givens
!c  transformation .
!c  ..
!c  ..scalar arguments..
      real piv,ww,cos,sin
!c  ..local scalars..
      real dd,one,store
!c  ..function references..
      real abs,sqrt
!c  ..
      one = 0.1e+01
      store = abs(piv)
      if(store.ge.ww) dd = store*sqrt(one+(ww/piv)**2)
      if(store.lt.ww) dd = ww*sqrt(one+(piv/ww)**2)
      cos = ww/dd
      sin = piv/dd
      ww = dd
      return
      end
      subroutine fpgrsp(ifsu,ifsv,ifbu,ifbv,iback,u,mu,v,mv,r,mr,dr,&
      iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm,mvnu,spu,spv,&
      right,q,au,av1,av2,bu,bv,a0,a1,b0,b1,c0,c1,cosi,nru,nrv)
!c  ..
!c  ..scalar arguments..
      real p,sq,fp
      integer ifsu,ifsv,ifbu,ifbv,iback,mu,mv,mr,iop0,iop1,nu,nv,nc,mm,mvnu
!c  ..array arguments..
      real u(mu),v(mv),r(mr),dr(6),tu(nu),tv(nv),c(nc),fpu(nu),fpv(nv),&
      spu(mu,4),spv(mv,4),right(mm),q(mvnu),au(nu,5),av1(nv,6),c0(nv),&
      av2(nv,4),a0(2,mv),b0(2,nv),cosi(2,nv),bu(nu,5),bv(nv,5),c1(nv),&
      a1(2,mv),b1(2,nv)
      integer nru(mu),nrv(mv)
!c  ..local scalars..
      real arg,co,dr01,dr02,dr03,dr11,dr12,dr13,fac,fac0,fac1,pinv,piv,si,term,one,three,half
      integer i,ic,ii,ij,ik,iq,irot,it,ir,i0,i1,i2,i3,j,jj,jk,jper,&
     j0,j1,k,k1,k2,l,l0,l1,l2,mvv,ncof,nrold,nroldu,nroldv,number,&
      numu,numu1,numv,numv1,nuu,nu4,nu7,nu8,nu9,nv11,nv4,nv7,nv8,n1
!c  ..local arrays..
      real h(5),h1(5),h2(4)
!c  ..function references..
      integer min0
      real cos,sin
!c  ..subroutine references..
!c    fpback,fpbspl,fpgivs,fpcyt1,fpcyt2,fpdisc,fpbacp,fprota
!c  ..
!c  let
!c               |     (spu)      |            |     (spv)      |
!c        (au) = | -------------- |     (av) = | -------------- |
!c               | sqrt(1/p) (bu) |            | sqrt(1/p) (bv) |
!c
!c                                | r  ' 0 |
!c                            q = | ------ |
!c                                | 0  ' 0 |
!c
!c  with c      : the (nu-4) x (nv-4) matrix which contains the b-spline
!c                coefficients.
!c       r      : the mu x mv matrix which contains the function values.
!c       spu,spv: the mu x (nu-4), resp. mv x (nv-4) observation matrices
!c                according to the least-squares problems in the u-,resp.
!c                v-direction.
!c       bu,bv  : the (nu-7) x (nu-4),resp. (nv-7) x (nv-4) matrices
!c                containing the discontinuity jumps of the derivatives
!c                of the b-splines in the u-,resp.v-variable at the knots
!c  the b-spline coefficients of the smoothing spline are then calculated
!c  as the least-squares solution of the following over-determined linear
!c  system of equations
!c
!c  (1)  (av) c (au)' = q
!c
!c  subject to the constraints
!c
!c  (2)  c(i,nv-3+j) = c(i,j), j=1,2,3 ; i=1,2,...,nu-4
!c
!c  (3)  if iop0 = 0  c(1,j) = dr(1)
!c          iop0 = 1  c(1,j) = dr(1)
!c                    c(2,j) = dr(1)+(dr(2)*cosi(1,j)+dr(3)*cosi(2,j))*
!c                            tu(5)/3. = c0(j) , j=1,2,...nv-4
!c
!c  (4)  if iop1 = 0  c(nu-4,j) = dr(4)
!c          iop1 = 1  c(nu-4,j) = dr(4)
!c                    c(nu-5,j) = dr(4)+(dr(5)*cosi(1,j)+dr(6)*cosi(2,j))
!c                                *(tu(nu-4)-tu(nu-3))/3. = c1(j)
!c
!c  set constants
      one = 1
      three = 3
      half = 0.5
!c  initialization
      nu4 = nu-4
      nu7 = nu-7
      nu8 = nu-8
      nu9 = nu-9
      nv4 = nv-4
      nv7 = nv-7
      nv8 = nv-8
      nv11 = nv-11
      nuu = nu4-iop0-iop1-2
      if(p.gt.0.) pinv = one/p
!c  it depends on the value of the flags ifsu,ifsv,ifbu,ifbv,iop0,iop1
!c  and on the value of p whether the matrices (spu), (spv), (bu), (bv),
!c  (cosi) still must be determined.
      if(ifsu.ne.0) go to 30
!c  calculate the non-zero elements of the matrix (spu) which is the ob-
!c  servation matrix according to the least-squares spline approximation
!c  problem in the u-direction.
      l = 4
      l1 = 5
      number = 0
      do 25 it=1,mu
        arg = u(it)
  10    if(arg.lt.tu(l1) .or. l.eq.nu4) go to 15
        l = l1
        l1 = l+1
        number = number+1
        go to 10
  15    call fpbspl(tu,nu,3,arg,l,h)
        do 20 i=1,4
          spu(it,i) = h(i)
  20    continue
        nru(it) = number
  25  continue
      ifsu = 1
!c  calculate the non-zero elements of the matrix (spv) which is the ob-
!c  servation matrix according to the least-squares spline approximation
!c  problem in the v-direction.
  30  if(ifsv.ne.0) go to 85
      l = 4
      l1 = 5
      number = 0
      do 50 it=1,mv
        arg = v(it)
  35    if(arg.lt.tv(l1) .or. l.eq.nv4) go to 40
        l = l1
        l1 = l+1
        number = number+1
        go to 35
  40    call fpbspl(tv,nv,3,arg,l,h)
        do 45 i=1,4
          spv(it,i) = h(i)
  45    continue
        nrv(it) = number
  50  continue
      ifsv = 1
      if(iop0.eq.0 .and. iop1.eq.0) go to 85
!c  calculate the coefficients of the interpolating splines for cos(v)
!c  and sin(v).
      do 55 i=1,nv4
         cosi(1,i) = 0.
         cosi(2,i) = 0.
  55  continue
      if(nv7.lt.4) go to 85
      do 65 i=1,nv7
         l = i+3
         arg = tv(l)
         call fpbspl(tv,nv,3,arg,l,h)
         do 60 j=1,3
            av1(i,j) = h(j)
  60     continue
         cosi(1,i) = cos(arg)
         cosi(2,i) = sin(arg)
  65  continue
      call fpcyt1(av1,nv7,nv)
      do 80 j=1,2
         do 70 i=1,nv7
            right(i) = cosi(j,i)
  70     continue
         call fpcyt2(av1,nv7,right,right,nv)
         do 75 i=1,nv7
            cosi(j,i+1) = right(i)
  75     continue
         cosi(j,1) = cosi(j,nv7+1)
         cosi(j,nv7+2) = cosi(j,2)
         cosi(j,nv4) = cosi(j,3)
  80  continue
  85  if(p.le.0.) go to  150
!c  calculate the non-zero elements of the matrix (bu).
      if(ifbu.ne.0 .or. nu8.eq.0) go to 90
      call fpdisc(tu,nu,5,bu,nu)
      ifbu = 1
!c  calculate the non-zero elements of the matrix (bv).
  90  if(ifbv.ne.0 .or. nv8.eq.0) go to 150
      call fpdisc(tv,nv,5,bv,nv)
      ifbv = 1
!c  substituting (2),(3) and (4) into (1), we obtain the overdetermined
!c  system
!c         (5)  (avv) (cc) (auu)' = (qq)
!c  from which the nuu*nv7 remaining coefficients
!c         c(i,j) , i=2+iop0,3+iop0,...,nu-5-iop1,j=1,2,...,nv-7.
!c  the elements of (cc), are then determined in the least-squares sense.
!c  simultaneously, we compute the resulting sum of squared residuals sq.
 150  dr01 = dr(1)
      dr11 = dr(4)
      do 155 i=1,mv
         a0(1,i) = dr01
         a1(1,i) = dr11
 155  continue
      if(nv8.eq.0 .or. p.le.0.) go to 165
      do 160 i=1,nv8
         b0(1,i) = 0.
         b1(1,i) = 0.
 160  continue
 165  mvv = mv
      if(iop0.eq.0) go to 195
      fac = (tu(5)-tu(4))/three
      dr02 = dr(2)*fac
      dr03 = dr(3)*fac
      do 170 i=1,nv4
         c0(i) = dr01+dr02*cosi(1,i)+dr03*cosi(2,i)
 170  continue
      do 180 i=1,mv
         number = nrv(i)
         fac = 0.
         do 175 j=1,4
            number = number+1
            fac = fac+c0(number)*spv(i,j)
 175     continue
         a0(2,i) = fac
 180  continue
      if(nv8.eq.0 .or. p.le.0.) go to 195
      do 190 i=1,nv8
         number = i
         fac = 0.
         do 185 j=1,5
            fac = fac+c0(number)*bv(i,j)
            number = number+1
 185     continue
         b0(2,i) = fac*pinv
 190  continue
      mvv = mv+nv8
 195  if(iop1.eq.0) go to 225
      fac = (tu(nu4)-tu(nu4+1))/three
      dr12 = dr(5)*fac
      dr13 = dr(6)*fac
      do 200 i=1,nv4
         c1(i) = dr11+dr12*cosi(1,i)+dr13*cosi(2,i)
 200  continue
      do 210 i=1,mv
         number = nrv(i)
         fac = 0.
         do 205 j=1,4
            number = number+1
            fac = fac+c1(number)*spv(i,j)
 205     continue
         a1(2,i) = fac
 210  continue
      if(nv8.eq.0 .or. p.le.0.) go to 225
      do 220 i=1,nv8
         number = i
         fac = 0.
         do 215 j=1,5
            fac = fac+c1(number)*bv(i,j)
            number = number+1
 215     continue
         b1(2,i) = fac*pinv
 220  continue
      mvv = mv+nv8
!c  we first determine the matrices (auu) and (qq). then we reduce the
!c  matrix (auu) to an unit upper triangular form (ru) using givens
!c  rotations without square roots. we apply the same transformations to
!c  the rows of matrix qq to obtain the mv x nuu matrix g.
!c  we store matrix (ru) into au and g into q.
 225  l = mvv*nuu
!c  initialization.
      sq = 0.
      if(l.eq.0) go to 245
      do 230 i=1,l
        q(i) = 0.
 230  continue
      do 240 i=1,nuu
        do 240 j=1,5
          au(i,j) = 0.
 240  continue
      l = 0
 245  nrold = 0
      n1 = nrold+1
      do 420 it=1,mu
        number = nru(it)
!c  find the appropriate column of q.
 250    do 260 j=1,mvv
           right(j) = 0.
 260    continue
        if(nrold.eq.number) go to 280
        if(p.le.0.) go to 410
!c  fetch a new row of matrix (bu).
        do 270 j=1,5
          h(j) = bu(n1,j)*pinv
 270    continue
        i0 = 1
        i1 = 5
        go to 310
!c  fetch a new row of matrix (spu).
 280    do 290 j=1,4
          h(j) = spu(it,j)
 290    continue
!c  find the appropriate column of q.
        do 300 j=1,mv
          l = l+1
          right(j) = r(l)
 300    continue
        i0 = 1
        i1 = 4
 310    j0 = n1
        j1 = nu7-number
!c  take into account that we eliminate the constraints (3)
 315     if(j0-1.gt.iop0) go to 335
         fac0 = h(i0)
         do 320 j=1,mv
            right(j) = right(j)-fac0*a0(j0,j)
 320     continue
         if(mv.eq.mvv) go to 330
         j = mv
         do 325 jj=1,nv8
            j = j+1
            right(j) = right(j)-fac0*b0(j0,jj)
 325     continue
 330     j0 = j0+1
         i0 = i0+1
         go to 315
!c  take into account that we eliminate the constraints (4)
 335     if(j1-1.gt.iop1) go to 360
         fac1 = h(i1)
         do 340 j=1,mv
            right(j) = right(j)-fac1*a1(j1,j)
 340     continue
         if(mv.eq.mvv) go to 350
         j = mv
         do 345 jj=1,nv8
            j = j+1
            right(j) = right(j)-fac1*b1(j1,jj)
 345     continue
 350     j1 = j1+1
         i1 = i1-1
         go to 335
 360     irot = nrold-iop0-1
         if(irot.lt.0) irot = 0
!c  rotate the new row of matrix (auu) into triangle.
        if(i0.gt.i1) go to 390
        do 385 i=i0,i1
          irot = irot+1
          piv = h(i)
          if(piv.eq.0.) go to 385
!c  calculate the parameters of the givens transformation.
          call fpgivs(piv,au(irot,1),co,si)
!c  apply that transformation to the rows of matrix (qq).
          iq = (irot-1)*mvv
          do 370 j=1,mvv
            iq = iq+1
            call fprota(co,si,right(j),q(iq))
 370      continue
!c  apply that transformation to the columns of (auu).
          if(i.eq.i1) go to 385
          i2 = 1
          i3 = i+1
          do 380 j=i3,i1
            i2 = i2+1
            call fprota(co,si,h(j),au(irot,i2))
 380      continue
 385    continue
!c  we update the sum of squared residuals.
 390    do 395 j=1,mvv
          sq = sq+right(j)**2
 395    continue
 400    if(nrold.eq.number) go to 420
 410    nrold = n1
        n1 = n1+1
        go to 250
 420  continue
      if(nuu.eq.0) go to 800
!c  we determine the matrix (avv) and then we reduce her to an unit
!c  upper triangular form (rv) using givens rotations without square
!c  roots. we apply the same transformations to the columns of matrix
!c  g to obtain the (nv-7) x (nu-6-iop0-iop1) matrix h.
!c  we store matrix (rv) into av1 and av2, h into c.
!c  the nv7 x nv7 triangular unit upper matrix (rv) has the form
!c              | av1 '     |
!c       (rv) = |     ' av2 |
!c              |  0  '     |
!c  with (av2) a nv7 x 4 matrix and (av1) a nv11 x nv11 unit upper
!c  triangular matrix of bandwidth 5.
      ncof = nuu*nv7
!c  initialization.
      do 430 i=1,ncof
        c(i) = 0.
 430  continue
      do 440 i=1,nv4
        av1(i,5) = 0.
        do 440 j=1,4
          av1(i,j) = 0.
          av2(i,j) = 0.
 440  continue
      jper = 0
      nrold = 0
      do 770 it=1,mv
        number = nrv(it)
 450    if(nrold.eq.number) go to 480
        if(p.le.0.) go to 760
!c  fetch a new row of matrix (bv).
        n1 = nrold+1
        do 460 j=1,5
          h(j) = bv(n1,j)*pinv
 460    continue
!c  find the appropiate row of g.
        do 465 j=1,nuu
          right(j) = 0.
 465    continue
        if(mv.eq.mvv) go to 510
        l = mv+n1
        do 470 j=1,nuu
          right(j) = q(l)
          l = l+mvv
 470    continue
        go to 510
!c  fetch a new row of matrix (spv)
 480    h(5) = 0.
        do 490 j=1,4
          h(j) = spv(it,j)
 490    continue
!c  find the appropiate row of g.
        l = it
        do 500 j=1,nuu
          right(j) = q(l)
          l = l+mvv
 500    continue
!c  test whether there are non-zero values in the new row of (avv)
!c  corresponding to the b-splines n(j;v),j=nv7+1,...,nv4.
 510     if(nrold.lt.nv11) go to 710
         if(jper.ne.0) go to 550
!c  initialize the matrix (av2).
         jk = nv11+1
         do 540 i=1,4
            ik = jk
            do 520 j=1,5
               if(ik.le.0) go to 530
               av2(ik,i) = av1(ik,j)
               ik = ik-1
 520        continue
 530        jk = jk+1
 540     continue
         jper = 1
!c  if one of the non-zero elements of the new row corresponds to one of
!c  the b-splines n(j;v),j=nv7+1,...,nv4, we take account of condition
!c  (2) for setting up this row of (avv). the row is stored in h1( the
!c  part with respect to av1) and h2 (the part with respect to av2).
 550     do 560 i=1,4
            h1(i) = 0.
            h2(i) = 0.
 560     continue
         h1(5) = 0.
         j = nrold-nv11
         do 600 i=1,5
            j = j+1
            l0 = j
 570        l1 = l0-4
            if(l1.le.0) go to 590
            if(l1.le.nv11) go to 580
            l0 = l1-nv11
            go to 570
 580        h1(l1) = h(i)
            go to 600
 590        h2(l0) = h2(l0) + h(i)
 600     continue
!c  rotate the new row of (avv) into triangle.
         if(nv11.le.0) go to 670
!c  rotations with the rows 1,2,...,nv11 of (avv).
         do 660 j=1,nv11
            piv = h1(1)
            i2 = min0(nv11-j,4)
            if(piv.eq.0.) go to 640
!c  calculate the parameters of the givens transformation.
            call fpgivs(piv,av1(j,1),co,si)
!c  apply that transformation to the columns of matrix g.
            ic = j
            do 610 i=1,nuu
               call fprota(co,si,right(i),c(ic))
               ic = ic+nv7
 610        continue
!c  apply that transformation to the rows of (avv) with respect to av2.
            do 620 i=1,4
               call fprota(co,si,h2(i),av2(j,i))
 620        continue
!c  apply that transformation to the rows of (avv) with respect to av1.
            if(i2.eq.0) go to 670
            do 630 i=1,i2
               i1 = i+1
               call fprota(co,si,h1(i1),av1(j,i1))
 630        continue
 640        do 650 i=1,i2
               h1(i) = h1(i+1)
 650        continue
            h1(i2+1) = 0.
 660     continue
!c  rotations with the rows nv11+1,...,nv7 of avv.
 670     do 700 j=1,4
            ij = nv11+j
            if(ij.le.0) go to 700
            piv = h2(j)
            if(piv.eq.0.) go to 700
!c  calculate the parameters of the givens transformation.
            call fpgivs(piv,av2(ij,j),co,si)
!c  apply that transformation to the columns of matrix g.
            ic = ij
            do 680 i=1,nuu
               call fprota(co,si,right(i),c(ic))
               ic = ic+nv7
 680        continue
            if(j.eq.4) go to 700
!c  apply that transformation to the rows of (avv) with respect to av2.
            j1 = j+1
            do 690 i=j1,4
               call fprota(co,si,h2(i),av2(ij,i))
 690        continue
 700     continue
!c  we update the sum of squared residuals.
         do 705 i=1,nuu
           sq = sq+right(i)**2
 705     continue
         go to 750
!c  rotation into triangle of the new row of (avv), in case the elements
!c  corresponding to the b-splines n(j;v),j=nv7+1,...,nv4 are all zero.
 710     irot =nrold
         do 740 i=1,5
            irot = irot+1
            piv = h(i)
            if(piv.eq.0.) go to 740
!c  calculate the parameters of the givens transformation.
            call fpgivs(piv,av1(irot,1),co,si)
!c  apply that transformation to the columns of matrix g.
            ic = irot
            do 720 j=1,nuu
               call fprota(co,si,right(j),c(ic))
               ic = ic+nv7
 720        continue
!c  apply that transformation to the rows of (avv).
            if(i.eq.5) go to 740
            i2 = 1
            i3 = i+1
            do 730 j=i3,5
               i2 = i2+1
               call fprota(co,si,h(j),av1(irot,i2))
 730        continue
 740     continue
!c  we update the sum of squared residuals.
         do 745 i=1,nuu
           sq = sq+right(i)**2
 745     continue
 750     if(nrold.eq.number) go to 770
 760     nrold = nrold+1
         go to 450
 770  continue
!c  test whether the b-spline coefficients must be determined.
      if(iback.ne.0) return
!c  backward substitution to obtain the b-spline coefficients as the
!c  solution of the linear system    (rv) (cr) (ru)' = h.
!c  first step: solve the system  (rv) (c1) = h.
      k = 1
      do 780 i=1,nuu
         call fpbacp(av1,av2,c(k),nv7,4,c(k),5,nv)
         k = k+nv7
 780  continue
!c  second step: solve the system  (cr) (ru)' = (c1).
      k = 0
      do 795 j=1,nv7
        k = k+1
        l = k
        do 785 i=1,nuu
          right(i) = c(l)
          l = l+nv7
 785    continue
        call fpback(au,right,nuu,5,right,nu)
        l = k
        do 790 i=1,nuu
           c(l) = right(i)
           l = l+nv7
 790    continue
 795  continue
!c  calculate from the conditions (2)-(3)-(4), the remaining b-spline
!c  coefficients.
 800  ncof = nu4*nv4
      j = ncof
      do 805 l=1,nv4
         q(l) = dr01
         q(j) = dr11
         j = j-1
 805  continue
      i = nv4
      j = 0
      if(iop0.eq.0) go to 815
      do 810 l=1,nv4
         i = i+1
         q(i) = c0(l)
 810  continue
 815  if(nuu.eq.0) go to 835
      do 830 l=1,nuu
         ii = i
         do 820 k=1,nv7
            i = i+1
            j = j+1
            q(i) = c(j)
 820     continue
         do 825 k=1,3
            ii = ii+1
            i = i+1
            q(i) = q(ii)
 825     continue
 830  continue
 835  if(iop1.eq.0) go to 845
      do 840 l=1,nv4
         i = i+1
         q(i) = c1(l)
 840  continue
 845  do 850 i=1,ncof
         c(i) = q(i)
 850  continue
!c  calculate the quantities
!c    res(i,j) = (r(i,j) - s(u(i),v(j)))**2 , i=1,2,..,mu;j=1,2,..,mv
!c    fp = sumi=1,mu(sumj=1,mv(res(i,j)))
!c    fpu(r) = sum''i(sumj=1,mv(res(i,j))) , r=1,2,...,nu-7
!c                  tu(r+3) <= u(i) <= tu(r+4)
!c    fpv(r) = sumi=1,mu(sum''j(res(i,j))) , r=1,2,...,nv-7
!c                  tv(r+3) <= v(j) <= tv(r+4)
      fp = 0.
      do 890 i=1,nu
        fpu(i) = 0.
 890  continue
      do 900 i=1,nv
        fpv(i) = 0.
 900  continue
      ir = 0
      nroldu = 0
!c  main loop for the different grid points.
      do 950 i1=1,mu
        numu = nru(i1)
        numu1 = numu+1
        nroldv = 0
        do 940 i2=1,mv
          numv = nrv(i2)
          numv1 = numv+1
          ir = ir+1
!c  evaluate s(u,v) at the current grid point by making the sum of the
!c  cross products of the non-zero b-splines at (u,v), multiplied with
!c  the appropiate b-spline coefficients.
          term = 0.
          k1 = numu*nv4+numv
          do 920 l1=1,4
            k2 = k1
            fac = spu(i1,l1)
            do 910 l2=1,4
              k2 = k2+1
              term = term+fac*spv(i2,l2)*c(k2)
 910        continue
            k1 = k1+nv4
 920      continue
!c  calculate the squared residual at the current grid point.
          term = (r(ir)-term)**2
!c  adjust the different parameters.
          fp = fp+term
          fpu(numu1) = fpu(numu1)+term
          fpv(numv1) = fpv(numv1)+term
          fac = term*half
          if(numv.eq.nroldv) go to 930
          fpv(numv1) = fpv(numv1)-fac
          fpv(numv) = fpv(numv)+fac
 930      nroldv = numv
          if(numu.eq.nroldu) go to 940
          fpu(numu1) = fpu(numu1)-fac
          fpu(numu) = fpu(numu)+fac
 940    continue
        nroldu = numu
 950  continue
      return
      end
      subroutine fpknot(x,m,t,n,fpint,nrdata,nrint,nest,istart)
!c  subroutine fpknot locates an additional knot for a spline of degree
!c  k and adjusts the corresponding parameters,i.e.
!c    t     : the position of the knots.
!c    n     : the number of knots.
!c    nrint : the number of knotintervals.
!c    fpint : the sum of squares of residual right hand sides
!c            for each knot interval.
!c    nrdata: the number of data points inside each knot interval.
!c  istart indicates that the smallest data point at which the new knot
!c  may be added is x(istart+1)
!c  ..
!c  ..scalar arguments..
      integer m,n,nrint,nest,istart
!c  ..array arguments..
      real x(m),t(nest),fpint(nest)
      integer nrdata(nest)
!c  ..local scalars..
      real an,am,fpmax
      integer ihalf,j,jbegin,jj,jk,jpoint,k,maxbeg,maxpt, next,nrx,number
!c  ..
      k = (n-nrint-1)/2
!c  search for knot interval t(number+k) <= x <= t(number+k+1) where
!c  fpint(number) is maximal on the condition that nrdata(number)
!c  not equals zero.
      fpmax = 0.
      jbegin = istart
      do 20 j=1,nrint
        jpoint = nrdata(j)
        if(fpmax.ge.fpint(j) .or. jpoint.eq.0) go to 10
        fpmax = fpint(j)
        number = j
        maxpt = jpoint
        maxbeg = jbegin
  10    jbegin = jbegin+jpoint+1
  20  continue
!c  let coincide the new knot t(number+k+1) with a data point x(nrx)
!c  inside the old knot interval t(number+k) <= x <= t(number+k+1).
      ihalf = maxpt/2+1
      nrx = maxbeg+ihalf
      next = number+1
      if(next.gt.nrint) go to 40
!c  adjust the different parameters.
      do 30 j=next,nrint
        jj = next+nrint-j
        fpint(jj+1) = fpint(jj)
        nrdata(jj+1) = nrdata(jj)
        jk = jj+k
        t(jk+1) = t(jk)
  30  continue
  40  nrdata(number) = ihalf-1
      nrdata(next) = maxpt-ihalf
      am = maxpt
      an = nrdata(number)
      fpint(number) = fpmax*an/am
      an = nrdata(next)
      fpint(next) = fpmax*an/am
      jk = next+k
      t(jk) = x(nrx)
      n = n+1
      nrint = nrint+1
      return
      end
      subroutine fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,mr,r0,r1,dr,&
      iopt,ider,tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpu,fpv,&
      nru,nrv,wrk,lwrk)
!c  given the set of function values r(i,j) defined on the rectangular
!c  grid (u(i),v(j)),i=1,2,...,mu;j=1,2,...,mv, fpopsp determines a
!c  smooth bicubic spline approximation with given knots tu(i),i=1,..,nu
!c  in the u-direction and tv(j),j=1,2,...,nv in the v-direction. this
!c  spline sp(u,v) will be periodic in the variable v and will satisfy
!c  the following constraints
!c
!c     s(tu(1),v) = dr(1) , tv(4) <=v<= tv(nv-3)
!c
!c     s(tu(nu),v) = dr(4) , tv(4) <=v<= tv(nv-3)
!c
!c  and (if iopt(2) = 1)
!c
!c     d s(tu(1),v)
!c     ------------ =  dr(2)*cos(v)+dr(3)*sin(v) , tv(4) <=v<= tv(nv-3)
!c     d u
!c
!c  and (if iopt(3) = 1)
!c
!c     d s(tu(nu),v)
!c     ------------- =  dr(5)*cos(v)+dr(6)*sin(v) , tv(4) <=v<= tv(nv-3)
!c     d u
!c
!c  where the parameters dr(i) correspond to the derivative values at the
!c  poles as defined in subroutine spgrid.
!c
!c  the b-spline coefficients of sp(u,v) are determined as the least-
!c  squares solution  of an overdetermined linear system which depends
!c  on the value of p and on the values dr(i),i=1,...,6. the correspond-
!c  ing sum of squared residuals sq is a simple quadratic function in
!c  the variables dr(i). these may or may not be provided. the values
!c  dr(i) which are not given will be determined so as to minimize the
!c  resulting sum of squared residuals sq. in that case the user must
!c  provide some initial guess dr(i) and some estimate (dr(i)-step,
!c  dr(i)+step) of the range of possible values for these latter.
!c
!c  sp(u,v) also depends on the parameter p (p>0) in such a way that
!c    - if p tends to infinity, sp(u,v) becomes the least-squares spline
!c      with given knots, satisfying the constraints.
!c    - if p tends to zero, sp(u,v) becomes the least-squares polynomial,
!c      satisfying the constraints.
!c    - the function  f(p)=sumi=1,mu(sumj=1,mv((r(i,j)-sp(u(i),v(j)))**2)
!c      is continuous and strictly decreasing for p>0.
!c
!c  ..scalar arguments..
      integer ifsu,ifsv,ifbu,ifbv,mu,mv,mr,nu,nv,nuest,nvest, nc,lwrk
      real r0,r1,p,fp
!c  ..array arguments..
      integer ider(4),nru(mu),nrv(mv),iopt(3)
      real u(mu),v(mv),r(mr),dr(6),tu(nu),tv(nv),c(nc),fpu(nu),fpv(nv), wrk(lwrk),step(2)
!c  ..local scalars..
      real res,sq,sqq,sq0,sq1,step1,step2,three
      integer i,id0,iop0,iop1,i1,j,l,lau,lav1,lav2,la0,la1,lbu,lbv,lb0,&
        lb1,lc0,lc1,lcs,lq,lri,lsu,lsv,l1,l2,mm,mvnu,number
!c  ..local arrays..
      integer nr(6)
      real delta(6),drr(6),sum(6),a(6,6),g(6)
!c  ..function references..
      integer max0
!c  ..subroutine references..
!c    fpgrsp,fpsysy
!c  ..
!c  set constant
      three = 3
!c  we partition the working space
      lsu = 1
      lsv = lsu+4*mu
      lri = lsv+4*mv
      mm = max0(nuest,mv+nvest)
      lq = lri+mm
      mvnu = nuest*(mv+nvest-8)
      lau = lq+mvnu
      lav1 = lau+5*nuest
      lav2 = lav1+6*nvest
      lbu = lav2+4*nvest
      lbv = lbu+5*nuest
      la0 = lbv+5*nvest
      la1 = la0+2*mv
      lb0 = la1+2*mv
      lb1 = lb0+2*nvest
      lc0 = lb1+2*nvest
      lc1 = lc0+nvest
      lcs = lc1+nvest
!c  we calculate the smoothing spline sp(u,v) according to the input
!c  values dr(i),i=1,...,6.
      iop0 = iopt(2)
      iop1 = iopt(3)
      id0 = ider(1)
      id1 = ider(3)
      call fpgrsp(ifsu,ifsv,ifbu,ifbv,0,u,mu,v,mv,r,mr,dr,&
      iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm,mvnu,&
      wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),&
      wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0),&
      wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
      sq0 = 0.
      sq1 = 0.
      if(id0.eq.0) sq0 = (r0-dr(1))**2
      if(id1.eq.0) sq1 = (r1-dr(4))**2
      sq = sq+sq0+sq1
!c in case all derivative values dr(i) are given (step<=0) or in case
!c we have spline interpolation, we accept this spline as a solution.
      if(sq.le.0.) return
      if(step(1).le.0. .and. step(2).le.0.) return
      do 10 i=1,6
        drr(i) = dr(i)
  10  continue
!c number denotes the number of derivative values dr(i) that still must
!c be optimized. let us denote these parameters by g(j),j=1,...,number.
      number = 0
      if(id0.gt.0) go to 20
      number = 1
      nr(1) = 1
      delta(1) = step(1)
  20  if(iop0.eq.0) go to 30
      if(ider(2).ne.0) go to 30
      step2 = step(1)*three/(tu(5)-tu(4))
      nr(number+1) = 2
      nr(number+2) = 3
      delta(number+1) = step2
      delta(number+2) = step2
      number = number+2
  30  if(id1.gt.0) go to 40
      number = number+1
      nr(number) = 4
      delta(number) = step(2)
  40  if(iop1.eq.0) go to 50
      if(ider(4).ne.0) go to 50
      step2 = step(2)*three/(tu(nu)-tu(nu-4))
      nr(number+1) = 5
      nr(number+2) = 6
      delta(number+1) = step2
      delta(number+2) = step2
      number = number+2
  50  if(number.eq.0) return
!c the sum of squared residulas sq is a quadratic polynomial in the
!c parameters g(j). we determine the unknown coefficients of this
!c polymomial by calculating (number+1)*(number+2)/2 different splines
!c according to specific values for g(j).
      do 60 i=1,number
         l = nr(i)
         step1 = delta(i)
         drr(l) = dr(l)+step1
         call fpgrsp(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,r,mr,drr,&
         iop0,iop1,tu,nu,tv,nv,p,c,nc,sum(i),fp,fpu,fpv,mm,mvnu,&
         wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),&
         wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0),&
         wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
         if(id0.eq.0) sq0 = (r0-drr(1))**2
         if(id1.eq.0) sq1 = (r1-drr(4))**2
         sum(i) = sum(i)+sq0+sq1
         drr(l) = dr(l)-step1
         call fpgrsp(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,r,mr,drr,&
         iop0,iop1,tu,nu,tv,nv,p,c,nc,sqq,fp,fpu,fpv,mm,mvnu,&
         wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),&
         wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0),&
         wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
         if(id0.eq.0) sq0 = (r0-drr(1))**2
         if(id1.eq.0) sq1 = (r1-drr(4))**2
         sqq = sqq+sq0+sq1
         drr(l) = dr(l)
         a(i,i) = (sum(i)+sqq-sq-sq)/step1**2
         if(a(i,i).le.0.) go to 110
         g(i) = (sqq-sum(i))/(step1+step1)
  60  continue
      if(number.eq.1) go to 90
      do 80 i=2,number
         l1 = nr(i)
         step1 = delta(i)
         drr(l1) = dr(l1)+step1
         i1 = i-1
         do 70 j=1,i1
            l2 = nr(j)
            step2 = delta(j)
            drr(l2) = dr(l2)+step2
            call fpgrsp(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,r,mr,drr,&
           iop0,iop1,tu,nu,tv,nv,p,c,nc,sqq,fp,fpu,fpv,mm,mvnu,&
            wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),&
            wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0),&
            wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
            if(id0.eq.0) sq0 = (r0-drr(1))**2
            if(id1.eq.0) sq1 = (r1-drr(4))**2
            sqq = sqq+sq0+sq1
            a(i,j) = (sq+sqq-sum(i)-sum(j))/(step1*step2)
            drr(l2) = dr(l2)
  70     continue
         drr(l1) = dr(l1)
  80  continue
!c the optimal values g(j) are found as the solution of the system
!c d (sq) / d (g(j)) = 0 , j=1,...,number.
  90  call fpsysy(a,number,g)
      do 100 i=1,number
         l = nr(i)
         dr(l) = dr(l)+g(i)
 100  continue
!c we determine the spline sp(u,v) according to the optimal values g(j).
 110  call fpgrsp(ifsu,ifsv,ifbu,ifbv,0,u,mu,v,mv,r,mr,dr,&
      iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm,mvnu,&
      wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),&
      wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0),&
      wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
      if(id0.eq.0) sq0 = (r0-dr(1))**2
      if(id1.eq.0) sq1 = (r1-dr(4))**2
      sq = sq+sq0+sq1
      return
      end
      subroutine fporde(x,y,m,kx,ky,tx,nx,ty,ny,nummer,index,nreg)
!c  subroutine fporde sorts the data points (x(i),y(i)),i=1,2,...,m
!c  according to the panel tx(l)<=x<tx(l+1),ty(k)<=y<ty(k+1), they belong
!c  to. for each panel a stack is constructed  containing the numbers
!c  of data points lying inside; index(j),j=1,2,...,nreg points to the
!c  first data point in the jth panel while nummer(i),i=1,2,...,m gives
!c  the number of the next data point in the panel.
!c  ..
!c  ..scalar arguments..
      integer m,kx,ky,nx,ny,nreg
!c  ..array arguments..
      real x(m),y(m),tx(nx),ty(ny)
      integer nummer(m),index(nreg)
!c  ..local scalars..
      real xi,yi
      integer i,im,k,kx1,ky1,k1,l,l1,nk1x,nk1y,num,nyy
!c  ..
      kx1 = kx+1
      ky1 = ky+1
      nk1x = nx-kx1
      nk1y = ny-ky1
      nyy = nk1y-ky
      do 10 i=1,nreg
        index(i) = 0
  10  continue
      do 60 im=1,m
        xi = x(im)
        yi = y(im)
        l = kx1
        l1 = l+1
  20    if(xi.lt.tx(l1) .or. l.eq.nk1x) go to 30
        l = l1
        l1 = l+1
        go to 20
  30    k = ky1
        k1 = k+1
  40    if(yi.lt.ty(k1) .or. k.eq.nk1y) go to 50
        k = k1
        k1 = k+1
        go to 40
  50    num = (l-kx1)*nyy+k-ky
        nummer(im) = index(num)
        index(num) = im
  60  continue
      return
      end

      subroutine fprank(a,f,n,m,na,tol,c,sq,rank,aa,ff,h)
!c  subroutine fprank finds the minimum norm solution of a least-
!c  squares problem in case of rank deficiency.
!c
!c  input parameters:
!c    a : array, which contains the non-zero elements of the observation
!c        matrix after triangularization by givens transformations.
!c    f : array, which contains the transformed right hand side.
!c    n : integer,wich contains the dimension of a.
!c    m : integer, which denotes the bandwidth of a.
!c  tol : real value, giving a threshold to determine the rank of a.
!c
!c  output parameters:
!c    c : array, which contains the minimum norm solution.
!c   sq : real value, giving the contribution of reducing the rank
!c        to the sum of squared residuals.
!c rank : integer, which contains the rank of matrix a.
!c
!c  ..scalar arguments..
      integer n,m,na,rank
      real tol,sq
!c  ..array arguments..
      real a(na,m),f(n),c(n),aa(n,m),ff(n),h(m)
!c  ..local scalars..
      integer i,ii,ij,i1,i2,j,jj,j1,j2,j3,k,kk,m1,nl
      real cos,fac,piv,sin,yi
      double precision store,stor1,stor2,stor3
!c  ..function references..
      integer min0
!c  ..subroutine references..
!c    fpgivs,fprota
!c  ..
      m1 = m-1
!c  the rank deficiency nl is considered to be the number of sufficient
!c  small diagonal elements of a.
      nl = 0
      sq = 0.
      do 90 i=1,n
        if(a(i,1).gt.tol) go to 90
!c  if a sufficient small diagonal element is found, we put it to
!c  zero. the remainder of the row corresponding to that zero diagonal
!c  element is then rotated into triangle by givens rotations .
!c  the rank deficiency is increased by one.
        nl = nl+1
        if(i.eq.n) go to 90
        yi = f(i)
        do 10 j=1,m1
          h(j) = a(i,j+1)
  10    continue
        h(m) = 0.
        i1 = i+1
        do 60 ii=i1,n
          i2 = min0(n-ii,m1)
          piv = h(1)
          if(piv.eq.0.) go to 30
          call fpgivs(piv,a(ii,1),cos,sin)
          call fprota(cos,sin,yi,f(ii))
          if(i2.eq.0) go to 70
          do 20 j=1,i2
            j1 = j+1
            call fprota(cos,sin,h(j1),a(ii,j1))
            h(j) = h(j1)
  20      continue
          go to 50
  30      if(i2.eq.0) go to 70
          do 40 j=1,i2
            h(j) = h(j+1)
  40      continue
  50      h(i2+1) = 0.
  60    continue
!c  add to the sum of squared residuals the contribution of deleting
!c  the row with small diagonal element.
  70    sq = sq+yi**2
  90  continue
!c  rank denotes the rank of a.
      rank = n-nl
!c  let b denote the (rank*n) upper trapezoidal matrix which can be
!c  obtained from the (n*n) upper triangular matrix a by deleting
!c  the rows and interchanging the columns corresponding to a zero
!c  diagonal element. if this matrix is factorized using givens
!c  transformations as  b = (r) (u)  where
!c    r is a (rank*rank) upper triangular matrix,
!c    u is a (rank*n) orthonormal matrix
!c  then the minimal least-squares solution c is given by c = b' v,
!c  where v is the solution of the system  (r) (r)' v = g  and
!c  g denotes the vector obtained from the old right hand side f, by
!c  removing the elements corresponding to a zero diagonal element of a.
!c  initialization.
      do 100 i=1,rank
        do 100 j=1,m
          aa(i,j) = 0.
 100  continue
!c  form in aa the upper triangular matrix obtained from a by
!c  removing rows and columns with zero diagonal elements. form in ff
!c  the new right hand side by removing the elements of the old right
!c  hand side corresponding to a deleted row.
      ii = 0
      do 120 i=1,n
        if(a(i,1).le.tol) go to 120
        ii = ii+1
        ff(ii) = f(i)
        aa(ii,1) = a(i,1)
        jj = ii
        kk = 1
        j = i
        j1 = min0(j-1,m1)
        if(j1.eq.0) go to 120
        do 110 k=1,j1
          j = j-1
          if(a(j,1).le.tol) go to 110
          kk = kk+1
          jj = jj-1
          aa(jj,kk) = a(j,k+1)
 110    continue
 120  continue
!c  form successively in h the columns of a with a zero diagonal element.
      ii = 0
      do 200 i=1,n
        ii = ii+1
        if(a(i,1).gt.tol) go to 200
        ii = ii-1
        if(ii.eq.0) go to 200
        jj = 1
        j = i
        j1 = min0(j-1,m1)
        do 130 k=1,j1
          j = j-1
          if(a(j,1).le.tol) go to 130
          h(jj) = a(j,k+1)
          jj = jj+1
 130    continue
        do 140 kk=jj,m
          h(kk) = 0.
 140    continue
!c  rotate this column into aa by givens transformations.
        jj = ii
        do 190 i1=1,ii
          j1 = min0(jj-1,m1)
          piv = h(1)
          if(piv.ne.0.) go to 160
          if(j1.eq.0) go to 200
          do 150 j2=1,j1
            j3 = j2+1
            h(j2) = h(j3)
 150      continue
          go to 180
 160      call fpgivs(piv,aa(jj,1),cos,sin)
          if(j1.eq.0) go to 200
          kk = jj
          do 170 j2=1,j1
            j3 = j2+1
            kk = kk-1
            call fprota(cos,sin,h(j3),aa(kk,j3))
            h(j2) = h(j3)
 170      continue
 180      jj = jj-1
          h(j3) = 0.
 190    continue
 200  continue
!c  solve the system (aa) (f1) = ff
      ff(rank) = ff(rank)/aa(rank,1)
      i = rank-1
      if(i.eq.0) go to 230
      do 220 j=2,rank
        store = ff(i)
        i1 = min0(j-1,m1)
        k = i
        do 210 ii=1,i1
          k = k+1
          stor1 = ff(k)
          stor2 = aa(i,ii+1)
          store = store-stor1*stor2
 210    continue
        stor1 = aa(i,1)
        ff(i) = store/stor1
        i = i-1
 220  continue
!c  solve the system  (aa)' (f2) = f1
 230  ff(1) = ff(1)/aa(1,1)
      if(rank.eq.1) go to 260
      do 250 j=2,rank
        store = ff(j)
        i1 = min0(j-1,m1)
        k = j
        do 240 ii=1,i1
          k = k-1
          stor1 = ff(k)
          stor2 = aa(k,ii+1)
          store = store-stor1*stor2
 240    continue
        stor1 = aa(j,1)
        ff(j) = store/stor1
 250  continue
!c  premultiply f2 by the transpoze of a.
 260  k = 0
      do 280 i=1,n
        store = 0.
        if(a(i,1).gt.tol) k = k+1
        j1 = min0(i,m)
        kk = k
        ij = i+1
        do 270 j=1,j1
          ij = ij-1
          if(a(ij,1).le.tol) go to 270
          stor1 = a(ij,j)
          stor2 = ff(kk)
          store = store+stor1*stor2
          kk = kk-1
 270    continue
        c(i) = store
 280  continue
!c  add to the sum of squared residuals the contribution of putting
!c  to zero the small diagonal elements of matrix (a).
      stor3 = 0.
      do 310 i=1,n
        if(a(i,1).gt.tol) go to 310
        store = f(i)
        i1 = min0(n-i,m1)
        if(i1.eq.0) go to 300
        do 290 j=1,i1
          ij = i+j
          stor1 = c(ij)
          stor2 = a(i,j+1)
          store = store-stor1*stor2
 290    continue
 300    fac = a(i,1)*c(i)
        stor1 = a(i,1)
        stor2 = c(i)
        stor1 = stor1*stor2
        stor3 = stor3+stor1*(stor1-store-store)
 310  continue
      fac = stor3
      sq = sq+fac
      return
      end

      real function fprati(p1,f1,p2,f2,p3,f3)
!c  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
!c  gives the value of p such that the rational interpolating function
!c  of the form r(p) = (u*p+v)/(p+w) equals zero at p.
!c  ..
!c  ..scalar arguments..
      real p1,f1,p2,f2,p3,f3
!c  ..local scalars..
      real h1,h2,h3,p
!c  ..
      if(p3.gt.0.) go to 10
!c  value of p in case p3 = infinity.
      p = (p1*(f1-f3)*f2-p2*(f2-f3)*f1)/((f1-f2)*f3)
      go to 20
!c  value of p in case p3 ^= infinity.
  10  h1 = f1*(f2-f3)
      h2 = f2*(f3-f1)
      h3 = f3*(f1-f2)
      p = -(p1*p2*h3+p2*p3*h1+p3*p1*h2)/(p1*h1+p2*h2+p3*h3)
!c  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
  20  if(f2.lt.0.) go to 30
      p1 = p2
      f1 = f2
      go to 40
  30  p3 = p2
      f3 = f2
  40  fprati = p
      return
      end
      subroutine fprota(cos,sin,a,b)
!c  subroutine fprota applies a givens rotation to a and b.
!c  ..
!c  ..scalar arguments..
      real cos,sin,a,b
!c ..local scalars..
      real stor1,stor2
!c  ..
      stor1 = a
      stor2 = b
      b = cos*stor2+sin*stor1
      a = cos*stor1-sin*stor2
      return
      end
      subroutine fprpsp(nt,np,co,si,c,f,ncoff)
!c  given the coefficients of a spherical spline function, subroutine
!c  fprpsp calculates the coefficients in the standard b-spline re-
!c  presentation of this bicubic spline.
!c  ..
!c  ..scalar arguments
      integer nt,np,ncoff
!c  ..array arguments
      real co(np),si(np),c(ncoff),f(ncoff)
!c  ..local scalars
      real cn,c1,c2,c3
      integer i,ii,j,k,l,ncof,npp,np4,nt4
!c  ..
      nt4 = nt-4
      np4 = np-4
      npp = np4-3
      ncof = 6+npp*(nt4-4)
      c1 = c(1)
      cn = c(ncof)
      j = ncoff
      do 10 i=1,np4
         f(i) = c1
         f(j) = cn
         j = j-1
  10  continue
      i = np4
      j=1
      do 70 l=3,nt4
         ii = i
         if(l.eq.3 .or. l.eq.nt4) go to 30
         do 20 k=1,npp
            i = i+1
            j = j+1
            f(i) = c(j)
  20     continue
         go to 50
  30     if(l.eq.nt4) c1 = cn
         c2 = c(j+1)
         c3 = c(j+2)
         j = j+2
         do 40 k=1,npp
            i = i+1
            f(i) = c1+c2*co(k)+c3*si(k)
  40     continue
  50     do 60 k=1,3
            ii = ii+1
            i = i+1
            f(i) = f(ii)
  60     continue
  70  continue
      do 80 i=1,ncoff
         c(i) = f(i)
  80  continue
      return
      end
      subroutine fpspgr(iopt,ider,u,mu,v,mv,r,mr,r0,r1,s,nuest,nvest,&
      tol,maxit,nc,nu,tu,nv,tv,c,fp,fp0,fpold,reducu,reducv,fpintu,&
      fpintv,dr,step,lastdi,nplusu,nplusv,lastu0,lastu1,nru,nrv,&
      nrdatu,nrdatv,wrk,lwrk,ier)
!c  ..
!c  ..scalar arguments..
      integer mu,mv,mr,nuest,nvest,maxit,nc,nu,nv,lastdi,nplusu,nplusv,&
      lastu0,lastu1,lwrk,ier
      real r0,r1,s,tol,fp,fp0,fpold,reducu,reducv
!c  ..array arguments..
      integer iopt(3),ider(4),nrdatu(nuest),nrdatv(nvest),nru(mu), nrv(mv)
      real u(mu),v(mv),r(mr),tu(nuest),tv(nvest),c(nc),fpintu(nuest),fpintv(nvest),dr(6),wrk(lwrk),step(2)
!c  ..local scalars..
      real acc,fpms,f1,f2,f3,p,per,pi,p1,p2,p3,vb,ve,rmax,rmin,rn,one,con1,con4,con9
      integer i,ich1,ich3,ifbu,ifbv,ifsu,ifsv,istart,iter,i1,i2,j,ju,&
      ktu,l,l1,l2,l3,l4,mpm,mumin,mu0,mu1,nn,nplu,nplv,npl1,nrintu,&
      nrintv,nue,numax,nve,nvmax
!c  ..local arrays..
      integer idd(4)
      real drr(6)
!c  ..function references..
      real abs,atan2,fprati
      integer max0,min0
!c  ..subroutine references..
!c    fpknot,fpopsp
!c  ..
!c   set constants
      one = 1
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
!c   initialization
      ifsu = 0
      ifsv = 0
      ifbu = 0
      ifbv = 0
      p = -one
      mumin = 4
      if(ider(1).ge.0) mumin = mumin-1
      if(iopt(2).eq.1 .and. ider(2).eq.1) mumin = mumin-1
      if(ider(3).ge.0) mumin = mumin-1
      if(iopt(3).eq.1 .and. ider(4).eq.1) mumin = mumin-1
      if(mumin.eq.0) mumin = 1
      pi = atan2(0.,-one)
      per = pi+pi
      vb = v(1)
      ve = vb+per
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c part 1: determination of the number of knots and their position.     c
!c ****************************************************************     c
!c  given a set of knots we compute the least-squares spline sinf(u,v)  c
!c  and the corresponding sum of squared residuals fp = f(p=inf).       c
!c  if iopt(1)=-1  sinf(u,v) is the requested approximation.            c
!c  if iopt(1)>=0  we check whether we can accept the knots:            c
!c    if fp <= s we will continue with the current set of knots.        c
!c    if fp >  s we will increase the number of knots and compute the   c
!c       corresponding least-squares spline until finally fp <= s.      c
!c    the initial choice of knots depends on the value of s and iopt.   c
!c    if s=0 we have spline interpolation; in that case the number of   c
!c     knots in the u-direction equals nu=numax=mu+6+iopt(2)+iopt(3)    c
!c     and in the v-direction nv=nvmax=mv+7.                            c
!c    if s>0 and                                                        c
!c      iopt(1)=0 we first compute the least-squares polynomial,i.e. a  c
!c       spline without interior knots : nu=8 ; nv=8.                   c
!c      iopt(1)=1 we start with the set of knots found at the last call c
!c       of the routine, except for the case that s > fp0; then we      c
!c       compute the least-squares polynomial directly.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(iopt(1).lt.0) go to 120
!c  acc denotes the absolute tolerance for the root of f(p)=s.
      acc = tol*s
!c  numax and nvmax denote the number of knots needed for interpolation.
      numax = mu+6+iopt(2)+iopt(3)
      nvmax = mv+7
      nue = min0(numax,nuest)
      nve = min0(nvmax,nvest)
      if(s.gt.0.) go to 100
!c  if s = 0, s(u,v) is an interpolating spline.
      nu = numax
      nv = nvmax
!c  test whether the required storage space exceeds the available one.
      if(nu.gt.nuest .or. nv.gt.nvest) go to 420
!c  find the position of the knots in the v-direction.
      do 10 l=1,mv
        tv(l+3) = v(l)
  10  continue
      tv(mv+4) = ve
      l1 = mv-2
      l2 = mv+5
      do 20 i=1,3
         tv(i) = v(l1)-per
         tv(l2) = v(i+1)+per
         l1 = l1+1
         l2 = l2+1
  20  continue
!c  if not all the derivative values g(i,j) are given, we will first
!c  estimate these values by computing a least-squares spline
      idd(1) = ider(1)
      if(idd(1).eq.0) idd(1) = 1
      if(idd(1).gt.0) dr(1) = r0
      idd(2) = ider(2)
      idd(3) = ider(3)
      if(idd(3).eq.0) idd(3) = 1
      if(idd(3).gt.0) dr(4) = r1
      idd(4) = ider(4)
      if(ider(1).lt.0 .or. ider(3).lt.0) go to 30
      if(iopt(2).ne.0 .and. ider(2).eq.0) go to 30
      if(iopt(3).eq.0 .or. ider(4).ne.0) go to 70
!c we set up the knots in the u-direction for computing the least-squares
!c spline.
  30  i1 = 3
      i2 = mu-2
      nu = 4
      do 40 i=1,mu
         if(i1.gt.i2) go to 50
         nu = nu+1
         tu(nu) = u(i1)
         i1 = i1+2
  40  continue
  50  do 60 i=1,4
         tu(i) = 0.
         nu = nu+1
         tu(nu) = pi
  60  continue
!c we compute the least-squares spline for estimating the derivatives.
      call fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,mr,r0,r1,dr,iopt,idd,&
       tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv,&
       wrk,lwrk)
      ifsu = 0
!c if all the derivatives at the origin are known, we compute the
!c interpolating spline.
!c we set up the knots in the u-direction, needed for interpolation.
  70  nn = numax-8
      if(nn.eq.0) go to 95
      ju = 2-iopt(2)
      do 80 l=1,nn
        tu(l+4) = u(ju)
        ju = ju+1
  80  continue
      nu = numax
      l = nu
      do 90 i=1,4
         tu(i) = 0.
         tu(l) = pi
         l = l-1
  90  continue
!c we compute the interpolating spline.
  95  call fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,mr,r0,r1,dr,iopt,idd,&
       tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv,&
       wrk,lwrk)
      go to 430
!c  if s>0 our initial choice of knots depends on the value of iopt(1).
 100  ier = 0
      if(iopt(1).eq.0) go to 115
      step(1) = -step(1)
      step(2) = -step(2)
      if(fp0.le.s) go to 115
!c  if iopt(1)=1 and fp0 > s we start computing the least-squares spline
!c  according to the set of knots found at the last call of the routine.
!c  we determine the number of grid coordinates u(i) inside each knot
!c  interval (tu(l),tu(l+1)).
      l = 5
      j = 1
      nrdatu(1) = 0
      mu0 = 2-iopt(2)
      mu1 = mu-1+iopt(3)
      do 105 i=mu0,mu1
        nrdatu(j) = nrdatu(j)+1
        if(u(i).lt.tu(l)) go to 105
        nrdatu(j) = nrdatu(j)-1
        l = l+1
        j = j+1
        nrdatu(j) = 0
 105  continue
!c  we determine the number of grid coordinates v(i) inside each knot
!c  interval (tv(l),tv(l+1)).
      l = 5
      j = 1
      nrdatv(1) = 0
      do 110 i=2,mv
        nrdatv(j) = nrdatv(j)+1
        if(v(i).lt.tv(l)) go to 110
        nrdatv(j) = nrdatv(j)-1
        l = l+1
        j = j+1
        nrdatv(j) = 0
 110  continue
      idd(1) = ider(1)
      idd(2) = ider(2)
      idd(3) = ider(3)
      idd(4) = ider(4)
      go to 120
!c  if iopt(1)=0 or iopt(1)=1 and s >= fp0,we start computing the least-
!c  squares polynomial (which is a spline without interior knots).
 115  ier = -2
      idd(1) = ider(1)
      idd(2) = 1
      idd(3) = ider(3)
      idd(4) = 1
      nu = 8
      nv = 8
      nrdatu(1) = mu-2+iopt(2)+iopt(3)
      nrdatv(1) = mv-1
      lastdi = 0
      nplusu = 0
      nplusv = 0
      fp0 = 0.
      fpold = 0.
      reducu = 0.
      reducv = 0.
!c  main loop for the different sets of knots.mpm=mu+mv is a save upper
!c  bound for the number of trials.
 120  mpm = mu+mv
      do 270 iter=1,mpm
!c  find nrintu (nrintv) which is the number of knot intervals in the
!c  u-direction (v-direction).
        nrintu = nu-7
        nrintv = nv-7
!c  find the position of the additional knots which are needed for the
!c  b-spline representation of s(u,v).
        i = nu
        do 125 j=1,4
          tu(j) = 0.
          tu(i) = pi
          i = i-1
 125    continue
        l1 = 4
        l2 = l1
        l3 = nv-3
        l4 = l3
        tv(l2) = vb
        tv(l3) = ve
        do 130 j=1,3
          l1 = l1+1
          l2 = l2-1
          l3 = l3+1
          l4 = l4-1
          tv(l2) = tv(l4)-per
          tv(l3) = tv(l1)+per
 130    continue
!c  find an estimate of the range of possible values for the optimal
!c  derivatives at the origin.
        ktu = nrdatu(1)+2-iopt(2)
        if(ktu.lt.mumin) ktu = mumin
        if(ktu.eq.lastu0) go to 140
         rmin = r0
         rmax = r0
         l = mv*ktu
         do 135 i=1,l
            if(r(i).lt.rmin) rmin = r(i)
            if(r(i).gt.rmax) rmax = r(i)
 135     continue
         step(1) = rmax-rmin
         lastu0 = ktu
 140    ktu = nrdatu(nrintu)+2-iopt(3)
        if(ktu.lt.mumin) ktu = mumin
        if(ktu.eq.lastu1) go to 150
         rmin = r1
         rmax = r1
         l = mv*ktu
         j = mr
         do 145 i=1,l
            if(r(j).lt.rmin) rmin = r(j)
            if(r(j).gt.rmax) rmax = r(j)
            j = j-1
 145     continue
         step(2) = rmax-rmin
         lastu1 = ktu
!c  find the least-squares spline sinf(u,v).
 150    call fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,mr,r0,r1,dr,iopt,&
        idd,tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,&
        nrv,wrk,lwrk)
        if(step(1).lt.0.) step(1) = -step(1)
        if(step(2).lt.0.) step(2) = -step(2)
        if(ier.eq.(-2)) fp0 = fp
!c  test whether the least-squares spline is an acceptable solution.
        if(iopt(1).lt.0) go to 440
        fpms = fp-s
        if(abs(fpms) .lt. acc) go to 440
!c  if f(p=inf) < s, we accept the choice of knots.
        if(fpms.lt.0.) go to 300
!c  if nu=numax and nv=nvmax, sinf(u,v) is an interpolating spline
        if(nu.eq.numax .and. nv.eq.nvmax) go to 430
!c  increase the number of knots.
!c  if nu=nue and nv=nve we cannot further increase the number of knots
!c  because of the storage capacity limitation.
        if(nu.eq.nue .and. nv.eq.nve) go to 420
        if(ider(1).eq.0) fpintu(1) = fpintu(1)+(r0-dr(1))**2
        if(ider(3).eq.0) fpintu(nrintu) = fpintu(nrintu)+(r1-dr(4))**2
        ier = 0
!c  adjust the parameter reducu or reducv according to the direction
!c  in which the last added knots were located.
        if(lastdi) 160,155,170
 155     nplv = 3
         idd(2) = ider(2)
         idd(4) = ider(4)
         fpold = fp
         go to 230
 160    reducu = fpold-fp
        go to 175
 170    reducv = fpold-fp
!c  store the sum of squared residuals for the current set of knots.
 175    fpold = fp
!c  find nplu, the number of knots we should add in the u-direction.
        nplu = 1
        if(nu.eq.8) go to 180
        npl1 = nplusu*2
        rn = nplusu
        if(reducu.gt.acc) npl1 = rn*fpms/reducu
        nplu = min0(nplusu*2,max0(npl1,nplusu/2,1))
!c  find nplv, the number of knots we should add in the v-direction.
 180    nplv = 3
        if(nv.eq.8) go to 190
        npl1 = nplusv*2
        rn = nplusv
        if(reducv.gt.acc) npl1 = rn*fpms/reducv
        nplv = min0(nplusv*2,max0(npl1,nplusv/2,1))
!c  test whether we are going to add knots in the u- or v-direction.
 190    if(nplu-nplv) 210,200,230
 200    if(lastdi.lt.0) go to 230
 210    if(nu.eq.nue) go to 230
!c  addition in the u-direction.
        lastdi = -1
        nplusu = nplu
        ifsu = 0
        istart = 0
        if(iopt(2).eq.0) istart = 1
        do 220 l=1,nplusu
!c  add a new knot in the u-direction
          call fpknot(u,mu,tu,nu,fpintu,nrdatu,nrintu,nuest,istart)
!c  test whether we cannot further increase the number of knots in the
!c  u-direction.
          if(nu.eq.nue) go to 270
 220    continue
        go to 270
 230    if(nv.eq.nve) go to 210
!c  addition in the v-direction.
        lastdi = 1
        nplusv = nplv
        ifsv = 0
        do 240 l=1,nplusv
!c  add a new knot in the v-direction.
          call fpknot(v,mv,tv,nv,fpintv,nrdatv,nrintv,nvest,1)
!c  test whether we cannot further increase the number of knots in the
!c  v-direction.
          if(nv.eq.nve) go to 270
 240    continue
!c  restart the computations with the new set of knots.
 270  continue
!c  test whether the least-squares polynomial is a solution of our
!c  approximation problem.
 300  if(ier.eq.(-2)) go to 440
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c part 2: determination of the smoothing spline sp(u,v)                c
!c *****************************************************                c
!c  we have determined the number of knots and their position. we now   c
!c  compute the b-spline coefficients of the smoothing spline sp(u,v).  c
!c  this smoothing spline depends on the parameter p in such a way that c
!c    f(p) = sumi=1,mu(sumj=1,mv((z(i,j)-sp(u(i),v(j)))**2)             c
!c  is a continuous, strictly decreasing function of p. moreover the    c
!c  least-squares polynomial corresponds to p=0 and the least-squares   c
!c  spline to p=infinity. then iteratively we have to determine the     c
!c  positive value of p such that f(p)=s. the process which is proposed c
!c  here makes use of rational interpolation. f(p) is approximated by a c
!c  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  c
!c  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)c
!c  are used to calculate the new value of p such that r(p)=s.          c
!c  convergence is guaranteed by taking f1 > 0 and f3 < 0.              c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = one
      do 305 i=1,6
        drr(i) = dr(i)
 305  continue
      ich1 = 0
      ich3 = 0
!c  iteration process to find the root of f(p)=s.
      do 350 iter = 1,maxit
!c  find the smoothing spline sp(u,v) and the corresponding sum f(p).
        call fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,mr,r0,r1,drr,iopt,&
        idd,tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,&
        nrv,wrk,lwrk)
!c  test whether the approximation sp(u,v) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
!c  test whether the maximum allowable number of iterations has been
!c  reached.
        if(iter.eq.maxit) go to 400
!c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 320
        if((f2-f3).gt.acc) go to 310
!c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 + p2*con1
        go to 350
 310    if(f2.lt.0.) ich3 = 1
 320    if(ich1.ne.0) go to 340
        if((f1-f2).gt.acc) go to 330
!c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 350
        if(p.ge.p3) p = p2*con1 + p3*con9
        go to 350
!c  test whether the iteration process proceeds as theoretically
!c  expected.
 330    if(f2.gt.0.) ich1 = 1
 340    if(f2.ge.f1 .or. f2.le.f3) go to 410
!c  find the new value of p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 350  continue
!c  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
      fp = 0.
 440  return
      end
      subroutine fpsphe(iopt,m,teta,phi,r,w,s,ntest,npest,eta,tol,maxit,&
      ib1,ib3,nc,ncc,intest,nrest,nt,tt,np,tp,c,fp,sup,fpint,coord,f,&
      ff,row,coco,cosi,a,q,bt,bp,spt,spp,h,index,nummer,wrk,lwrk,ier)
!c  ..
!c  ..scalar arguments..
      integer iopt,m,ntest,npest,maxit,ib1,ib3,nc,ncc,intest,nrest, nt,np,lwrk,ier
      real s,eta,tol,fp,sup
!c  ..array arguments..
      real teta(m),phi(m),r(m),w(m),tt(ntest),tp(npest),c(nc),&
      fpint(intest),coord(intest),f(ncc),ff(nc),row(npest),coco(npest),&
      cosi(npest),a(ncc,ib1),q(ncc,ib3),bt(ntest,5),bp(npest,5),&
      spt(m,4),spp(m,4),h(ib3),wrk(lwrk)
      integer index(nrest),nummer(m)
!c  ..local scalars..
      real aa,acc,arg,cn,co,c1,dmax,d1,d2,eps,facc,facs,fac1,fac2,fn,&
      fpmax,fpms,f1,f2,f3,hti,htj,p,pi,pinv,piv,pi2,p1,p2,p3,ri,si,&
      sigma,sq,store,wi,rn,one,con1,con9,con4,half,ten
      integer i,iband,iband1,iband3,iband4,ich1,ich3,ii,ij,il,in,irot,&
      iter,i1,i2,i3,j,jlt,jrot,j1,j2,l,la,lf,lh,ll,lp,lt,lwest,l1,l2,&
      l3,l4,ncof,ncoff,npp,np4,nreg,nrint,nrr,nr1,ntt,nt4,nt6,num,&
      num1,rank
!c  ..local arrays..
      real ht(4),hp(4)
!c  ..function references..
      real abs,atan,fprati,sqrt,cos,sin
      integer min0
!c  ..subroutine references..
!c   fpback,fpbspl,fpgivs,fpdisc,fporde,fprank,fprota,fprpsp
!c  ..
!c  set constants
      one = 0.1e+01
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
      half = 0.5e0
      ten = 0.1e+02
      pi = atan(one)*4
      pi2 = pi+pi
      eps = sqrt(eta)
      if(iopt.lt.0) go to 70
!c  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
      if(iopt.eq.0) go to 10
      if(s.lt.sup) if(np-11) 60,70,70
!c  if iopt=0 we begin by computing the weighted least-squares polynomial
!c  of the form
!c     s(teta,phi) = c1*f1(teta) + cn*fn(teta)
!c  where f1(teta) and fn(teta) are the cubic polynomials satisfying
!c     f1(0) = 1, f1(pi) = f1'(0) = f1'(pi) = 0 ; fn(teta) = 1-f1(teta).
!c  the corresponding weighted sum of squared residuals gives the upper
!c  bound sup for the smoothing factor s.
  10  sup = 0.
      d1 = 0.
      d2 = 0.
      c1 = 0.
      cn = 0.
      fac1 = pi*(one + half)
      fac2 = (one + one)/pi**3
      aa = 0.
      do 40 i=1,m
         wi = w(i)
         ri = r(i)*wi
         arg = teta(i)
         fn = fac2*arg*arg*(fac1-arg)
         f1 = (one-fn)*wi
         fn = fn*wi
         if(fn.eq.0.) go to 20
         call fpgivs(fn,d1,co,si)
         call fprota(co,si,f1,aa)
         call fprota(co,si,ri,cn)
 20      if(f1.eq.0.) go to 30
         call fpgivs(f1,d2,co,si)
         call fprota(co,si,ri,c1)
 30      sup = sup+ri*ri
 40   continue
      if(d2.ne.0.) c1 = c1/d2
      if(d1.ne.0.) cn = (cn-aa*c1)/d1
!c  find the b-spline representation of this least-squares polynomial
      nt = 8
      np = 8
      do 50 i=1,4
         c(i) = c1
         c(i+4) = c1
         c(i+8) = cn
         c(i+12) = cn
         tt(i) = 0.
         tt(i+4) = pi
         tp(i) = 0.
         tp(i+4) = pi2
  50  continue
      fp = sup
!c  test whether the least-squares polynomial is an acceptable solution
      fpms = sup-s
      if(fpms.lt.acc) go to 960
!c  test whether we cannot further increase the number of knots.
  60  if(npest.lt.11 .or. ntest.lt.9) go to 950
!c  find the initial set of interior knots of the spherical spline in
!c  case iopt = 0.
      np = 11
      tp(5) = pi*half
      tp(6) = pi
      tp(7) = tp(5)+pi
      nt = 9
      tt(5) = tp(5)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  part 1 : computation of least-squares spherical splines.            c
!c  ********************************************************            c
!c  if iopt < 0 we compute the least-squares spherical spline according c
!c  to the given set of knots.                                          c
!c  if iopt >=0 we compute least-squares spherical splines with increas-c
!c  ing numbers of knots until the corresponding sum f(p=inf)<=s.       c
!c  the initial set of knots then depends on the value of iopt:         c
!c    if iopt=0 we start with one interior knot in the teta-direction   c
!c              (pi/2) and three in the phi-direction (pi/2,pi,3*pi/2). c
!c    if iopt>0 we start with the set of knots found at the last call   c
!c              of the routine.                                         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  main loop for the different sets of knots. m is a save upper bound
!c  for the number of trials.
  70  do 570 iter=1,m
!c  find the position of the additional knots which are needed for the
!c  b-spline representation of s(teta,phi).
         l1 = 4
         l2 = l1
         l3 = np-3
         l4 = l3
         tp(l2) = 0.
         tp(l3) = pi2
         do 80 i=1,3
            l1 = l1+1
            l2 = l2-1
            l3 = l3+1
            l4 = l4-1
            tp(l2) = tp(l4)-pi2
            tp(l3) = tp(l1)+pi2
  80     continue
        l = nt
        do 90 i=1,4
          tt(i) = 0.
          tt(l) = pi
          l = l-1
  90    continue
!c  find nrint, the total number of knot intervals and nreg, the number
!c  of panels in which the approximation domain is subdivided by the
!c  intersection of knots.
        ntt = nt-7
        npp = np-7
        nrr = npp/2
        nr1 = nrr+1
        nrint = ntt+npp
        nreg = ntt*npp
!c  arrange the data points according to the panel they belong to.
        call fporde(teta,phi,m,3,3,tt,nt,tp,np,nummer,index,nreg)
!c  find the b-spline coefficients coco and cosi of the cubic spline
!c  approximations sc(phi) and ss(phi) for cos(phi) and sin(phi).
        do 100 i=1,npp
           coco(i) = 0.
           cosi(i) = 0.
           do 100 j=1,npp
              a(i,j) = 0.
 100    continue
!c  the coefficients coco and cosi are obtained from the conditions
!c  sc(tp(i))=cos(tp(i)),resp. ss(tp(i))=sin(tp(i)),i=4,5,...np-4.
        do 150 i=1,npp
           l2 = i+3
           arg = tp(l2)
           call fpbspl(tp,np,3,arg,l2,hp)
           do 110 j=1,npp
              row(j) = 0.
 110       continue
           ll = i
           do 120 j=1,3
              if(ll.gt.npp) ll= 1
              row(ll) = row(ll)+hp(j)
              ll = ll+1
 120       continue
           facc = cos(arg)
           facs = sin(arg)
           do 140 j=1,npp
              piv = row(j)
              if(piv.eq.0.) go to 140
              call fpgivs(piv,a(j,1),co,si)
              call fprota(co,si,facc,coco(j))
              call fprota(co,si,facs,cosi(j))
              if(j.eq.npp) go to 150
              j1 = j+1
              i2 = 1
              do 130 l=j1,npp
                 i2 = i2+1
                 call fprota(co,si,row(l),a(j,i2))
 130          continue
 140       continue
 150    continue
        call fpback(a,coco,npp,npp,coco,ncc)
        call fpback(a,cosi,npp,npp,cosi,ncc)
!c  find ncof, the dimension of the spherical spline and ncoff, the
!c  number of coefficients in the standard b-spline representation.
        nt4 = nt-4
        np4 = np-4
        ncoff = nt4*np4
        ncof = 6+npp*(ntt-1)
!c  find the bandwidth of the observation matrix a.
        iband = 4*npp
        if(ntt.eq.4) iband = 3*(npp+1)
        if(ntt.lt.4) iband = ncof
        iband1 = iband-1
!c  initialize the observation matrix a.
        do 160 i=1,ncof
          f(i) = 0.
          do 160 j=1,iband
            a(i,j) = 0.
 160    continue
!c  initialize the sum of squared residuals.
        fp = 0.
!c  fetch the data points in the new order. main loop for the
!c  different panels.
        do 340 num=1,nreg
!c  fix certain constants for the current panel; jrot records the column
!c  number of the first non-zero element in a row of the observation
!c  matrix according to a data point of the panel.
          num1 = num-1
          lt = num1/npp
          l1 = lt+4
          lp = num1-lt*npp+1
          l2 = lp+3
          lt = lt+1
          jrot = 0
          if(lt.gt.2) jrot = 3+(lt-3)*npp
!c  test whether there are still data points in the current panel.
          in = index(num)
 170      if(in.eq.0) go to 340
!c  fetch a new data point.
          wi = w(in)
          ri = r(in)*wi
!c  evaluate for the teta-direction, the 4 non-zero b-splines at teta(in)
          call fpbspl(tt,nt,3,teta(in),l1,ht)
!c  evaluate for the phi-direction, the 4 non-zero b-splines at phi(in)
          call fpbspl(tp,np,3,phi(in),l2,hp)
!c  store the value of these b-splines in spt and spp resp.
          do 180 i=1,4
            spp(in,i) = hp(i)
            spt(in,i) = ht(i)
 180      continue
!c  initialize the new row of observation matrix.
          do 190 i=1,iband
            h(i) = 0.
 190      continue
!c  calculate the non-zero elements of the new row by making the cross
!c  products of the non-zero b-splines in teta- and phi-direction and
!c  by taking into account the conditions of the spherical splines.
          do 200 i=1,npp
             row(i) = 0.
 200      continue
!c  take into account the condition (3) of the spherical splines.
          ll = lp
          do 210 i=1,4
             if(ll.gt.npp) ll=1
             row(ll) = row(ll)+hp(i)
             ll = ll+1
 210      continue
!c  take into account the other conditions of the spherical splines.
          if(lt.gt.2 .and. lt.lt.(ntt-1)) go to 230
          facc = 0.
          facs = 0.
          do 220 i=1,npp
             facc = facc+row(i)*coco(i)
             facs = facs+row(i)*cosi(i)
 220     continue
!c  fill in the non-zero elements of the new row.
 230     j1 = 0
         do 280 j =1,4
            jlt = j+lt
            htj = ht(j)
            if(jlt.gt.2 .and. jlt.le.nt4) go to 240
            j1 = j1+1
            h(j1) = h(j1)+htj
            go to 280
 240        if(jlt.eq.3 .or. jlt.eq.nt4) go to 260
            do 250 i=1,npp
               j1 = j1+1
               h(j1) = row(i)*htj
 250        continue
            go to 280
 260        if(jlt.eq.3) go to 270
            h(j1+1) = facc*htj
            h(j1+2) = facs*htj
            h(j1+3) = htj
            j1 = j1+2
            go to 280
 270        h(1) = h(1)+htj
            h(2) = facc*htj
            h(3) = facs*htj
            j1 = 3
 280      continue
          do 290 i=1,iband
            h(i) = h(i)*wi
 290      continue
!c  rotate the row into triangle by givens transformations.
          irot = jrot
          do 310 i=1,iband
            irot = irot+1
            piv = h(i)
            if(piv.eq.0.) go to 310
!c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(irot,1),co,si)
!c  apply that transformation to the right hand side.
            call fprota(co,si,ri,f(irot))
            if(i.eq.iband) go to 320
!c  apply that transformation to the left hand side.
            i2 = 1
            i3 = i+1
            do 300 j=i3,iband
              i2 = i2+1
              call fprota(co,si,h(j),a(irot,i2))
 300        continue
 310      continue
!c  add the contribution of the row to the sum of squares of residual
!c  right hand sides.
 320      fp = fp+ri**2
!c  find the number of the next data point in the panel.
 330      in = nummer(in)
          go to 170
 340    continue
!c  find dmax, the maximum value for the diagonal elements in the reduced
!c  triangle.
        dmax = 0.
        do 350 i=1,ncof
          if(a(i,1).le.dmax) go to 350
          dmax = a(i,1)
 350    continue
!c  check whether the observation matrix is rank deficient.
        sigma = eps*dmax
        do 360 i=1,ncof
          if(a(i,1).le.sigma) go to 370
 360    continue
!c  backward substitution in case of full rank.
        call fpback(a,f,ncof,iband,c,ncc)
        rank = ncof
        do 365 i=1,ncof
          q(i,1) = a(i,1)/dmax
 365    continue
        go to 390
!c  in case of rank deficiency, find the minimum norm solution.
 370    lwest = ncof*iband+ncof+iband
        if(lwrk.lt.lwest) go to 925
        lf = 1
        lh = lf+ncof
        la = lh+iband
        do 380 i=1,ncof
          ff(i) = f(i)
          do 380 j=1,iband
            q(i,j) = a(i,j)
 380    continue
        call fprank(q,ff,ncof,iband,ncc,sigma,c,sq,rank,wrk(la), wrk(lf),wrk(lh))
        do 385 i=1,ncof
          q(i,1) = q(i,1)/dmax
 385    continue
!c  add to the sum of squared residuals, the contribution of reducing
!c  the rank.
        fp = fp+sq
!c  find the coefficients in the standard b-spline representation of
!c  the spherical spline.
 390    call fprpsp(nt,np,coco,cosi,c,ff,ncoff)
!c  test whether the least-squares spline is an acceptable solution.
        if(iopt.lt.0) if(fp) 970,970,980
        fpms = fp-s
        if(abs(fpms).le.acc) if(fp) 970,970,980
!c  if f(p=inf) < s, accept the choice of knots.
        if(fpms.lt.0.) go to 580
!c  test whether we cannot further increase the number of knots.
        if(ncof.gt.m) go to 935
!c  search where to add a new knot.
!c  find for each interval the sum of squared residuals fpint for the
!c  data points having the coordinate belonging to that knot interval.
!c  calculate also coord which is the same sum, weighted by the position
!c  of the data points considered.
 440    do 450 i=1,nrint
          fpint(i) = 0.
          coord(i) = 0.
 450    continue
        do 490 num=1,nreg
          num1 = num-1
          lt = num1/npp
          l1 = lt+1
          lp = num1-lt*npp
          l2 = lp+1+ntt
          jrot = lt*np4+lp
          in = index(num)
 460      if(in.eq.0) go to 490
          store = 0.
          i1 = jrot
          do 480 i=1,4
            hti = spt(in,i)
            j1 = i1
            do 470 j=1,4
              j1 = j1+1
              store = store+hti*spp(in,j)*c(j1)
 470        continue
            i1 = i1+np4
 480      continue
          store = (w(in)*(r(in)-store))**2
          fpint(l1) = fpint(l1)+store
          coord(l1) = coord(l1)+store*teta(in)
          fpint(l2) = fpint(l2)+store
          coord(l2) = coord(l2)+store*phi(in)
          in = nummer(in)
          go to 460
 490    continue
!c  find the interval for which fpint is maximal on the condition that
!c  there still can be added a knot.
        l1 = 1
        l2 = nrint
        if(ntest.lt.nt+1) l1=ntt+1
        if(npest.lt.np+2) l2=ntt
!c  test whether we cannot further increase the number of knots.
        if(l1.gt.l2) go to 950
 500    fpmax = 0.
        l = 0
        do 510 i=l1,l2
          if(fpmax.ge.fpint(i)) go to 510
          l = i
          fpmax = fpint(i)
 510    continue
        if(l.eq.0) go to 930
!c  calculate the position of the new knot.
        arg = coord(l)/fpint(l)
!c  test in what direction the new knot is going to be added.
        if(l.gt.ntt) go to 530
!c  addition in the teta-direction
        l4 = l+4
        fpint(l) = 0.
        fac1 = tt(l4)-arg
        fac2 = arg-tt(l4-1)
        if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 500
        j = nt
        do 520 i=l4,nt
          tt(j+1) = tt(j)
          j = j-1
 520    continue
        tt(l4) = arg
        nt = nt+1
        go to 570
!c  addition in the phi-direction
 530    l4 = l+4-ntt
        if(arg.lt.pi) go to 540
        arg = arg-pi
        l4 = l4-nrr
 540    fpint(l) = 0.
        fac1 = tp(l4)-arg
        fac2 = arg-tp(l4-1)
        if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 500
        ll = nrr+4
        j = ll
        do 550 i=l4,ll
          tp(j+1) = tp(j)
          j = j-1
 550    continue
        tp(l4) = arg
        np = np+2
        nrr = nrr+1
        do 560 i=5,ll
          j = i+nrr
          tp(j) = tp(i)+pi
 560    continue
!c  restart the computations with the new set of knots.
 570  continue
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c part 2: determination of the smoothing spherical spline.             c
!c ********************************************************             c
!c we have determined the number of knots and their position. we now    c
!c compute the coefficients of the smoothing spline sp(teta,phi).       c
!c the observation matrix a is extended by the rows of a matrix, expres-c
!c sing that sp(teta,phi) must be a constant function in the variable   c
!c phi and a cubic polynomial in the variable teta. the corresponding   c
!c weights of these additional rows are set to 1/(p). iteratively       c
!c we than have to determine the value of p such that f(p) = sum((w(i)* c
!c (r(i)-sp(teta(i),phi(i))))**2)  be = s.                              c
!c we already know that the least-squares polynomial corresponds to p=0,c
!c and that the least-squares spherical spline corresponds to p=infin.  c
!c the iteration process makes use of rational interpolation. since f(p)c
!c is a convex and strictly decreasing function of p, it can be approx- c
!c imated by a rational function of the form r(p) = (u*p+v)/(p+w).      c
!c three values of p (p1,p2,p3) with corresponding values of f(p) (f1=  c
!c f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the new value   c
!c of p such that r(p)=s. convergence is guaranteed by taking f1>0,f3<0.c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  evaluate the discontinuity jumps of the 3-th order derivative of
!c  the b-splines at the knots tt(l),l=5,...,nt-4.
 580  call fpdisc(tt,nt,5,bt,ntest)
!c  evaluate the discontinuity jumps of the 3-th order derivative of
!c  the b-splines at the knots tp(l),l=5,...,np-4.
      call fpdisc(tp,np,5,bp,npest)
!c  initial value for p.
      p1 = 0.
      f1 = sup-s
      p3 = -one
      f3 = fpms
      p = 0.
      do 585 i=1,ncof
        p = p+a(i,1)
 585  continue
      rn = ncof
      p = rn/p
!c  find the bandwidth of the extended observation matrix.
      iband4 = iband+3
      if(ntt.le.4) iband4 = ncof
      iband3 = iband4 -1
      ich1 = 0
      ich3 = 0
!c  iteration process to find the root of f(p)=s.
      do 920 iter=1,maxit
        pinv = one/p
!c  store the triangularized observation matrix into q.
        do 600 i=1,ncof
          ff(i) = f(i)
          do 590 j=1,iband4
            q(i,j) = 0.
 590      continue
          do 600 j=1,iband
            q(i,j) = a(i,j)
 600    continue
!c  extend the observation matrix with the rows of a matrix, expressing
!c  that for teta=cst. sp(teta,phi) must be a constant function.
        nt6 = nt-6
        do 720 i=5,np4
          ii = i-4
          do 610 l=1,npp
             row(l) = 0.
 610      continue
          ll = ii
          do 620  l=1,5
             if(ll.gt.npp) ll=1
             row(ll) = row(ll)+bp(ii,l)
             ll = ll+1
 620      continue
          facc = 0.
          facs = 0.
          do 630 l=1,npp
             facc = facc+row(l)*coco(l)
             facs = facs+row(l)*cosi(l)
 630      continue
          do 720 j=1,nt6
!c  initialize the new row.
            do 640 l=1,iband
              h(l) = 0.
 640        continue
!c  fill in the non-zero elements of the row. jrot records the column
!c  number of the first non-zero element in the row.
            jrot = 4+(j-2)*npp
            if(j.gt.1 .and. j.lt.nt6) go to 650
            h(1) = facc
            h(2) = facs
            if(j.eq.1) jrot = 2
            go to 670
 650        do 660 l=1,npp
               h(l)=row(l)
 660        continue
 670        do 675 l=1,iband
               h(l) = h(l)*pinv
 675        continue
            ri = 0.
!c  rotate the new row into triangle by givens transformations.
            do 710 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband1,ncof-irot)
              if(piv.eq.0.) if(i2) 720,720,690
!c  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),co,si)
!c  apply that givens transformation to the right hand side.
              call fprota(co,si,ri,ff(irot))
              if(i2.eq.0) go to 720
!c  apply that givens transformation to the left hand side.
              do 680 l=1,i2
                l1 = l+1
                call fprota(co,si,h(l1),q(irot,l1))
 680          continue
 690          do 700 l=1,i2
                h(l) = h(l+1)
 700          continue
              h(i2+1) = 0.
 710        continue
 720    continue
!c  extend the observation matrix with the rows of a matrix expressing
!c  that for phi=cst. sp(teta,phi) must be a cubic polynomial.
        do 810 i=5,nt4
          ii = i-4
          do 810 j=1,npp
!c  initialize the new row
            do 730 l=1,iband4
              h(l) = 0.
 730        continue
!c  fill in the non-zero elements of the row. jrot records the column
!c  number of the first non-zero element in the row.
            j1 = 1
            do 760 l=1,5
               il = ii+l
               ij = npp
               if(il.ne.3 .and. il.ne.nt4) go to 750
               j1 = j1+3-j
               j2 = j1-2
               ij = 0
               if(il.ne.3) go to 740
               j1 = 1
               j2 = 2
               ij = j+2
 740           h(j2) = bt(ii,l)*coco(j)
               h(j2+1) = bt(ii,l)*cosi(j)
 750           h(j1) = h(j1)+bt(ii,l)
               j1 = j1+ij
 760        continue
            do 765 l=1,iband4
               h(l) = h(l)*pinv
 765        continue
            ri = 0.
            jrot = 1
            if(ii.gt.2) jrot = 3+j+(ii-3)*npp
!c  rotate the new row into triangle by givens transformations.
            do 800 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband3,ncof-irot)
              if(piv.eq.0.) if(i2) 810,810,780
!c  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),co,si)
!c  apply that givens transformation to the right hand side.
              call fprota(co,si,ri,ff(irot))
              if(i2.eq.0) go to 810
!c  apply that givens transformation to the left hand side.
              do 770 l=1,i2
                l1 = l+1
                call fprota(co,si,h(l1),q(irot,l1))
 770          continue
 780          do 790 l=1,i2
                h(l) = h(l+1)
 790          continue
              h(i2+1) = 0.
 800        continue
 810    continue
!c  find dmax, the maximum value for the diagonal elements in the
!c  reduced triangle.
        dmax = 0.
        do 820 i=1,ncof
          if(q(i,1).le.dmax) go to 820
          dmax = q(i,1)
 820    continue
!c  check whether the matrix is rank deficient.
        sigma = eps*dmax
        do 830 i=1,ncof
          if(q(i,1).le.sigma) go to 840
 830    continue
!c  backward substitution in case of full rank.
        call fpback(q,ff,ncof,iband4,c,ncc)
        rank = ncof
        go to 845
!c  in case of rank deficiency, find the minimum norm solution.
 840    lwest = ncof*iband4+ncof+iband4
        if(lwrk.lt.lwest) go to 925
        lf = 1
        lh = lf+ncof
        la = lh+iband4
        call fprank(q,ff,ncof,iband4,ncc,sigma,c,sq,rank,wrk(la), wrk(lf),wrk(lh))
 845    do 850 i=1,ncof
           q(i,1) = q(i,1)/dmax
 850    continue
!c  find the coefficients in the standard b-spline representation of
!c  the spherical spline.
        call fprpsp(nt,np,coco,cosi,c,ff,ncoff)
!c  compute f(p).
        fp = 0.
        do 890 num = 1,nreg
          num1 = num-1
          lt = num1/npp
          lp = num1-lt*npp
          jrot = lt*np4+lp
          in = index(num)
 860      if(in.eq.0) go to 890
          store = 0.
          i1 = jrot
          do 880 i=1,4
            hti = spt(in,i)
            j1 = i1
            do 870 j=1,4
              j1 = j1+1
              store = store+hti*spp(in,j)*c(j1)
 870        continue
            i1 = i1+np4
 880      continue
          fp = fp+(w(in)*(r(in)-store))**2
          in = nummer(in)
          go to 860
 890    continue
!c  test whether the approximation sp(teta,phi) is an acceptable solution
        fpms = fp-s
        if(abs(fpms).le.acc) go to 980
!c  test whether the maximum allowable number of iterations has been
!c  reached.
        if(iter.eq.maxit) go to 940
!c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 900
        if((f2-f3).gt.acc) go to 895
!c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 + p2*con1
        go to 920
 895    if(f2.lt.0.) ich3 = 1
 900    if(ich1.ne.0) go to 910
        if((f1-f2).gt.acc) go to 905
!c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 920
        if(p.ge.p3) p = p2*con1 +p3*con9
        go to 920
 905    if(f2.gt.0.) ich1 = 1
!c  test whether the iteration process proceeds as theoretically
!c  expected.
 910    if(f2.ge.f1 .or. f2.le.f3) go to 945
!c  find the new value of p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 920  continue
!c  error codes and messages.
 925  ier = lwest
      go to 990
 930  ier = 5
      go to 990
 935  ier = 4
      go to 990
 940  ier = 3
      go to 990
 945  ier = 2
      go to 990
 950  ier = 1
      go to 990
 960  ier = -2
      go to 990
 970  ier = -1
      fp = 0.
 980  if(ncof.ne.rank) ier = -rank
 990  return
      end
      subroutine fpsysy(a,n,g)
!c subroutine fpsysy solves a linear n x n symmetric system
!c    (a) * (b) = (g)
!c on input, vector g contains the right hand side ; on output it will
!c contain the solution (b).
!c  ..
!c  ..scalar arguments..
      integer n
!c  ..array arguments..
      real a(6,6),g(6)
!c  ..local scalars..
      real fac
      integer i,i1,j,k
!c  ..
      g(1) = g(1)/a(1,1)
      if(n.eq.1) return
!c  decomposition of the symmetric matrix (a) = (l) * (d) *(l)'
!c  with (l) a unit lower triangular matrix and (d) a diagonal
!c  matrix
      do 10 k=2,n
         a(k,1) = a(k,1)/a(1,1)
  10  continue
      do 40 i=2,n
         i1 = i-1
         do 30 k=i,n
            fac = a(k,i)
            do 20 j=1,i1
               fac = fac-a(j,j)*a(k,j)*a(i,j)
  20        continue
            a(k,i) = fac
            if(k.gt.i) a(k,i) = fac/a(i,i)
  30     continue
  40  continue
!c  solve the system (l)*(d)*(l)'*(b) = (g).
!c  first step : solve (l)*(d)*(c) = (g).
      do 60 i=2,n
         i1 = i-1
         fac = g(i)
         do 50 j=1,i1
            fac = fac-g(j)*a(j,j)*a(i,j)
  50     continue
         g(i) = fac/a(i,i)
  60  continue
!c  second step : solve (l)'*(b) = (c)
      i = n
      do 80 j=2,n
         i1 = i
         i = i-1
         fac = g(i)
         do 70 k=i1,n
            fac = fac-g(k)*a(k,i)
  70     continue
         g(i) = fac
  80  continue
      return
      end
