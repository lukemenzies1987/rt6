Subroutine monitr(neq,ldysav,t,hlast,hnext,y,ydot,ysav,r,acor,imon,inln, &
  hmin,hmax,nqu)

!       .. Use Statements ..
  Use nag_library, Only: d02zaf,nag_wp
!       .. Scalar Arguments ..
  Real (Kind=nag_wp), Intent (In)      :: hlast, t
  Real (Kind=nag_wp), Intent (Inout)   :: hmax, hmin, hnext
  Integer, Intent (Inout)              :: imon
  Integer, Intent (Out)                :: inln
  Integer, Intent (In)                 :: ldysav, neq, nqu
!       .. Array Arguments ..
  Real (Kind=nag_wp), Intent (In)      :: acor(neq,2), r(neq),           &
                                          ydot(neq), ysav(ldysav,*)
  Real (Kind=nag_wp), Intent (Inout)   :: y(neq)
!       .. Local Scalars ..
  Real (Kind=nag_wp)                   :: errloc
  Integer                              :: i, ifail,nout
!       .. Executable Statements ..
  inln = 3
  nout=78
  If (imon==-1) Then

    ifail = -1

    do i=1,neq
      errloc = d02zaf(neq,acor(i,2),acor(i,1),ifail)
      Write (nout,99999) t, y(i),i, errloc
    end do
  imon=2
  End If


  Return 
99999   Format (2(EN20.10,X),I3,' ** WARNING scaled local error = ', &
          (EN20.10,X))
99998   Format (1X,F10.6,3(F13.7,2X))
End Subroutine monitr
