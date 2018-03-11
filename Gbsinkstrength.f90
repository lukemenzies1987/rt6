function gbss( sss,rgr )
  use kindmod
  implicit none
  real(kind=dp) ::rkr,sss,fun,rgr,gbss


  rkr  = rgr*dsqrt(abs(sss))
  fun  = 3.*(rkr/dtanh(rkr)-1.)
  gbss = rgr**(-2._dp)*rkr**2._dp*fun/(rkr**2._dp-fun)
  return
end

