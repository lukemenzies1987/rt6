subroutine jac (neq, t, y, ml, mu, pd, nrowpd)
  integer neq, ml, mu, nrowpd
  double precision t
  real, dimension(neq) :: y
  real, dimension(nrowpd,neq)::pd
  !pd=0.0
  return
  end
