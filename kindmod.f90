module kindmod
  Use nag_library, Only: nag_wp
  implicit none
  integer, parameter :: dp=selected_real_kind(15,300)!dp=selected_real_kind(8)
  integer, parameter :: ip=selected_int_kind(16)
  integer, parameter :: rp=selected_real_kind(8) 
 
end module kindmod
