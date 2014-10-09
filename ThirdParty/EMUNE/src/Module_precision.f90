module precision
  ! define the precision to be used, and other constant values

  implicit none

  integer, parameter :: dp = selected_real_kind(15,307)

  ! user constants
  real(dp), parameter :: zero = 0.0_dp, One = 1.0_dp 
  real(dp), parameter :: pi = 3.14159265358979_dp
  real(dp), parameter, dimension(3) :: e1 = (/One, zero, zero/), &
                       e2 = (/zero, One, zero/), &
                       e3 = (/zero, zero, One/) ! unit vectors 
  

end module precision
