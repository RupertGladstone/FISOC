
program FISh_caller

  use FISh_main

  implicit none 

  integer :: tc

  call FISh_initialize()

  do tc = 1,tl
     call FISh_run(tc,h,u)
  end do

  call FISh_finalize()

end program FISh_caller

