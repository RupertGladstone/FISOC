
program FISh_caller

  use FISh_MOD

  implicit none 

  integer :: tc_count

  call FISh_initialize()

print*,tl," tl"
  do tc_count = 1,tl
     call FISh_run()
  end do

  call FISh_finalize()

end program FISh_caller

