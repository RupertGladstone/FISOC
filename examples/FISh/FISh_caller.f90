
program FISh_caller

  use FISh_MOD

  implicit none 

  call FISh_initialize()

  do tc = 1,tl
     call FISh_run()
  end do

  call FISh_finalize()

end program FISh_caller

