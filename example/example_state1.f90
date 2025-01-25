program example_state1
  use la_state
  implicit none
  type(linalg_state) :: err

  ! Create a state flag
  err = linalg_state(LINALG_VALUE_ERROR,'just an example with scalar ', &
  &                  'integer=',1,'real=',2.0,'complex=',(3.0,1.0),'and array ',[1,2,3],'inputs')

  ! Print flag
  print *, err%print()

  ! Check success
  print *, 'Check error: ',err%error()
  print *, 'Check flag : ',err == LINALG_VALUE_ERROR

end program example_state1
