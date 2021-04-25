program main
   use omp_lib

   print*,"I am the master, and I am alone"

!$OMP PARALLEL
   call say_hello(omp_get_thread_num(), omp_get_num_threads())
!$OMP END PARALLEL

end program main

subroutine say_hello(my_thread_num, total_threads)
   integer :: my_thread_num, total_threads
   print*,"Hello, I am thread ",my_thread_num,"among ",total_threads,"threads"
end subroutine say_hello
