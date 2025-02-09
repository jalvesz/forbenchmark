program benchmark_dot_coarray

#if defined(USE_COARRAY)

   use kinds
   use fordot
   use forbenchmark

   implicit none

   type(benchmark)       :: bench
   real(rk), allocatable :: u(:)
   real(rk), allocatable :: v(:)
   real(rk)              :: a
   integer(ik)           :: p
   integer               :: nl

   call bench%init(4, 'Fordot_product','benchmarks/dot/results/dot', 1000)

   do p = 5000_ik,100000_ik, 5000_ik

      !===============================================================================
      if (allocated(u)) deallocate(u)
      if (allocated(v)) deallocate(v)
      allocate(u(p))
      allocate(v(p))
      call random_number(u)
      call random_number(v)
      !===============================================================================


      !===============================================================================
      call bench%start_benchmark(1,'dot_product','a = dot_product(u,v)',[p])
      do nl = 1,bench%nloops
         u = u + real(nl,rk) ! to prevent compiler from optimizing (loop-invariant)
         v = v + real(nl,rk) ! to prevent compiler from optimizing (loop-invariant)
         a = dot_product(u,v)
      end do
      call bench%stop_benchmark(cmp_gflops)
      !===============================================================================


      !===============================================================================
      call bench%start_benchmark(2,'m1_co', "a = dot_product(u,v,'coarray','m1')",[p])
      do nl = 1,bench%nloops
         u = u + real(nl,rk) ! to prevent compiler from optimizing (loop-invariant)
         v = v + real(nl,rk) ! to prevent compiler from optimizing (loop-invariant)
         a = dot_product(u,v,'coarray','m1')
      end do
      call bench%stop_benchmark(cmp_gflops)
      !===============================================================================


      !===============================================================================
      call bench%start_benchmark(3,'m2_co', "a = dot_product(u,v,'coarray','m2')",[p])
      do nl = 1,bench%nloops
         u = u + real(nl,rk) ! to prevent compiler from optimizing (loop-invariant)
         v = v + real(nl,rk) ! to prevent compiler from optimizing (loop-invariant)
         a = dot_product(u,v,'coarray','m2')
      end do
      call bench%stop_benchmark(cmp_gflops)
      !===============================================================================


      !===============================================================================
      call bench%start_benchmark(4,'m3_co', "a = dot_product(u,v,'coarray','m3')",[p])
      do nl = 1,bench%nloops
         u = u + real(nl,rk) ! to prevent compiler from optimizing (loop-invariant)
         v = v + real(nl,rk) ! to prevent compiler from optimizing (loop-invariant)
         a = dot_product(u,v,'coarray','m3')
      end do
      call bench%stop_benchmark(cmp_gflops)
      !===============================================================================


   end do

   call bench%finalize()

contains

   !===============================================================================
   function cmp_gflops(argi,argr) result(gflops)
      integer(ik), dimension(:), intent(in), optional :: argi
      real(rk),    dimension(:), intent(in), optional :: argr
      real(rk)                                        :: gflops

      gflops = real(argi(1),kind=rk)*1.0e-9_rk
   end function cmp_gflops
   !===============================================================================

#endif

end program benchmark_dot_coarray
