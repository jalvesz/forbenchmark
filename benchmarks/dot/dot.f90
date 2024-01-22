program benchmark_dot

   use kinds
   use fordot
   use fast_math, only: fprod, fprod_kahan
   use forbenchmark

   implicit none

   type(benchmark)       :: bench
   real(rk), allocatable :: u(:)
   real(rk), allocatable :: v(:)
   real(rk)              :: a
   integer(ik)           :: p
   integer               :: nl
   integer               :: seed_size
   integer, allocatable  :: seed_array(:)

   call random_seed(size = seed_size)
   allocate(seed_array(seed_size))
   seed_array = 123456789

   call bench%init(7,'Benchmark dot_product','benchmarks/dot/results/dot', 1000)

   do p = 50000_ik,1000000_ik, 50000_ik

      !===============================================================================
      if (allocated(u)) deallocate(u)
      if (allocated(v)) deallocate(v)
      allocate(u(p))
      allocate(v(p))
      call random_seed(put = seed_array)
      !===============================================================================


      !===============================================================================
      call random_number(u)
      call random_number(v)
      call bench%start_benchmark(1,'dot_product','a = dot_product(u,v)',[p])
      a = 0._rk
      do nl = 1,bench%nloops
         a = a + dot_product(u,v)
      end do
      call bench%stop_benchmark(cmp_gflops)
      print *, a
      !===============================================================================


      !===============================================================================
      call random_number(u)
      call random_number(v)
      call bench%start_benchmark(2,'m1', "a = dot_product(u,v,'m1')",[p])
      a = 0._rk
      do nl = 1,bench%nloops
         a = a + dot_product(u,v,'m1')
      end do
      call bench%stop_benchmark(cmp_gflops)
      print *, a
      !===============================================================================


      !===============================================================================
      call random_number(u)
      call random_number(v)
      call bench%start_benchmark(3,'blas', "a = dot_product(u,v,'m2')",[p])
      a = 0._rk
      do nl = 1,bench%nloops
         a = a + dot_product(u,v,'m2')
      end do
      call bench%stop_benchmark(cmp_gflops)
      print *, a
      !===============================================================================


      !===============================================================================
      call random_number(u)
      call random_number(v)
      call bench%start_benchmark(4,'m3', "a = dot_product(u,v,'m3')",[p])
      a = 0._rk
      do nl = 1,bench%nloops
         a = a + dot_product(u,v,'m3')
      end do
      call bench%stop_benchmark(cmp_gflops)
      print *, a
      !===============================================================================


      !===============================================================================
      call random_number(u)
      call random_number(v)
      call bench%start_benchmark(5,'m4', "a = dot_product(u,v,'m4')",[p])
      a = 0._rk
      do nl = 1,bench%nloops
         a = a + dot_product(u,v,'m4')
      end do
      call bench%stop_benchmark(cmp_gflops)
      print *, a
      !===============================================================================

      !===============================================================================
      call random_number(u)
      call random_number(v)
      call bench%start_benchmark(6,'chunks', "a = fprod(u,v)",[p])
      a = 0._rk
      do nl = 1,bench%nloops
         a = a + fprod(u,v)
      end do
      call bench%stop_benchmark(cmp_gflops)
      print *, a
      !===============================================================================


      !===============================================================================
      call random_number(u)
      call random_number(v)
      call bench%start_benchmark(7,'kahan', "a = fprod_kahan(u,v)",[p])
      a = 0._rk
      do nl = 1,bench%nloops
         a = a + fprod_kahan(u,v)
      end do
      call bench%stop_benchmark(cmp_gflops)
      print *, a
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

end program benchmark_dot
