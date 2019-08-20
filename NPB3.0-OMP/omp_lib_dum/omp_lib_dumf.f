      integer function omp_get_num_threads()
      omp_get_num_threads = 1
      return
      end

      integer function omp_get_thread_num()
      omp_get_thread_num = 0
      return
      end

      integer function omp_get_max_threads()
      omp_get_max_threads = 1
      return
      end
