      module sia2_gpu

          implicit none

          integer, parameter :: gpu_tpb = 192  ! threads per block
          ! 268435456=256 402653184=384 469762048=448 805306368=768 1073741824=1024 1879048192=1.75Gb
          integer, parameter :: gpu_mem = 1879048192 !  total memory requestd on GPU

          integer, parameter :: size_of_int = 4
          integer, parameter :: size_of_float = 4

          real, parameter :: log9m(9) = (/0.102e0, 0.272e0, 0.427e0, &
                 0.532e0, 0.721e0, 0.952e0, 1.31e0, 1.74e0, 3.31e0/)

          real, parameter :: log6m(6) = (/0.145e0, 0.385e0, 0.585e0, &
                 0.860e0, 1.345e0, 2.70e0/)

          real, save ::            &
                 h_inv_c_chl,      &  ! host inverse caron/chla ratio
                 h_inv_iced,       &  ! host inverse pure ice densiry
                 h_a_ice_ir,       &  ! host near-ir abosroption
                 h_inv_ad_denom,   &  ! host inverse ad_denom
                 h_a_factor           ! host snow extra absorption factor

          integer, save, pointer :: &
              c_i32_in(:), &   !
              mi_gpu_index(:)

          real, save, pointer :: &
              c_f32_in(:),  & !
              c_f32_scr(:), &   !
              c_f32_out(:)   !

          integer, save :: &
              g_i32_in,  & ! {nilyr(1), nslyr(1), srftyp(1)}
              g_f32_in,  & ! { am_sum(1), snow_d/bv(ice_z), smalg(ice_z),poc(ice_z),coszen(1),g(ice_z),r2st(ice_z),dth(ice_z),ed0_dir(wavl+1),ed0_dif(wavl_1)}
              g_f32_scr, & ! { 11*(layers+1) scatch memory for k2st, edd params, aops }
              g_f32_out    ! { PUR (ice_z),Ed_W(ice_z) }

          integer, save :: &
              i32_in_p_th,     &
              f32_in_p_th,     &
              f32_scr_p_th,    &
              f32_out_p_th,    &
              i32_in_size,     &
              f32_in_size,     &
              f32_scr_size,    &
              f32_out_size,    &
              gpu_max_threads, &
              gpu_max_blocks

          contains


          function find_logm(i,lda_n)

              integer, intent(IN) :: i,lda_n
              real :: find_logm,logm

              if (lda_n .eq. 1) then
                  logm = 1
              elseif (lda_n .eq. 9) then
                  logm = log9m(i)
              elseif(lda_n .eq. 6) then
                  logm = log6m(i)
              endif

              find_logm = logm

          end function find_logm


          subroutine gpu_allocate(layers)

              integer :: layers,thread_size

              i32_in_p_th = 3                            ! thread total in ints
              f32_in_p_th = 1 + 6*(layers) + 2*(32)      ! thread total in floats
              f32_scr_p_th = 11*(layers+1)               ! thread total scratch floats
              f32_out_p_th = 2*(layers)                  ! thread total out floats

              thread_size = i32_in_p_th*size_of_int + &
                (f32_in_p_th + f32_scr_p_th + f32_out_p_th)*size_of_float

              gpu_max_blocks = gpu_mem/(gpu_tpb * thread_size)

              gpu_max_threads = gpu_max_blocks * gpu_tpb

              i32_in_size = gpu_max_threads*i32_in_p_th
              print *,'Allocating CUDA array: c_i32_in ', &
                  (size_of_int*i32_in_size)/1024**2,'MB'
              allocate(c_i32_in(i32_in_size))

              f32_in_size = gpu_max_threads*f32_in_p_th
              print *,'Allocating CUDA array: c_f32_in ', &
                 (size_of_int*f32_in_size)/1024**2,'MB'
              allocate(c_f32_in(f32_in_size))

              f32_scr_size = gpu_max_threads*f32_scr_p_th
              print *,'Allocating CUDA array: c_f32_scr ', &
                 (size_of_int*f32_scr_size)/1024**2,'MB'
              allocate(c_f32_scr(f32_scr_size))

              f32_out_size = gpu_max_threads*f32_out_p_th
              print *,'Allocating CUDA array: c_f32_out ', &
                  (size_of_int*f32_out_size)/1024**2,'MB'
              allocate(c_f32_out(f32_out_size))

              print *,'Allocating CUDA index array: ',gpu_max_threads,'= max gpu threads'
              allocate(mi_gpu_index(gpu_max_threads))

          end subroutine gpu_allocate


      end module sia2_gpu
