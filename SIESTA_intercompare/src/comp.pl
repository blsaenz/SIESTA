

# init main variables
our $makefile,$machine,$comp,$opt,$omp,$opt1,$opt2,$libpath,$incl_path,$optlink,$endian_switch,$openmp_switch,$extra_shit;

# grab line arguments
$machine = $ARGV[0];
$comp = $ARGV[1];
$opt = $ARGV[2];
$omp = $ARGV[3];
$gpu = $ARGV[4];

if ($omp eq "gpu") {
    $gpu = "gpu";
}

# create makefile
&load_make_parts();
&create_makefile();
open(MAKEFILE, ">Makefile") || warn "Can't open makefile!!\n";
print MAKEFILE $makefile;
close(MAKEFILE);

# build
system("make clean");
system("make pre");
system("make");



sub load_make_parts () {

    $openmp_switch = "";
    $endian_switch = "";
    $gpu_make_objects = "";
    $gpu_make_code = "";
    $gpu_cpp_code = "";
    $gpu_cpp_switch = "";
    $gpulink = "";
    $comp_ifort = "";
    $comp_gfort = "";

		if($comp eq "ifort") { $comp_ifort="-DIFORT"; }
		if($comp eq "gfortan") { $comp_gfort="-DGFORT";	}



# ----------------------------------------------------------------------
# MAC
# ----------------------------------------------------------------------
#LIB = -L../../support/netcdf-3.6.3-gfortran_64_mac_10.11/lib -lnetcdf
#INC = -I../../support/netcdf-3.6.3-gfortran_64_mac_10.11/include
#LAPACK = -L/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/Current/ -lBLAS -lLAPACK

    if ($machine eq "mac") {
        if ($comp eq "gfortran") {
            $lib_path = "~/Dropbox/Projects/KPP_ECO_ICE/support/netcdf-3.6.3-gfortran_64_mac_10.11/lib";
            $incl_path = "~/Dropbox/Projects/KPP_ECO_ICE/support/netcdf-3.6.3-gfortran_64_mac_10.11/include";
            $optlink =  "-L/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/Current -lBLAS -lLAPACK \$(lib_path)/libnetcdf.a ";
            if ($opt == 0) {
                $opt1 = "-g -O0 -fbounds-check -m64";
                $opt2 = $opt1
            } elsif ($opt == 1) {
                $opt1 = "-g -O0  -m64";
                $opt2 = $opt1
            } elsif ($opt == 2) {
                $opt1 = "-g -m64 -O2";
                $opt2 = "-g -m64 -O1"
            } elsif ($opt == 3) {
                $opt1 = "-g -O3 -m64 -ftree-vectorize";
                $opt2 = "-g -m64 -O1"
            }
            if ($omp eq "omp") {
                $opt1 = $opt1 . " -fopenmp";
            }
        }
    }

# ----------------------------------------------------------------------
# PETREL-old
# ----------------------------------------------------------------------


    if ($machine eq "ptold") {
        if ($comp eq "ifort") {

            $lib_path = "/Users/blsaenz/opt/netcdf-3.6.3-ifort64/lib";
            $incl_path = "/Users/blsaenz/opt/netcdf-3.6.3-ifort64/include";
            $optlink =  "-L\$(MKL_LIB) -I\$(MKL_INC) -lmkl_lapack95_lp64 \$(MKL_LIB)/libmkl_intel_lp64.a \$(MKL_LIB)/libmkl_sequential.a \$(MKL_LIB)/libmkl_core.a \$(MKL_LIB)/libmkl_sequential.a \$(MKL_LIB)/libmkl_core.a \$(MKL_LIB)/libmkl_sequential.a \$(MKL_LIB)/libmkl_core.a -lpthread \$(lib_path)/libnetcdf.a ";
#            $optlink =  "-L/Developer/SDKs/MacOSX10.5.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/ -lBLAS -lLAPACK \$(lib_path)/libnetcdf.a ";
            if ($opt == 0) {
                $opt1 = "-g -m64  -use_asm -O0 -fp-model source -assume minus0 -extend_source -CB -traceback -fpe-all=0 -debug all -heap-arrays -fstack-protector";
            } elsif ($opt == 1) {
                $opt1 = "-g -m64  -use_asm -O0 -fp-model source -assume minus0 -extend_source -traceback -fpe-all=0 -debug all -heap-arrays -fstack-protector";
            } elsif ($opt == 2) {
                $opt1 = "-g -m64  -use_asm -O1 -fp-model source -assume minus0 -extend_source -traceback -fpe-all=0 -debug all -heap-arrays -fstack-protector";
            } elsif ($opt == 3) {
                $opt1 = "-g -m64  -use_asm -vec-report0 -O2 -xS -ip -fp-model source -assume minus0 -extend_source -align -heap-arrays -fpe-all=0 -traceback -fstack-protector";
            } elsif ($opt == 4) {
                $opt1 = "-g -m64  -use_asm -vec-report0 -O3 -xS -ip -fp-model source -assume minus0 -extend_source -align -heap-arrays -fpe-all=0 -traceback -fstack-protector";
            } elsif ($opt == 5) {
                $opt1 = "-m64  -use_asm -vec-report0 -O3 -xS -ip -fp-model source -assume minus0 -extend_source -align -heap-arrays -fpe-all=0 -fstack-protector";
            } elsif ($opt == 6) {
                $opt1 = "-g -m64 -use_asm -vec-report0 -O3 -xS -ip -no-prec-div -extend_source -align -heap-arrays -fpe-all=0 -traceback -fstack-protector";
            } elsif ($opt == 7) {
                $opt1 = "-m64 -use_asm -vec-report0 -O3 -xS -ip -no-prec-div -extend_source -align -heap-arrays -fpe-all=0 -fstack-protector";
            }
            if ($omp eq "omp") {
                $opt1 = $opt1 . " -openmp";
            }
            $openmp_switch = "-Domp_on";
            $extra_shit = <<EOT
MKL_INC = /opt/intel/composerxe/mkl/include
MKL_LIB = /opt/intel/composerxe/mkl/lib
EOT
;
            if ($gpu eq "gpu") {
                $gpu_cpp_switch = "-DGPU";
                $gpu_make_objects = "sia2_gpu.o sia2_edd_gpu.o sia2_edd_gpu_calc.o ";
                $gpu_make_code = <<EOT
sia2_edd_gpu_calc.o:
	$comp -c \$(opt1) sia2_edd_gpu_calc.f90
sia2_gpu.o:
	$comp -c \$(opt1) sia2_gpu.f90
sia2_edd_gpu.o:
	/usr/local/cuda/bin/nvcc -g -c -Xptxas="-v" -m64 -O3 -arch compute_11 -code compute_11 sia2_edd_gpu.cu -I\"/Developer/GPU Computing/C/common/inc"
EOT
;
# put -deviceemu in nvcc compile line for emulation mode
# nvcc -c -m64 -Xptxas="-v" sia2_edd_gpu.cu -I"/Developer/GPU Computing/C/common/inc"
                $gpu_cpp_code = "";
                $gpulink = "-L/usr/local/cuda/lib -lcuda -lcudart -L\"/Developer/GPU Computing/C/lib\" -lcutil_x86_64 -L/usr/lib -lstdc++";
#                $gpulink = "-L/usr/local/cuda/lib -lcudart -L/usr/lib/ -lstdc++";
             }

        } elsif ($comp eq "gfortran") {
            $lib_path = "/Users/blsaenz/opt/lib";
            $incl_path = "/Users/blsaenz/opt/include";
            $optlink =  "-L/Developer/SDKs/MacOSX10.7.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/ -lBLAS -lLAPACK \$(lib_path)/libnetcdf.a ";
            if ($opt == 0) {
                $opt1 = "-g -O0 -fbounds-check -fbacktrace -m64";
            } elsif ($opt == 1) {
                $opt1 = "-g -O0 -fbacktrace -m64";
            } elsif ($opt == 2) {
                $opt1 = "-g -m64 -O2 -fbacktrace";
            } elsif ($opt == 3) {
                $opt1 = "-g -O3 -fbounds-check -mtune=native -march=native -m64 -mssse3 -ftree-vectorize -fbacktrace";
            }
            if ($omp eq "omp") {
                $opt1 = $opt1 . " -fopenmp";
            }
        }
    }

# ----------------------------------------------------------------------
# ALBEDO
# ----------------------------------------------------------------------

    if ($machine eq "albedo") {
      if ($comp eq "gfortran") {
            $lib_path = "/opt/netcdf-3.6.3-gfortran64/lib";
            $incl_path = "/opt/netcdf-3.6.3-gfortran64/include";
            $optlink =  "-L/usr/lib64 -lBLAS -lLAPACK \$(lib_path)/libnetcdf.a ";
            if ($opt == 0) {
                $opt1 = "-g -O0 -fbounds-check -fbacktrace -m64";
            } elsif ($opt == 1) {
                $opt1 = "-g -O0 -fbacktrace -m64";
            } elsif ($opt == 2) {
                $opt1 = "-g -m64 -O2 -fbacktrace";
            } elsif ($opt == 3) {
                $opt1 = "-g -O3 -fbounds-check -mtune=native -march=native -m64 -mssse3 -ftree-vectorize -fbacktrace";
            }
            if ($omp eq "omp") {
                $opt1 = $opt1 . " -fopenmp";
            }
        }
    }

# ----------------------------------------------------------------------
# ICY
# ----------------------------------------------------------------------

    if ($machine eq "i") {
        if ($comp eq "ifort") {
            $lib_path = "/opt/netcdf-3.6.3-ifort64/lib";
            $incl_path = "/opt/netcdf-3.6.3-ifort64/include";
            $optlink =  "-L\$(MKL_LIB) -I\$(MKL_INC) -lmkl_lapack95_lp64 \$(MKL_LIB)/libmkl_intel_lp64.a \$(MKL_LIB)/libmkl_sequential.a \$(MKL_LIB)/libmkl_core.a \$(MKL_LIB)/libmkl_sequential.a \$(MKL_LIB)/libmkl_core.a \$(MKL_LIB)/libmkl_sequential.a \$(MKL_LIB)/libmkl_core.a -lpthread \$(lib_path)/libnetcdf.a ";
            if ($opt == 0) {
                $opt1 = "-g -m64 -O0 -fp-model source -extend_source -CB -traceback -fpe-all=0 -debug all -heap-arrays -fstack-protector";
            } elsif ($opt == 1) {
                $opt1 = "-g -m64 -O0 -fp-model source -extend_source -traceback -fpe-all=0 -debug all -heap-arrays -fstack-protector";
            } elsif ($opt == 2) {
                $opt1 = "-g -m64 -O1 -fp-model source -extend_source -traceback -fpe-all=0 -debug all -heap-arrays -fstack-protector";
            } elsif ($opt == 3) {
                $opt1 = "-g -m64  -vec-report0 -O2 -xSSE4.2 -ip -fp-model source -align -assume minus0 -xHost -extend_source -align -heap-arrays -fpe-all=0 -traceback -fstack-protector";
            } elsif ($opt == 4) {
                $opt1 = "-g -m64  -vec-report0 -O3 -xSSE4.2 -ip -fp-model source -assume minus0 -xHost -extend_source -align -heap-arrays -fpe-all=0 -traceback -fstack-protector";
            } elsif ($opt == 5) {
                $opt1 = "-g -m64  -vec-report0 -O3 -xSSE4.2 -ip -no-prec-div -xHost -extend_source -align -heap-arrays -fpe-all=0 -traceback -fstack-protector";
            }
            if ($omp eq "omp") {
                $opt1 = $opt1 . " -openmp";
                $openmp_switch = "-Domp_on";
            }
            $extra_shit = <<EOT
MKL_INC = /opt/intel/Compiler/11.1/072/mkl/include
MKL_LIB = /opt/intel/Compiler/11.1/072/mkl/lib/em64t
EOT
;
            if ($gpu eq "gpu") {
                $gpu_cpp_switch = "-DGPU";
                $gpu_make_objects = "sia2_gpu.o sia2_edd_gpu.o sia2_edd_gpu_calc.o ";
                $gpu_make_code = <<EOT
sia2_edd_gpu_calc.o:
	$comp -c \$(opt1) sia2_edd_gpu_calc.f90
sia2_gpu.o:
	$comp -c \$(opt1) sia2_gpu.f90
sia2_edd_gpu.o:
	/usr/local/cuda/bin/nvcc -g -c -Xptxas="-v" -m64 -O3 -arch compute_20 -code compute_20 sia2_edd_gpu.cu -I\"/home/blsaenz/NVIDIA_GPU_Computing_SDK/C/common/inc"
EOT
;
# put -deviceemu in nvcc compile line for emulation mode
# nvcc -c -m64 -Xptxas="-v" sia2_edd_gpu.cu -I"/Developer/GPU Computing/C/common/inc"
                $gpu_cpp_code = "";
                $gpulink = "-L/usr/local/cuda/lib64 -lcuda -lcudart -L\"/home/blsaenz/NVIDIA_GPU_Computing_SDK/C/lib\" -lcutil_x86_64 -L/usr/lib64 -lstdc++";
#                $gpulink = "-L/usr/local/cuda/lib -lcudart -L/usr/lib/ -lstdc++";
             }
        } elsif ($comp eq "gfortran") {
            $lib_path = "/opt/netcdf-3.6.3-gfort64/lib";
            $incl_path = "/opt/netcdf-3.6.3-gfort64/include";
            $optlink =  "-L/opt/lapack-lite-3.1.1-gfort64 -lBLAS -lLAPACK \$(lib_path)/libnetcdf.a ";
            if ($opt == 0) {
                $opt1 = "-g -O0 -fbounds-check -fbacktrace -m64";
            } elsif ($opt == 1) {
                $opt1 = "-g -O0 -fbacktrace -m64";
            } elsif ($opt == 2) {
                $opt1 = "-g -m64 -O2 -fbacktrace";
            } elsif ($opt == 3) {
                $opt1 = "-g -O3 -fbounds-check -mtune=native -march=native -m64 -mssse3 -ftree-vectorize -fbacktrace";
            }
            if ($omp eq "omp") {
                $opt1 = $opt1 . " -fopenmp";
            }
        } elsif ($comp eq "f95") {
            $lib_path = "/opt/netcdf-3.6.3-f95_64/lib";
            $incl_path = "/opt/netcdf-3.6.3-f95_64/include";
            $optlink =  "\$(lib_path)/libnetcdf.a ";
            if ($opt == 0) {
                $opt1 = "-g -C -dalign -m64 -xarch=native -xchip=native -xcache=native -e -xcheck=stkovf -xlic_lib=sunperf ";
            } elsif ($opt == 1) {
                $opt1 = "-g -dalign -m64 -xarch=native -xchip=native -xcache=native -e -xcheck=stkovf -xlic_lib=sunperf ";
            } elsif ($opt == 2) {
                $opt1 = "-g -dalign -m64 -xarch=native -xchip=native -xcache=native -e -xcheck=stkovf -xlic_lib=sunperf ";
            } elsif ($opt == 3) {
                $opt1 = "-g -dalign -fast -m64 -xarch=native -xchip=native -xcache=native -e -xlic_lib=sunperf ";
            }
            if ($omp eq "omp") {
                if ($opt < 3) {
                    $opt1 = $opt1 . " -xopenmp=noopt";
                } else {
                    $opt1 = $opt1 . " -xopenmp";
                }
            }
        }
    }

# ----------------------------------------------------------------------
# SNOW
# ----------------------------------------------------------------------

    if ($machine eq "s") {
        if ($comp eq "ifort") {
            $lib_path = "/opt/netcdf-4.2.1-ifort64-system-zlib/lib";
            $incl_path = "/opt/netcdf-4.2.1-ifort64-system-zlib/include";
            $optlink =  "-L\$(MKL_LIB) -I\$(MKL_INC) -lmkl_lapack95_lp64 \$(MKL_LIB)/libmkl_intel_lp64.a \$(MKL_LIB)/libmkl_sequential.a \$(MKL_LIB)/libmkl_core.a \$(MKL_LIB)/libmkl_sequential.a \$(MKL_LIB)/libmkl_core.a \$(MKL_LIB)/libmkl_sequential.a \$(MKL_LIB)/libmkl_core.a -lpthread -L\$(lib_path) -lnetcdff ";
            if ($opt == 0) {
                $opt1 = "-g -m64 -O0 -fp-model source -extend_source -CB -traceback -fpe-all=0 -debug all -heap-arrays -fstack-protector";
								$opt2 = $opt1;
            } elsif ($opt == 1) {
                $opt1 = "-g -m64 -O0 -fp-model source -extend_source -traceback -fpe-all=0 -debug all -heap-arrays -fstack-protector";
								$opt2 = $opt1;
            } elsif ($opt == 2) {
                $opt1 = "-g -m64 -O1 -fp-model source -extend_source -traceback -fpe-all=0 -debug all -heap-arrays -fstack-protector";
								$opt2 = $opt1;
            } elsif ($opt == 3) {
                $opt1 = "-g -m64  -vec-report0 -O2 -xSSE4.2 -ip -fp-model source -align -assume minus0 -xHost -extend_source -align -heap-arrays -traceback -fpe-all=0 -debug all -fstack-protector";
                $opt2 = "-g -m64 -O1 -fp-model source -extend_source -traceback -fpe-all=0 -debug all -heap-arrays -fstack-protector";
            } elsif ($opt == 4) {
                $opt1 = "-m64  -vec-report0 -O2 -xSSE4.2 -ip -fp-model source -align -assume minus0 -xHost -extend_source -align -heap-arrays -fpe-all=0 -fstack-protector";
                $opt2 = "-m64 -O1 -fp-model source -extend_source -fpe-all=0 -heap-arrays -fstack-protector";
            } elsif ($opt == 5) {
                $opt1 = "-g -m64  -vec-report0 -O3 -xSSE4.2 -ip -fp-model source -assume minus0 -xHost -extend_source -align -heap-arrays -fpe-all=0 -traceback -fstack-protector";
                $opt2 = "-g -m64  -vec-report0 -O1 -xSSE4.2 -ip -fp-model source -assume minus0 -xHost -extend_source -align -heap-arrays -fpe-all=0 -traceback -fstack-protector";
            } elsif ($opt == 6) {
                $opt1 = "-g -m64  -vec-report0 -O3 -xSSE4.2 -ip -no-prec-div -xHost -extend_source -align -heap-arrays -fpe-all=0 -traceback -fstack-protector";
								$opt2 = $opt1;
            }
            if ($omp eq "omp") {
                $opt1 = $opt1 . " -openmp";
                $opt2 = $opt2 . " -openmp";
               $openmp_switch = "-Domp_on";
            }
            $extra_shit = <<EOT
MKL_INC = /opt/intel/Compiler/11.1/073/mkl/include
MKL_LIB = /opt/intel/Compiler/11.1/073/mkl/lib/em64t
EOT
;
            if ($gpu eq "gpu") {
                $gpu_cpp_switch = "-DGPU";
                $gpu_make_objects = "sia2_gpu.o sia2_edd_gpu.o sia2_edd_gpu_calc.o ";
                $gpu_make_code = <<EOT
sia2_edd_gpu_calc.o:
	$comp -c \$(opt1) sia2_edd_gpu_calc.f90
sia2_gpu.o:
	$comp -c \$(opt1) sia2_gpu.f90
sia2_edd_gpu.o:
	/usr/local/cuda/bin/nvcc -g -c -Xptxas="-v" -m64 -O3 -arch compute_20 -code compute_20 sia2_edd_gpu.cu -I\"/home/blsaenz/CUDA_SDK/C/common/inc"
EOT
;
# put -deviceemu in nvcc compile line for emulation mode
# nvcc -c -m64 -Xptxas="-v" sia2_edd_gpu.cu -I"/Developer/GPU Computing/C/common/inc"
                $gpu_cpp_code = "";
                $gpulink = "-L/usr/local/cuda/lib64 -lcuda -lcudart -L\"/home/blsaenz/CUDA_SDK/C/lib\" -lcutil_x86_64 -L/usr/lib64 -lstdc++";
#                $gpulink = "-L/usr/local/cuda/lib -lcudart -L/usr/lib/ -lstdc++";
             }
        } elsif ($comp eq "gfortran") {
            $lib_path = "/opt/netcdf-3.6.3-gf64/lib";
            $incl_path = "/opt/netcdf-3.6.3-gf64/include";
            $optlink =  "-L/opt/lapack-lite-3.1.1-gf64 -lBLAS -lLAPACK \$(lib_path)/libnetcdf.a ";
            if ($opt == 0) {
                $opt1 = "-g -O0 -fbounds-check -m64";
            } elsif ($opt == 1) {
                $opt1 = "-g -O0  -m64";
            } elsif ($opt == 2) {
                $opt1 = "-g -m64 -O2 ";
            } elsif ($opt == 3) {
                $opt1 = "-g -O3 -fbounds-check -mtune=native -march=native -m64 -mssse3 -ftree-vectorize ";
            }
            if ($omp eq "omp") {
                $opt1 = $opt1 . " -fopenmp";
            }
        } elsif ($comp eq "f95") {
            $lib_path = "/opt/netcdf-3.6.3-f95_64/lib";
            $incl_path = "/opt/netcdf-3.6.3-f95_64/include";
            $optlink =  "\$(lib_path)/libnetcdf.a ";
            if ($opt == 0) {
                $opt1 = "-g -C -dalign -m64 -xarch=native -xchip=native -xcache=native -e -xcheck=stkovf -xlic_lib=sunperf ";
            } elsif ($opt == 1) {
                $opt1 = "-g -dalign -m64 -xarch=native -xchip=native -xcache=native -e -xcheck=stkovf -xlic_lib=sunperf ";
            } elsif ($opt == 2) {
                $opt1 = "-g -dalign -m64 -xarch=native -xchip=native -xcache=native -e -xcheck=stkovf -xlic_lib=sunperf ";
            } elsif ($opt == 3) {
                $opt1 = "-g -dalign -fast -m64 -xarch=native -xchip=native -xcache=native -e -xlic_lib=sunperf ";
            }
            if ($omp eq "omp") {
                if ($opt < 3) {
                    $opt1 = $opt1 . " -xopenmp=noopt";
                } else {
                    $opt1 = $opt1 . " -xopenmp";
                }
            }
        }
    }


# ----------------------------------------------------------------------
# RAFT
# ----------------------------------------------------------------------

    if ($machine eq "r") {
        if ($comp eq "ifort") {
            $lib_path = "/opt/netcdf-4.2.1-ifort64-system-zlib/lib";
            $incl_path = "/opt/netcdf-4.2.1-ifort64-system-zlib/include";
            $optlink =  "-L\$(MKL_LIB) -I\$(MKL_INC) -lmkl_lapack95_lp64 \$(MKL_LIB)/libmkl_intel_lp64.a \$(MKL_LIB)/libmkl_sequential.a \$(MKL_LIB)/libmkl_core.a \$(MKL_LIB)/libmkl_sequential.a \$(MKL_LIB)/libmkl_core.a \$(MKL_LIB)/libmkl_sequential.a \$(MKL_LIB)/libmkl_core.a -lpthread -L\$(lib_path) -lnetcdff ";
            if ($opt == 0) {
                $opt1 = "-g -m64 -O0 -fp-model source -extend_source -CB -traceback -fpe-all=0 -debug all -heap-arrays -fstack-protector";
								$opt2 = $opt1;
            } elsif ($opt == 1) {
                $opt1 = "-g -m64 -O0 -fp-model source -extend_source -traceback -fpe-all=0 -debug all -heap-arrays -fstack-protector";
								$opt2 = $opt1;
            } elsif ($opt == 2) {
                $opt1 = "-g -m64 -O1 -fp-model source -extend_source -traceback -fpe-all=0 -debug all -heap-arrays -fstack-protector";
								$opt2 = $opt1;
            } elsif ($opt == 3) {
                $opt1 = "-g -m64  -vec-report0 -O2 -xSSE4.2 -ip -fp-model source -align -assume minus0 -xHost -extend_source -align -heap-arrays -traceback -fpe-all=0 -debug all -fstack-protector";
                $opt2 = "-g -m64 -O1 -fp-model source -extend_source -traceback -fpe-all=0 -debug all -heap-arrays -fstack-protector";
            } elsif ($opt == 4) {
                $opt1 = "-m64  -vec-report0 -O2 -xSSE4.2 -ip -fp-model source -align -assume minus0 -xHost -extend_source -align -heap-arrays -fpe-all=0 -fstack-protector";
                $opt2 = "-m64 -O1 -fp-model source -extend_source -fpe-all=0 -heap-arrays -fstack-protector";
            } elsif ($opt == 5) {
                $opt1 = "-g -m64  -vec-report0 -O3 -xSSE4.2 -ip -fp-model source -assume minus0 -xHost -extend_source -align -heap-arrays -fpe-all=0 -traceback -fstack-protector";
                $opt2 = "-g -m64  -vec-report0 -O1 -xSSE4.2 -ip -fp-model source -assume minus0 -xHost -extend_source -align -heap-arrays -fpe-all=0 -traceback -fstack-protector";
            } elsif ($opt == 6) {
                $opt1 = "-g -m64  -vec-report0 -O3 -xSSE4.2 -ip -no-prec-div -xHost -extend_source -align -heap-arrays -fpe-all=0 -traceback -fstack-protector";
								$opt2 = $opt1;
            }
            if ($omp eq "omp") {
                $opt1 = $opt1 . " -openmp";
                $opt2 = $opt2 . " -openmp";
               $openmp_switch = "-Domp_on";
            }
            $extra_shit = <<EOT
MKL_INC = /opt/intel/Compiler/11.1/073/mkl/include
MKL_LIB = /opt/intel/Compiler/11.1/073/mkl/lib/em64t
EOT
;
            if ($gpu eq "gpu") {
                $gpu_cpp_switch = "-DGPU";
                $gpu_make_objects = "sia2_gpu.o sia2_edd_gpu.o sia2_edd_gpu_calc.o ";
                $gpu_make_code = <<EOT
sia2_edd_gpu_calc.o:
	$comp -c \$(opt1) sia2_edd_gpu_calc.f90
sia2_gpu.o:
	$comp -c \$(opt1) sia2_gpu.f90
sia2_edd_gpu.o:
	/usr/local/cuda/bin/nvcc -g -c -Xptxas="-v" -m64 -O3 -arch compute_20 -code compute_20 sia2_edd_gpu.cu -I\"/home/blsaenz/CUDA_SDK/C/common/inc"
EOT
;
# put -deviceemu in nvcc compile line for emulation mode
# nvcc -c -m64 -Xptxas="-v" sia2_edd_gpu.cu -I"/Developer/GPU Computing/C/common/inc"
                $gpu_cpp_code = "";
                $gpulink = "-L/usr/local/cuda/lib64 -lcuda -lcudart -L\"/home/blsaenz/CUDA_SDK/C/lib\" -lcutil_x86_64 -L/usr/lib64 -lstdc++";
#                $gpulink = "-L/usr/local/cuda/lib -lcudart -L/usr/lib/ -lstdc++";
             }
        } elsif ($comp eq "gfortran") {
            $lib_path = "/opt/netcdf-3.6.3-gf64/lib";
            $incl_path = "/opt/netcdf-3.6.3-gf64/include";
            $optlink =  "-L/opt/lapack-lite-3.1.1-gf64 -lBLAS -lLAPACK \$(lib_path)/libnetcdf.a ";
            if ($opt == 0) {
                $opt1 = "-g -O0 -fbounds-check -m64";
            } elsif ($opt == 1) {
                $opt1 = "-g -O0  -m64";
            } elsif ($opt == 2) {
                $opt1 = "-g -m64 -O2 ";
            } elsif ($opt == 3) {
                $opt1 = "-g -O3 -fbounds-check -mtune=native -march=native -m64 -mssse3 -ftree-vectorize ";
            }
            if ($omp eq "omp") {
                $opt1 = $opt1 . " -fopenmp";
            }
        } elsif ($comp eq "f95") {
            $lib_path = "/opt/netcdf-3.6.3-f95_64/lib";
            $incl_path = "/opt/netcdf-3.6.3-f95_64/include";
            $optlink =  "\$(lib_path)/libnetcdf.a ";
            if ($opt == 0) {
                $opt1 = "-g -C -dalign -m64 -xarch=native -xchip=native -xcache=native -e -xcheck=stkovf -xlic_lib=sunperf ";
            } elsif ($opt == 1) {
                $opt1 = "-g -dalign -m64 -xarch=native -xchip=native -xcache=native -e -xcheck=stkovf -xlic_lib=sunperf ";
            } elsif ($opt == 2) {
                $opt1 = "-g -dalign -m64 -xarch=native -xchip=native -xcache=native -e -xcheck=stkovf -xlic_lib=sunperf ";
            } elsif ($opt == 3) {
                $opt1 = "-g -dalign -fast -m64 -xarch=native -xchip=native -xcache=native -e -xlic_lib=sunperf ";
            }
            if ($omp eq "omp") {
                if ($opt < 3) {
                    $opt1 = $opt1 . " -xopenmp=noopt";
                } else {
                    $opt1 = $opt1 . " -xopenmp";
                }
            }
        }
    }


# ----------------------------------------------------------------------
# CEES-SPARC
# ----------------------------------------------------------------------

    if ($machine eq "c") {
        if  ($comp eq "f95") {
            $lib_path = "/home/blsaenz/netcdf-3.6.3/lib";
            $incl_path = "/home/blsaenz/netcdf-3.6.3/include";
            $optlink =  "\$(lib_path)/libnetcdf.a ";
            if ($opt == 0) {
                $opt1 = "-g -C -dalign -m64 -xarch=sparcvis2 -xchip=native -xcache=native -e -xcheck=stkovf -xlic_lib=sunperf ";
            } elsif ($opt == 1) {
                $opt1 = "-g -C -dalign -m64 -xarch=sparcvis2 -xchip=native -xcache=native -e -xcheck=stkovf -xlic_lib=sunperf ";
            } elsif ($opt == 2) {
                $opt1 = "-g -C -dalign -m64 -xarch=sparcvis2 -xchip=native -xcache=native -e -xcheck=stkovf -xlic_lib=sunperf ";
            } elsif ($opt == 3) {
                $opt1 = "-g -dalign -fast -m64 -xarch=sparcvis2 -xchip=native -xcache=native -e -xlic_lib=sunperf ";
            }
            if ($omp eq "omp") {
                if ($opt < 3) {
                    $opt1 = $opt1 . " -xopenmp=noopt";
                } else {
                    $opt1 = $opt1 . " -xopenmp";
                }
            }
            # sparc need big-endian flag, and can also parallelize advection
            $openmp_switch = "-Domp_on";
            $endian_switch = "-DBigendian";
        }
    }

# ----------------------------------------------------------------------
# CEES-TOOL
# ----------------------------------------------------------------------

    if ($machine eq "t") {
        if ($comp eq "ifort") {
            $lib_path = "/usr/local/netcdf/lib";
            $incl_path = "/usr/local/netcdf/include";
            $optlink = "-L$(lib_path) -lnetcdf -L/usr/local/INTEL/mkl/9.1/lib/em64t -I/usr/local/INTEL/mkl/9.1/include -lmkl_lapack -lmkl_em64t -lguide -lpthread ";
            if ($opt == 0) {
                $opt1 = "-g -m64 -O0 -extend_source -heap-arrays -CB -traceback -static -fpe:0 -assume 2underscores";
            } elsif ($opt == 1) {
                $opt1 = "-g -m64 -O0 -extend_source -heap-arrays -traceback -static -fpe:0 -assume 2underscores";
            } elsif ($opt == 2) {
                $opt1 = "-g -m64 -O1 -extend_source -heap-arrays -traceback -static -fpe:0 -assume 2underscores";
            } elsif ($opt == 3) {
                $opt1 = "-vec-report0 -O3 -ip -no-prec-div -m64 -extend_source -align -heap-arrays -fpp -static -assume 2underscores";
            }
            if ($omp eq "omp") {
                $opt1 = $opt1 . " -openmp";
            }
        }
    }


}    # End of load_make_parts subroutine



sub create_makefile () {


    if ($comp eq "f95" && $machine eq "s") {
        $comp = "/opt/sun/sunstudio12.1/bin/f95";
    }


	$makefile = <<EOT
lib_path = $lib_path
incl_path = $incl_path
$extra_shit

opt1 = $opt1
opt2 = $opt2

optlink = $optlink

ofile = sia2_globals.o $gpu_make_objects sia2_adv_3_setup.o sia2_adv_3_tracers.o sia2_adv_3_grid_ridge.o sia2_env_create_ice_sub.o sia2_env_edd_solution.o sia2_env.o sia2_netcdf_funcs.o sia2_interpolate_forcing.o sia2_forcing.o sia2_write_output.o sia2_initialize.o sia2_read_parameters.o sia2_read_pj.o sia2_timing.o seaice2.o

seaice2:		\$(ofile)
	$comp -o sia2temp \$(opt1) \$(ofile) \$(optlink) $gpulink
seaice2.o:
	$comp -c \$(opt1) seaice2.f90
$gpu_make_code
sia2_timing.o:
	$comp -c \$(opt1) sia2_timing.f90
sia2_read_pj.o:
	$comp -c \$(opt1) sia2_read_pj.f90
sia2_read_parameters.o:
	$comp -c \$(opt1) sia2_read_parameters.f90
sia2_initialize.o:
	$comp -c \$(opt1) sia2_initialize.f90
sia2_forcing.o:
	$comp -c \$(opt1) sia2_forcing.f90
sia2_interpolate_forcing.o:
	$comp -c \$(opt1) sia2_interpolate_forcing.f90
sia2_adv_3_setup.o:
	$comp -c \$(opt1) sia2_adv_3_setup.f90
sia2_adv_3_tracers.o:
	$comp -c \$(opt1) sia2_adv_3_tracers.f90
sia2_adv_3_grid_ridge.o:
	$comp -c \$(opt1) sia2_adv_3_grid_ridge.f90
sia2_env.o:
	$comp -c \$(opt1) sia2_env.f90
sia2_env_create_ice_sub.o:
	$comp -c \$(opt1) sia2_env_create_ice_sub.f90
sia2_env_edd_solution.o:
	$comp -c \$(opt2) sia2_env_edd_solution.f90
sia2_write_output.o:
	$comp -c \$(opt1) sia2_write_output.f90
sia2_netcdf_funcs.o:
	$comp -c \$(opt1) sia2_netcdf_funcs.f90
sia2_globals.o:
	$comp -c \$(opt1) sia2_globals.f90

pre:
	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) sia2_env_update_skeletal.inc.F sia2_env_update_skeletal.inc.f90

	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) sia2_env_ice_remap.inc.F sia2_env_ice_remap.inc.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) sia2_env_snow_remap.inc.F sia2_env_snow_remap.inc.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) sia2_env_pond_mix3.inc.F sia2_env_pond_mix3.inc.f90

	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) sia2_env_redist_new_ice.inc.F sia2_env_redist_new_ice.inc.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) sia2_env_desalinate_7.inc.F sia2_env_desalinate_7.inc.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) sia2_env_flood.inc.F sia2_env_flood.inc.f90

	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) sia2_env_grid_remap.inc.F sia2_env_grid_remap.inc.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) sia2_env_grid_calcs.inc.F sia2_env_grid_calcs.inc.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) sia2_env_flux_heat.inc.F sia2_env_flux_heat.inc.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) sia2_env_redist.inc.F sia2_env_redist.inc.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) sia2_adv_3_ridge.inc.F sia2_adv_3_ridge.inc.f90

	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) sia2_env_create_ice_sub.F sia2_env_create_ice_sub.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) $gpu_cpp_switch sia2_env_icemodel_2.inc.F sia2_env_icemodel_2.inc.f90

	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) $openmp_switch sia2_edd_gpu_calc.F sia2_edd_gpu_calc.f90

	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) sia2_netcdf_funcs.F sia2_netcdf_funcs.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) sia2_write_output.F sia2_write_output.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) sia2_env.F sia2_env.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) $openmp_switch sia2_adv_3_grid_ridge.F sia2_adv_3_grid_ridge.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) sia2_adv_3_setup.F sia2_adv_3_setup.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) $openmp_switch sia2_adv_3_tracers.F sia2_adv_3_tracers.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) sia2_forcing.F sia2_forcing.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) $endian_switch sia2_read_pj.F sia2_read_pj.f90

$gpu_cpp_code

	/opt/local/bin/cpp-mp-5 -w -P -C -I\$(incl_path) $gpu_cpp_switch $comp_ifort seaice2.F seaice2.f90

clean:
	rm -f *.o *.mod seaice2.f90 sia2_env_flood.inc.f90 sia2_env_redist_new_ice.inc.f90 sia2_env_redist.inc.f90 sia2_adv_3_ridge.inc.f90 sia2_env_pond_mix3.inc.f90 sia2_env_create_ice_sub.f90 sia2_env_desalinate_7.inc.f90 sia2_netcdf_funcs.f90 sia2_write_output.f90 sia2temp sia2_env.f90 sia2_env_icemodel_2.inc.f90 sia2_adv_3_grid_ridge.f90 sia2_adv_3_tracers.f90 sia2_adv_3_setup.f90 sia2_env_create_ice.inc.f90 sia2_env_flux_heat.inc.f90 sia2_env_grid_remap.inc.f90 sia2_env_grid_calcs.inc.f90 sia2_env_ice_remap.inc.f90 sia2_env_snow_remap.inc.f90 sia2_env_desalinate_01.inc.f90 sia2_env_desalinate_2.inc.f90 sia2_env_desalinate_3.inc.f90 sia2_forcing.f90 sia2_read_pj.f90 sia2_env_update_skeletal.inc.f90 sia2_edd_gpu_calc.f90

EOT
;

}


