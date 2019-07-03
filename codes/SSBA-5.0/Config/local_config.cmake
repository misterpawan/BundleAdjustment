# Other optional libraries
enable_feature (V3DLIB_ENABLE_SUITESPARSE)
enable_feature_inc_path (V3DLIB_ENABLE_SUITESPARSE /usr/include/suitesparse)
#enable_feature_inc_path (V3DLIB_ENABLE_SUITESPARSE ~/SuiteSparse/include)

#enable_feature_lib_path (V3DLIB_ENABLE_SUITESPARSE ~/SuiteSparse/lib)
enable_feature_lib_path (V3DLIB_ENABLE_SUITESPARSE /usr/lib/x86_64-linux-gnu)
enable_feature_libraries (V3DLIB_ENABLE_SUITESPARSE colamd umfpack cxsparse)

#Intel MKL libraries
enable_feature (ENABLE_INTEL_MKL)
enable_feature_inc_path (ENABLE_INTEL_MKL /home/siddhant/intel/compilers_and_libraries_2019.4.243/linux/mkl/include/)
enable_feature_lib_path (ENABLE_INTEL_MKL /home/siddhant/intel/mkl/lib/intel64_lin)
enable_feature_libraries (ENABLE_INTEL_MKL mkl_intel_lp64 mkl_sequential mkl_core pthread m dl)

# Debug/release mode selection
set (CMAKE_BUILD_TYPE Release)
#set (CMAKE_BUILD_TYPE Debug)
