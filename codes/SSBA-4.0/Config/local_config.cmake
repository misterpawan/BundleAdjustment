# Other optional libraries
enable_feature (V3DLIB_ENABLE_SUITESPARSE)
enable_feature_inc_path (V3DLIB_ENABLE_SUITESPARSE ~/SuiteSparse/include)
enable_feature_lib_path (V3DLIB_ENABLE_SUITESPARSE ~/SuiteSparse/lib)
enable_feature_libraries (V3DLIB_ENABLE_SUITESPARSE colamd umfpack cxsparse)

#Intel MKL libraries
enable_feature (ENABLE_INTEL_MKL)
enable_feature_inc_path (ENABLE_INTEL_MKL /opt/intel/compilers_and_libraries_2019.3.199/linux/mkl/include/)
enable_feature_lib_path (ENABLE_INTEL_MKL /opt/intel/mkl/lib/intel64_lin)
#enable_feature_libraries (ENABLE_INTEL_MKL mkl_intel_lp64 mkl_sequential mkl_core pthread m dl)
#for using OpenMP(GNU libgomp)
enable_feature_libraries (ENABLE_INTEL_MKL mkl_intel_lp64 mkl_gnu_thread gomp mkl_core pthread m dl)

# Debug/release mode selection
set (CMAKE_BUILD_TYPE Release)
#set (CMAKE_BUILD_TYPE Debug)
