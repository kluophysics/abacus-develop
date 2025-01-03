list(APPEND objects
# Manifolds/CFixedRankQ2F.cpp
Manifolds/CStiefel.cpp
# Manifolds/CSymFixedRankQ.cpp
Manifolds/Element.cpp
Manifolds/Euclidean.cpp
# Manifolds/FixedRankE.cpp
# Manifolds/FixedRankQ2F.cpp
Manifolds/Grassmann.cpp
Manifolds/Manifold.cpp
Manifolds/MultiManifolds.cpp
# Manifolds/SPDManifold.cpp
Manifolds/SmartSpace.cpp
# Manifolds/Sphere.cpp
# Manifolds/SphereTx.cpp
Manifolds/Stiefel.cpp
# Manifolds/SymFixedRankQ.cpp 

Others/BlasLapackCppWrapper.cpp
Others/ForDebug.cpp
# Others/MinPNormConHull.cpp
Others/SparseMatrix.cpp
Others/Spline.cpp
Others/Timer.cpp
Others/randgen.cpp 
Others/SparseBLAS/nist_spblas.cpp 
Others/wavelet/wavelet.cpp 

# Solvers/AManPG.cpp 
Solvers/LRBFGS.cpp 
# Solvers/LRBFGSSub.cpp 
Solvers/LRTRSR1.cpp 
# Solvers/ManPG.cpp 
Solvers/RBFGS.cpp 
# Solvers/RBFGSSub.cpp 
Solvers/RBroydenFamily.cpp 
Solvers/RCG.cpp 
# Solvers/RGS.cpp 
Solvers/RNewton.cpp 
Solvers/RSD.cpp 
Solvers/RTRNewton.cpp 
Solvers/RTRSD.cpp 
Solvers/RTRSR1.cpp 
Solvers/RWRBFGS.cpp 
Solvers/Solvers.cpp 
# Solvers/SolversNSM.cpp 
# Solvers/SolversNSMPGLS.cpp 
# Solvers/SolversNSMSub.cpp 
# Solvers/SolversNSMSubLS.cpp 
Solvers/SolversSMLS.cpp 
Solvers/SolversSMTR.cpp 
Solvers/SolversSM.cpp 

directmin.cpp
    )

# find_package(Armadillo REQUIRED)

 add_library(
    directmin
    OBJECT
    ${objects}
)

target_link_libraries(directmin
  # ${ARMADILLO_LIBRARIES}
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(
cwrapper/blas/ 
cwrapper/lapack/
Others/SparseBLAS/
Others/fftw/
Others/wavelet/
)

target_include_directories(directmin PUBLIC
    #  Manifolds
    #  Others
    #  Problems
    #  Solvers
)

if(ENABLE_COVERAGE)
  add_coverage(directmin)
endif()