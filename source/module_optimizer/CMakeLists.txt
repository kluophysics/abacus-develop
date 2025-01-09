list(APPEND objects
# manifolds/matrix.cpp
manifolds/stiefel.cpp
manifolds/euclidean.cpp
# manifolds/stiefel_manifold.cpp

manifolds/manifold.cpp
# manifolds/product_manifold.cpp
problems/problem.cpp


optimizer_base.cpp
optimizer_ls_base.cpp
  )

find_package(Armadillo REQUIRED)

 add_library(
    optimizer
    OBJECT
    ${objects}
)

target_link_libraries(optimizer
  ${ARMADILLO_LIBRARIES}
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(
  manifolds/
# cwrapper/blas/ 
# cwrapper/lapack/
# Others/SparseBLAS/
# Others/fftw/
# Others/wavelet/
)

target_include_directories(optimizer PUBLIC
    #  Manifolds
    #  Others
    #  Problems
    #  Solvers
)

if(ENABLE_COVERAGE)
  add_coverage(optimizer)
endif()