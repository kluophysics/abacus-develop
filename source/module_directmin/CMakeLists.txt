
list(APPEND objects

    composite.cpp 
    directmin.cpp

    linesearch_options.cpp
    linesearch_helper.cpp
    linesearch_base.cpp

    occupation.cpp
    options.cpp
    
    problem.cpp

    stiefel.cpp

    # optimizer/options_ls.cpp
    # optimizer/base.cpp
    # optimizer/ls_base.cpp
    # optimizer/ls_rcg.cpp
    # optimizer/ls_rsd.cpp
    # optimizer/ls_rbfgs.cpp
    
    # objective/problem.cpp
    # objective/problem_test.cpp
    # objective/problem_eig.cpp
    # objective/problem_ks.cpp
    # manifold/stiefel.cpp
    )

find_package(Armadillo REQUIRED)


add_library(
    directmin
    OBJECT
    ${objects}
)
target_link_libraries(directmin
  ${ARMADILLO_LIBRARIES}
)
# find_package(FFTW3 REQUIRED)
# find_package(LAPACK REQUIRED)
# include_directories(${FFTW3_INCLUDE_DIRS})

# target_include_directories(directmin ROPTLIB)

# set(DEFAULT_COMPILE_OPTIONS ${DEFAULT_COMPILE_OPTIONS} -O3 -ffastmath -march=native -ggdb3
# -std=c++0x -Wall -g -DDRIVERCPP -DROPTLIB_WITH_FFTW -llapack -lblas -lfftw3 -lm)

# target_link_libraries(directmin 
#   PRIVATE
#     math_libs
# )


if(ENABLE_COVERAGE)
  add_coverage(directmin)
endif()

# if(BUILD_TESTING)
#   add_subdirectory(test)
# endif()