list(APPEND objects
    optim.h
    )

find_package(Armadillo REQUIRED)

 add_library(
    optim
    OBJECT
    ${objects}
)

target_link_libraries(optim
  ${ARMADILLO_LIBRARIES}
)

if(ENABLE_COVERAGE)
  add_coverage(optim)
endif()