project(mpc)
set(CMAKE_BUILD_TYPE Debug)
cmake_minimum_required (VERSION 3.5)

# -g allows for gdb debugging
# turn on -03 for best performance
add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS "${CXX_FLAGS}")

# set(sources src/custom_MPC.cpp src/test_custom.cpp)
# set(sources src/Ipopt_MPC.cpp src/test_ipopt_MPC.cpp)
# set(sources src/Ipopt_MPC_one_step.cpp src/test_ipopt_MPC_single_step.cpp)
# set(sources src/ipopt_without_cppad_example.cpp)
# set(sources src/Ipopt_MPC_naive.cpp src/test_naive_ipopt_MPC.cpp)
# set(sources src/apply_ipopt_example.cpp)
set(sources src/apply_ipopt_kinematic_model.cpp)

include_directories(src/Eigen-3.3)
include_directories(/usr/local/include/coin)
link_directories(/usr/local/lib)

add_executable(mpc ${sources})

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
target_link_libraries(mpc ipopt python2.7)
else()
target_link_libraries(mpc ipopt python)
endif()
