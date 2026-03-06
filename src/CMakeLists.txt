cmake_minimum_required(VERSION 3.16)
project(MC_simple_stage1 C)

set(DSFMT_MEXP 19937)

option(USE_MPI "Build with MPI support for independent-run parallelism" OFF)

add_executable(MC_simple
    main.c
    memory.c
    input_parser.c
    lattice.c
    mc_update.c
    dSFMT.c
)

if(USE_MPI)
    find_package(MPI REQUIRED C)
    target_compile_definitions(MC_simple PRIVATE DSFMT_MEXP=${DSFMT_MEXP})
    target_link_libraries(MC_simple PRIVATE m MPI::MPI_C)
else()
    target_compile_definitions(MC_simple PRIVATE NO_MPI DSFMT_MEXP=${DSFMT_MEXP})
    target_link_libraries(MC_simple PRIVATE m)
endif()

target_compile_options(MC_simple PRIVATE -O2)
