# --------------------------------------------------------------------------- #
#    CMake file for SatellitesBlock/test                                      #
#                                                                             #
#    In this file we define the unit tests for this library.                  #
#    Tests that should be supported by ctest executable must be added         #
#    with the add_test() command.                                             #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# Since we are using the block factory, objects from linked libraries
# may wrongly appear unused, and by default the linker does not include them,
# so we have to force the linking.
if (BUILD_SHARED_LIBS)
    if (UNIX AND (NOT APPLE))
        add_link_options("-Wl,--no-as-needed")
    endif ()
else ()
    if (MSVC)
        add_link_options("/WHOLEARCHIVE" "/FORCE:MULTIPLE")
    else () # Unix
        if (APPLE)
            add_link_options("-Wl,-all_load")
        else ()
            add_link_options("-Wl,--whole-archive,--allow-multiple-definition")
        endif ()
    endif ()
endif ()

# ----- SatellitesBlock_test -------------------------------------------------- #

add_executable(${modName}_test test.cpp)
target_link_libraries(${modName}_test PRIVATE ${modNamespace}::${modName})

if (WIN32 AND BUILD_SHARED_LIBS AND NOT DEFINED ENV{CONDA_BUILD})
    # On Windows with shared libraries enabled copies all runtime DLL
    # dependencies next to the executable in the build directory
    # (POST_BUILD), using CMake's TARGET_RUNTIME_DLLS generator expression
    add_custom_command(TARGET ${modName}_test POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_RUNTIME_DLLS:${modName}_test>
            $<TARGET_FILE_DIR:${modName}_test>
            COMMAND_EXPAND_LISTS
            COMMENT "Copying runtime DLLs next to ${modName}_test.exe")
endif ()

add_test(NAME ${modName}_test
         COMMAND ${modName}_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# The label lets the module CI select exactly the tests exercising this
# module with `ctest -L <module>` (see tests/cmake/TestLabels.cmake in the
# umbrella for the system-tests counterpart).
set_tests_properties(${modName}_test PROPERTIES LABELS "${modName}")

# --------------------------------------------------------------------------- #

# ----- SatellitesBlock_milp_test -------------------------------------------- #
# solves a small SCDP instance as a MILP; it needs a MILPSolver, so it is
# built only when the module is available
if (NOT TARGET ${modNamespace}::MILPSolver)
    find_package(MILPSolver QUIET)
endif ()
if (TARGET ${modNamespace}::MILPSolver)
    add_executable(${modName}_milp_test test_milp.cpp)
    target_link_libraries(${modName}_milp_test PRIVATE
            ${modNamespace}::${modName} ${modNamespace}::MILPSolver)

    add_test(NAME ${modName}_milp_test
             COMMAND ${modName}_milp_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    set_tests_properties(${modName}_milp_test PROPERTIES LABELS "${modName}")
endif ()

# ----- SatellitesBlock_ld_test ---------------------------------------------- #
# solves a small SCDP instance by Lagrangian decomposition, so it needs
# LagrangianDualSolver, BundleSolver and a *MILPSolver; it is built only
# when all the modules are available
# BundleSolver may have been built with libtorch (LearningBundle): in that
# case its imported target references torch's own CMake targets, that only
# exist after finding the Torch package
find_package(Torch QUIET)
if (NOT TARGET ${modNamespace}::BundleSolver)
    find_package(BundleSolver QUIET)
endif ()
if (NOT TARGET ${modNamespace}::LagrangianDualSolver)
    find_package(LagrangianDualSolver QUIET)
endif ()
if (TARGET ${modNamespace}::BundleSolver AND
    TARGET ${modNamespace}::LagrangianDualSolver AND
    TARGET ${modNamespace}::MILPSolver)
    add_executable(${modName}_ld_test test_ld.cpp)
    target_link_libraries(${modName}_ld_test PRIVATE
            ${modNamespace}::${modName}
            ${modNamespace}::LagrangianDualSolver
            ${modNamespace}::BundleSolver
            ${modNamespace}::MILPSolver)

    add_test(NAME ${modName}_ld_test
             COMMAND ${modName}_ld_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    set_tests_properties(${modName}_ld_test PROPERTIES LABELS "${modName}")
endif ()
