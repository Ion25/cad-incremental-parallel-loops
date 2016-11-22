ExternalProject_Add(
    GTest_EP
	GIT_REPOSITORY https://github.com/google/googletest.git
	GIT_TAG "release-${GTEST_VERSION}"
	INSTALL_COMMAND ""
	BUILD_BYPRODUCTS ${CMAKE_BINARY_DIR}/resources/src/googletest-build/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a ${CMAKE_BINARY_DIR}/resources/src/googletest-build/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main.a
)

ExternalProject_Get_Property(GTest_EP source_dir)
ExternalProject_Get_Property(GTest_EP binary_dir)

add_imported_library(GTESTCORE STATIC "${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a" "${source_dir}/include")
add_imported_library(GTESTMAIN STATIC "${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main.a" "${source_dir}/include")

set(GTEST_LIBRARIES GTESTCORE_STATIC GTESTMAIN_STATIC pthread dl)

add_dependencies(GTESTCORE_STATIC GTest_EP)
add_dependencies(GTESTMAIN_STATIC GTest_EP)
