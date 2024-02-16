set(HEADERS 
	all_files/Tridiagonal.hpp
)

set(SOURCES
	all_files/Tridiagonal.cpp
)

add_library(${PROJECT_NAME} STATIC ${HEADERS} ${SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)
