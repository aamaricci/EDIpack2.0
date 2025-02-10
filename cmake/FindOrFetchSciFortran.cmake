include(FetchContent)
FIND_PACKAGE(PkgConfig)

PKG_SEARCH_MODULE(SCIFOR QUIET scifor)

IF(SCIFOR_FOUND)
  MESSAGE(STATUS "SCIFOR found at: ${SCIFOR_PREFIX}")
ELSE()
    MESSAGE(STATUS "SciFortran not found. Fetching from GitHub...")
    include(ExternalProject)
    
    set(SCIFOR_SOURCE_DIR ${CMAKE_BINARY_DIR}/scifor/src)
    set(SCIFOR_BUILD_DIR ${CMAKE_BINARY_DIR}/scifor/build)
    set(SCIFOR_STAMP_DIR ${CMAKE_BINARY_DIR}/scifor/stamp)
    

    ExternalProject_Add(SciFortran
        PREFIX ${CMAKE_BINARY_DIR}/scifor
        GIT_REPOSITORY https://github.com/SciFortran/SciFortran.git
        GIT_PROGRESS 1
        GIT_TAG master
        SOURCE_DIR ${SCIFOR_SOURCE_DIR}
        BINARY_DIR ${SCIFOR_BUILD_DIR}
        STAMP_DIR ${SCIFOR_STAMP_DIR}
        CMAKE_ARGS 
            -DCMAKE_INSTALL_PREFIX:PATH=${SCIFOR_INSTALL_DIR}
            -DBUILD_SHARED_LIBS=OFF  # Ensure static build
               
        #Trick to override SciFortran's naming convention. We don't need it and it's a mess to make it work in here
        PATCH_COMMAND > ${SCIFOR_SOURCE_DIR}/cmake/MainConfig.cmake && > ${SCIFOR_SOURCE_DIR}/cmake/PostConfig.cmake
        BUILD_COMMAND $(MAKE)  # Parallel build
        INSTALL_COMMAND ""
    )

    # Extract include and lib paths
    set(SCIFOR_INCLUDE_DIRS "${SCIFOR_BUILD_DIR}/include")
    
ENDIF()

INCLUDE_DIRECTORIES(BEFORE ${SCIFOR_INCLUDE_DIRS})
