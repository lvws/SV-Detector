############################################################
# Create a library
############################################################


#增加调试信息  ${PROJECT_SOURCE_DIR}
add_definitions("-Wall -g")

include_directories("./boost/include")
include_directories("./htslib/include")

#Generate the SHARED library from the library sources
add_library(readfasta_library SHARED
    src/reference.cpp
)

add_library(bamread_library SHARED
    src/BamReader.cpp
)

add_library(fusion_library SHARED
    src/fusion.cpp
)

# add_library(mapping_library SHARED
#     src/mapping.cpp
# )

#add_library(bam::library ALIAS bam_library)
target_include_directories(readfasta_library
    PUBLIC 
        ./include
)

target_include_directories(bamread_library
    PUBLIC
        ./include
)

target_include_directories(fusion_library
    PUBLIC
        ./include
)

# target_include_directories(mapping_library
#     PUBLIC
#         ${PROJECT_SOURCE_DIR}/include
# )



############################################################
# Create an executable
############################################################

# Add an executable with the above sources
add_executable(sv-detector
    src/main.cpp
)

#使用相对路径
# Find shared libraries next to the executable
set_target_properties( sv-detector PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,$ORIGIN/")

# link the new hello_library target with the hello_binary target
target_link_libraries( sv-detector
    PRIVATE 
        readfasta_library
        bamread_library
        fusion_library
        # mapping_library
)


target_link_libraries(sv-detector "libhts.so")
target_link_libraries(sv-detector "libboost_serialization.so")
target_link_libraries(sv-detector z)