############################################################
# Create a library
############################################################
#增加调试信息
add_definitions("-Wall -g")

include_directories("/dssg02/home/lvws/programs/miniconda3/envs/conda_gcc_env//include")
include_directories("/dssg02/home/lvws/programs/htslib/include")

#Generate the shared library from the library sources
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
        ${PROJECT_SOURCE_DIR}/include
)

target_include_directories(bamread_library
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

target_include_directories(fusion_library
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
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

# link the new hello_library target with the hello_binary target
target_link_libraries( sv-detector
    PRIVATE 
        readfasta_library
        bamread_library
        fusion_library
        # mapping_library
)


target_link_libraries(sv-detector "/dssg02/home/lvws/programs/htslib/lib/libhts.so")
target_link_libraries(sv-detector "/dssg02/home/lvws/programs/miniconda3/envs/conda_gcc_env/lib/libboost_serialization.so")
target_link_libraries(sv-detector z)