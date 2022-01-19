############################################################
# Create a library
############################################################

include_directories("/dssg02/home/lvws/programs/miniconda3/envs/conda_gcc_env//include")
#Generate the shared library from the library sources
add_library(readfasta_library SHARED 
    reference.cpp
)
#add_library(bam::library ALIAS bam_library)

target_include_directories(readfasta_library
    PUBLIC 
        ${PROJECT_SOURCE_DIR}/include
)



############################################################
# Create an executable
############################################################

# Add an executable with the above sources
add_executable(test-index
    test-index.cpp
)

# link the new hello_library target with the hello_binary target
target_link_libraries( test-index
    PRIVATE 
        readfasta_library

)

target_link_libraries(test-index "/dssg02/home/lvws/programs/miniconda3/envs/conda_gcc_env/lib/libboost_serialization.so")
#target_link_libraries(readRefer z)