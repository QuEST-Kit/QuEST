add_subdirectory(src)

add_library(QuEST_lib
    ${QuEST_SRC}
)

target_include_directories(QuEST_lib 
    PRIVATE src
    PUBLIC include
)

target_compile_options(
    QuEST_lib 
    PRIVATE
    "-std=c99"
)


