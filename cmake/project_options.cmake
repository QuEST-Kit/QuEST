function(project_options_setup target)
    # C and C++ standards
    target_compile_features(${target} INTERFACE cxx_std_20)
    target_compile_features(${target} INTERFACE c_std_17)

    # Identifying compilers and languages
    set(clang_variants "AppleClang,ARMClang,Clang,CrayClang,FujitsuClang,IntelLLVM,XLClang,IBMClang")

    # C and CXX compiler warnings
    set(msvc_warnings
            /W4
            /permissive-
    )
    set(llvm_warnings
            -Wall
            -Wextra
            -Wshadow
            -Wpedantic
            -Wno-unused-parameter
    )
    set(gcc_warnings
            ${llvm_warnings}
            -Wlogical-op
            $<$<BOOL:${ENABLE_MULTITHREADING}>:-Wno-unknown-pragmas>
    )

    foreach(language IN ITEMS C CXX)
        target_compile_options(${target}
                INTERFACE
                $<$<COMPILE_LANG_AND_ID:${language},GNU>:$<BUILD_INTERFACE:${gcc_warnings}>>
                $<$<COMPILE_LANG_AND_ID:${language},MSVC>:$<BUILD_INTERFACE:${msvc_warnings}>>
                $<$<COMPILE_LANG_AND_ID:${language},${clang_variants}>:$<BUILD_INTERFACE:${llvm_warnings}>>
        )
    endforeach()

    # CUDA and HIP Warnings
    set(cuda_warnings
            -Wall
            -Wextra
            -Wunused
            -Wconversion
            -Wshadow
    )

    target_compile_options(
            ${target}
            INTERFACE
            $<$<COMPILE_LANGUAGE:CUDA,HIP>:$<BUILD_INTERFACE:${cuda_warnings}>>
    )

    # Complex Arithmetic Option
    foreach(language IN ITEMS C CXX)
        target_compile_options(${target}
                INTERFACE
                $<$<COMPILE_LANG_AND_ID:${language},GNU>:$<BUILD_INTERFACE:-fcx-method=fortran>>
                $<$<COMPILE_LANG_AND_ID:${language},MSVC>:$<BUILD_INTERFACE:/fp:precise>>
                $<$<COMPILE_LANG_AND_ID:${language},${clang_variants}>:$<BUILD_INTERFACE:-fcomplex-arithmetic=improved>>
        )
    endforeach()

endfunction(project_options_setup)