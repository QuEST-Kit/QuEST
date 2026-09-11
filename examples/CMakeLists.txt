# @author Oliver Thomson Brown
# @author Erich Essmann (updating rpath)
# @author Tyson Jones (refactored to functions + glob)



# compiles in_fn to ../examples/direc/in_fn_c

function(add_example direc in_fn)

  # fn.cpp -> filename=fn, filelang=cpp
  get_filename_component(filename ${in_fn} NAME_WE)
  get_filename_component(filelang ${in_fn} LAST_EXT)
  string(REPLACE "." "" filelang ${filelang})

  set(target "${filename}_${filelang}")
  set(out_fn "${filename}_${filelang}")
  set(out_dir "${CMAKE_INSTALL_BINDIR}/examples/${direc}/")
  
  add_executable(${target} ${in_fn})
  target_link_libraries(${target} PUBLIC QuEST)


  if (QUEST_ENABLE_INSTALL AND QUEST_INSTALL_BINARIES)
    install(
      TARGETS ${target}
      RUNTIME
      DESTINATION ${out_dir}
      COMPONENT Examples
    )
  endif ()

  set_target_properties(${target} 
    PROPERTIES 
      OUTPUT_NAME "${out_fn}"
  )

  setup_quest_rpath(${target} "${out_dir}")

endfunction()



# compiles all .c or .cpp files in the subdirecs below

function(add_all_local_examples)

  get_filename_component(direc ${CMAKE_CURRENT_LIST_DIR} NAME)

  file(GLOB files
    "${CMAKE_CURRENT_SOURCE_DIR}/*.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
  )

  foreach(fn ${files})
    add_example(${direc} ${fn})
  endforeach()

endfunction()



# all subdirecs below will invoke add_all_local_examples()

add_subdirectory(isolated)
add_subdirectory(extended)
add_subdirectory(automated)
