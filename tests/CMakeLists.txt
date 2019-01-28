#add_custom_target(tests
#    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/example.sh
#    DEPENDS QuEST
#)

add_custom_target(tests
    COMMAND 
        python ${CMAKE_CURRENT_SOURCE_DIR}/example.py ${QuEST_LIB_PATH}
    DEPENDS 
        QuEST
)
