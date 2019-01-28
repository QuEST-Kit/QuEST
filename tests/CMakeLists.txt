#add_custom_target(tests
#    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/example.sh
#    DEPENDS QuEST
#)

add_custom_target(tests
    COMMAND 
        ${CMAKE_CURRENT_SOURCE_DIR}/QuESTTest.py -Q ${QuEST_LIB_PATH} -p ${CMAKE_CURRENT_SOURCE_DIR}/unitPy
    DEPENDS 
        QuEST
)
