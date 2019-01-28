add_custom_target(tests
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/example.sh
    DEPENDS QuEST
)
