#########################################################
#                  Garfield Simulation                  #
#               Author : LAGARDE Francois               #
#########################################################

add_executable(GasFile gasfile.cpp)
target_link_libraries(GasFile ${ROOT_LIBRARIES} Garfield)
install(TARGETS GasFile RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)

add_executable(TownsendCof TownsendCof.cpp)
target_link_libraries(TownsendCof ${ROOT_LIBRARIES} Garfield)
install(TARGETS TownsendCof RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)

add_executable(gasfile_new gasfile_new.cpp)
target_link_libraries(gasfile_new ${ROOT_LIBRARIES} Garfield)
install(TARGETS gasfile_new RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)
