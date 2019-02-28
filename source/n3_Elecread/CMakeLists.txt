#########################################################
#                  Garfield Simulation                  #
#               Author : LAGARDE Francois               #
#########################################################

add_executable(GasFile2 gasfile.cpp)
target_link_libraries(GasFile2 ${ROOT_LIBRARIES} Garfield)
install(TARGETS GasFile2 RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)

add_executable(Elecread Elecread.cpp)
target_link_libraries(Elecread ${ROOT_LIBRARIES} Garfield)
install(TARGETS Elecread RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)
