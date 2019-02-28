#########################################################
#                  Garfield Simulation                  #
#               Author : LAGARDE Francois               #
#########################################################

add_executable(Signal_correct Signal_correct.cpp)
target_link_libraries(Signal_correct ${ROOT_LIBRARIES} Garfield)
install(TARGETS Signal_correct RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)

add_executable(Elecread Elecread.cpp)
target_link_libraries(Elecread ${ROOT_LIBRARIES} Garfield)
install(TARGETS Elecread RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)
