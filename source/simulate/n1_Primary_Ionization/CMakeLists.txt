#########################################################
#                  Garfield Simulation                  #
#               Author : LAGARDE Francois               #
#########################################################

add_executable(Ionization Ionization.C)
target_link_libraries(Ionization ${ROOT_LIBRARIES})
install(TARGETS Ionization RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)
