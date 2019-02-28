#########################################################
#                  Garfield Simulation                  #
#               Author : LAGARDE Francois               #
#########################################################

add_executable(Ionization Ionization.cpp)
target_link_libraries(Ionization ${ROOT_LIBRARIES} Garfield)
install(TARGETS Ionization RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)
