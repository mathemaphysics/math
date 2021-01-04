add_executable(test46 test46.c mbuild-nsymm.c funcs.c)
target_link_libraries(test46 math m)
target_include_directories(test46 PUBLIC ${PROJECT_SOURCE_DIR})