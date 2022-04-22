add_executable(detox main.cpp settings.cpp dbgraph.cpp dsnode.cpp
               kmernpp.cpp util.cpp debug.cpp coverage.cpp alignment.cpp
               graphaln.cpp refcomp.cpp correctgraph.cpp crfmult.cpp
               readaln.cpp)

target_link_libraries(detox readfile kmer pgm pthread libdai gmp)

if (ZLIB_FOUND)
        target_link_libraries(detox ${ZLIB_LIBRARY})
endif (ZLIB_FOUND)

install(TARGETS detox RUNTIME DESTINATION bin)
add_subdirectory(readfile)
add_subdirectory(kmer)
add_subdirectory(pgm)
add_subdirectory(libdai)

