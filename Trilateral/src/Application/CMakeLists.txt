file(GLOB_RECURSE INTERNAL
     Include/*.h
     Source/*.cpp
)
include_directories( ../External/Include)
add_library(INTERNAL_SOURCE ${INTERNAL})
