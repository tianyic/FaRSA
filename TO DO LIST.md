Current C++ code v2 is indeed much faster than v1 by the followings:

-- Load problems lightning fast now by loading sparse format data

-- Define SparseOperator, accelerate matrix operation  

-- Merge some operations in one loop

-- Accelerate sort algorithm


TO DO:

-- Accelerate line search

-- Code seems have some bug, since on some problems, v2 runs more iterations than Matlab implementation. Find and fix it.

-- Write a script to convert datasets scaled in Matlab into the sparse format which can be read by readProblem method in Datastream.h

-- Implement CDsolver()

-- Write make file




