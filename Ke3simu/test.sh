#! /bin/bash
root -l <<EOF
.L StHelixD.cxx+
.L StPhysicalHelixD.cxx+
.L simv2.C++
simv2(1e4, 0, 0 ,"out/test.gamma")
.q
EOF  
