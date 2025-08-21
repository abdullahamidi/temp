mkdir build
cd build
cmake -DTIOGA_HAS_NODEGID:BOOL=off -DTIOGA_ENABLE_TIMERS:BOOL=on ..
make
cd ..
cd case/;./run.sh 8;cd -
