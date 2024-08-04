#!/usr/bin/env bash

scriptDir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "$scriptDir"

if [ -d build ]; then
    if [ -f build/CMakeCache.txt ]; then rm build/CMakeCache.txt; fi
#    if [ -d build/CMakeFiles ]; then rm -rf build/CMakeFiles; fi
fi

if ! cmake                                      \
    -H.                                         \
    -Bbuild                                     \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON          \
    -DMCGS_BUILD_TESTS=ON                       \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo           \
    -DCMAKE_CXX_FLAGS="-fno-omit-frame-pointer"
    then exit 1
fi

cd build
cmake --build .
cmake --install . --prefix install
