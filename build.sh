#!/usr/bin/env bash

scriptDir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "$scriptDir"

if [ -d build ]; then rm -rf build; fi

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
