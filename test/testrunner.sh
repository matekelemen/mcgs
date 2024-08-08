#!/bin/bash

# Set working directory to the <repo-root>/test.
scriptDir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "$scriptDir"

# Check if test files are available and extract them if they're not.
for fileName in mok_structure.mm mok_structure_rhs.mm mok_fluid.mm mok_fluid_rhs.mm; do
    if ! [ -f $fileName ]; then
        if ! tar -xvzf input.tar.gz; then
            exit 1
        fi
        break
    fi
done

# Run tests
cd ..
for shrinkingFactor in 1 8 64 256; do
    for threadCount in 1 2 4; do
        for case in mok_structure mok_fluid; do
            echo "Running $case on $threadCount threads with $shrinkingFactor shrinking ..."
            export OMP_NUM_THREADS=$threadCount
            if ! build/bin/${1}mcgs_testrunner          \
                    test/${case}.mm test/${case}_rhs.mm \
                    -s $shrinkingFactor                 \
                    -i 1e2                              \
                    -o reordered.mm
                then exit 1
            fi
        done
    done
done
