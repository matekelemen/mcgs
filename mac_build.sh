#!/bin/bash
# authors: Máté Kelemen
# Run this script with the -h flag for info.

# Name of this script
scriptName="$(basename ${BASH_SOURCE[0]})"

# Function for printing usage info
print_help() {
    echo "$scriptName - Configure, build, and install MCGS."
    echo "-h                    : print this help and exit"
    echo "-C                    : clean build and install directories, then exit"
    echo "-b build_path         : path to the build directory (created if it does not exist yet)"
    echo "-i install_path       : path to the install directory (created if it does not exist yet)"
    echo "-t buildType         : build type [FullDebug, Debug, Release, RelWithDebInfo] (Default: Release)"
    echo
    echo "This script provides a build environment for MCGS targeting systems running on Apple Silicon."
    echo "The interface is minimal, so users seeking more control over the build process are invited to tweak this"
    echo "script to suit their requirements better."
    echo
    echo "Build requirements:"
    echo " - Homebrew"
    echo " - CMake (can be installed from homebrew via 'brew install cmake')"
    echo " - LLVM (can be installed from homebrew via 'brew install llvm')"
    echo
    echo "Caveats:"
    echo "The clang that gets shipped by default lacks OpenMP binaries, so one option is to use another version"
    echo "of the compiler without this shortcoming. Therefore, this script relies on LLVM installed from Homebrew."
    echo "Note that the version of LLVM provided by Homebrew is likely different from that of the system, and is not"
    echo "put on the PATH to avoid braking the standard build system."
}

# Utility variables
# Path to the directory containing this script
scriptDir="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
sourceDir="$scriptDir"

toolchainRoot=""                        # <== root path of the compiler package (llvm)
toolchainBin=""                         # <== directory containing compiler executables
toolchainLib=""                         # <== directory containing compiler libraries
toolchainInclude=""                     # <== directory containing compiler headers

generatorTarget="Unix Makefiles"        # <== name of the generator program in CMake
ccacheFlag=""                           # <== sets CXX_COMPILER_LAUNCHER in CMake to ccache if available

# Define default arguments
buildType="Release".                    # <== passed to CMAKE_BUILD_TYPE
buildDir="${sourceDir}/build"           # <== path to the build directory
installDir="${sourceDir}/install"       # <== path to install MCGS to
clean=0                                 # <== clean the build and install directories, then exit

# Parse command line arguments
while getopts "hCb:i:t:" arg; do
    case "$arg" in
        h)  # Print help and exit without doing anything
            print_help
            exit 0
            ;;
        C)  # Set clean flag
            clean=1
            ;;
        b)  # Set build directory
            buildDir="$OPTARG"
            ;;
        i)  # Set install directory
            installDir="$OPTARG"
            ;;
        t)  # Set build type
            buildType="$OPTARG"
            (("${buildType}" == "Debug" || "${buildType}" == "RelWithDebInfo" || "${buildType}" == "Release")) || (print_help && echo "Error: invalid build type: ${buildType}" && exit 1)
            ;;
        \?) # Unrecognized argumnet
            echo "Error: unrecognized argument: $arg"
            exit 1
    esac
done

# Check write access to the build directory
if [ -d "$buildDir" ]; then
    if ! [[ -w "$buildDir" ]]; then
        echo "Error: user '$(hostname)' has no write access to the build directory: '$buildDir'"
        exit 1
    fi
fi

# Check write access to the install dir
if [ -d "$installDir" ]; then
    if ! [[ -w "$installDir" ]]; then
        echo "Error: user '$(hostname)' has no write access to the install directory: '$installDir'"
        exit 1
    fi
fi

# If requested, clear build and install directories, then exit
if [ $clean -ne 0 ]; then
    if [ -d "$buildDir" ]; then
        for item in "$buildDir"; do
            rm -rf "$item"
        done
    fi
    if [ -d "$installDir/mcgs" ]; then
        for item in "$installDir/mcgs"; do
            rm -rf "$item"
        done
    fi
    exit 0
fi

# Check whether CMake is available
if ! command -v cmake &> /dev/null; then
    echo "Error: $scriptName requires CMake"
    echo "Consider running 'brew install cmake'"
    exit 1
fi

# OpenMP is not available on the default clang
# that Apple ships with its system, so another
# toolchain must be used.
# => this script relies on Homebrew to install
# and find necessary packages (such as llvm).
if ! command -v brew &> /dev/null; then
    echo "Error: $scriptName requires Homebrew"
    exit 1
fi

checkHomebrewPackage() {
    if ! brew list "$1" >/dev/null 2>&1; then
        echo "Error: missing dependency: $1"
        echo "Consider running 'brew install $1'"
        exit 1
    fi
}

# Check whether LLVM is installed, and populate related paths
checkHomebrewPackage llvm
toolchainRoot="$(brew --prefix llvm)"
toolchainBin="${toolchainRoot}/bin"
toolchainLib="${toolchainRoot}/lib"
toolchainInclude="${toolchainRoot}/include"

checkOptionalHomebrewPackage() {
    if ! command -v "$1" >/dev/null 2>&1; then
        return 1
    fi
    return 0
}

# Optional dependency - ninja
if checkOptionalHomebrewPackage ninja; then
    generatorTarget="Ninja"
fi

# Create the build directory if it does not exist yet
if [ ! -d "$buildDir" ]; then
    mkdir -p "$buildDir"
fi

cmakeArguments=()
for argument in "-H$sourceDir"                                          \
                "-B$buildDir"                                           \
                "-DCMAKE_INSTALL_PREFIX:STRING=$installDir"             \
                "-G$generatorTarget"                                    \
                "-DCMAKE_BUILD_TYPE:STRING=${buildType}"                \
                "-DCMAKE_C_COMPILER:STRING=${toolchainBin}/clang"       \
                "-DCMAKE_CXX_COMPILER:STRING=${toolchainBin}/clang++"   \
                "-DOPENMP_LIBRARIES:STRING=${toolchainLib}"             \
                "-DOPENMP_INCLUDES:STRING=${toolchainInclude}"          \
                "-DOPENMP_C:STRING=${toolchainBin}/clang"               \
                "-DOPENMP_CXX:STRING=${toolchainBin}/clang++"           \
                "-DCMAKE_COLOR_DIAGNOSTICS:BOOL=ON"                     \
                "-DMCGS_BUILD_TESTS=ON"                                \
                "$ccacheFlag"; do
    if ! [ -z "$argument" ]; then
        cmakeArguments+=("$argument")
    fi
done

# Configure
if ! cmake $(for argument in "${cmakeArguments[@]}"; do echo "$argument"; done); then
    exit $?
fi

# Build and install
if ! cmake --build "$buildDir" --target install -j; then
    exit $?
fi

exit 0
