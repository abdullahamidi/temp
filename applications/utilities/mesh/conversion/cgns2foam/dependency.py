import os
import sys


ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.append(ROOT)

import argparse
import shutil
import tools.build.boost
import tools.build.cmake
import tools.git
import tools.os


from dotenv import load_dotenv


def getenv(key):
    result = os.getenv(key)
    if result is None:
        exit(-1)
    return result


def install():
    DEPENDENCY_DIR = getenv("DEPENDENCY_DIR")
    DEPENDENCY_BUILD_DIR = getenv("DEPENDENCY_BUILD_DIR")
    DEPENDENCY_INSTALL_DIR = getenv("DEPENDENCY_INSTALL_DIR")

    DEPENDENCIES_PATH = os.path.abspath(f"{ROOT}/{DEPENDENCY_DIR}")
    BUILD_PATH = os.path.abspath(f"{ROOT}/{DEPENDENCY_BUILD_DIR}")
    INSTALL_PATH = os.path.abspath(f"{ROOT}/{DEPENDENCY_INSTALL_DIR}")

    # if os.path.isdir(BUILD_PATH):
    #     shutil.rmtree(BUILD_PATH)
    # if os.path.isdir(INSTALL_PATH):
    #     shutil.rmtree(INSTALL_PATH)

    tools.os.add_env("PATH", f"{INSTALL_PATH}/bin")
    tools.os.add_env("LD_LIBRARY_PATH", f"{INSTALL_PATH}/lib")

    tools.build.boost.do_install(
        f"{DEPENDENCIES_PATH}/boost",
        f"{BUILD_PATH}/boost",
        INSTALL_PATH,
        [
            "link=static",
            "--with-log",
            "--with-filesystem",
            "--with-thread",
            "--with-system",
            "--with-program_options",
            "--with-json",
        ],
    )

    tools.build.cmake.do_install(
        f"{DEPENDENCIES_PATH}/hdf5",
        f"{BUILD_PATH}/hdf5",
        INSTALL_PATH,
        HDF5_ENABLE_ZLIB_SUPPORT=True,
        # BUILD_SHARED_LIBS=False
    )

    tools.build.cmake.do_install(
        f"{DEPENDENCIES_PATH}/cgns",
        f"{BUILD_PATH}/cgns",
        INSTALL_PATH,
        HDF5_ROOT=INSTALL_PATH,
        # BUILD_SHARED_LIBS=False
    )


def bundle():
    BUNDLE_DIR = getenv("BUNDLE_DIR")

    BUNDLE_PATH = os.path.abspath(f"{ROOT}/{BUNDLE_DIR}")

    if os.path.isdir(BUNDLE_PATH):
        shutil.rmtree(BUNDLE_PATH)

    tools.git.bundle_submodules(ROOT, BUNDLE_PATH)


def unbundle():
    BUNDLE_DIR = getenv("BUNDLE_DIR")

    BUNDLE_PATH = os.path.abspath(f"{ROOT}/{BUNDLE_DIR}")

    tools.git.init_submodules(ROOT, False)
    tools.git.unbundle_submodules(ROOT, BUNDLE_PATH)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="dependency.py",
        description="Install, bundle and unbundle dependencies",
    )
    parser.add_argument("-i", "--install", action="store_true")
    parser.add_argument("-b", "--bundle", action="store_true")
    parser.add_argument("-u", "--unbundle", action="store_true")
    args = parser.parse_args()

    load_dotenv(f"{ROOT}/.config")

    if args.install:
        install()
    if args.bundle:
        bundle()
    if args.unbundle:
        unbundle()
