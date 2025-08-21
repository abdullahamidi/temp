import os
import tools.os


def do_install(deps_path, install_path):
    if os.name == "nt":
        old_cwd = os.getcwd()
        os.chdir(f"{deps_path}/openssl")
        do_install_nt(install_path)
        os.chdir(old_cwd)
    else:
        print("not installing openssl")


def do_install_nt(install_path):
    tools.os.add_env('PATH', 'C:/Program Files/Microsoft Visual Studio/2022/Enterprise/VC/Auxiliary/Build')
    tools.os.add_env('PATH', 'C:/Strawberry/perl/bin')
    tools.os.add_env('PATH', 'C:/Program Files/NASM')

    tools.os.check("vcvars64.bat", "no MSVC environment found")
    tools.os.check("perl", "please install perl")
    tools.os.check("nasm", "please install NASM")
    tools.os.check("nmake", "NMAKE not found")

    tools.os.run(f"vcvars64.bat & perl Configure VC-WIN64 --prefix=\"{install_path}\" & nmake & nmake install", "build and install")
