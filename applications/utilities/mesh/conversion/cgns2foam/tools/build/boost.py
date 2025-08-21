import os
import multiprocessing
import shutil
import tools.os
import tools.git


def do_install(src_path, build_path, install_path, options=[]):
    src_path = os.path.abspath(src_path)
    build_path = os.path.abspath(build_path)
    install_path = os.path.abspath(install_path)
    
    if not os.path.isdir(build_path):
        os.makedirs(build_path)

    tools.git.prepare_submodule(src_path)

    old_cwd = os.getcwd()
    os.chdir(build_path)
    if not os.path.isfile(f"./b2{tools.os.EXE_EXT}"):
        print("Bootstrapping boost...")
        tools.os.run(
            f'{src_path}{tools.os.PATH_SEP}bootstrap{tools.os.SHELL_EXT}',
            "bootstrap"
        )

    build_options = ' '.join([
        f'-j{multiprocessing.cpu_count()} ',
        f'--prefix={install_path}'
    ] + options)

    
    os.chdir(src_path)
    tools.os.run(
        f"{build_path}{tools.os.PATH_SEP}b2{tools.os.EXE_EXT} --build-dir={build_path}/.. headers {build_options}",
        "build"
    )
    tools.os.run(
        f"{build_path}{tools.os.PATH_SEP}b2{tools.os.EXE_EXT} --build-dir={build_path}/.. install {build_options}",
        "install"
    )

    os.chdir(old_cwd)
