import os
import multiprocessing
import tools.os
import tools.git


def do_install(src_path, build_path, install_path, options):
    src_path = os.path.abspath(src_path)
    build_path = os.path.abspath(build_path)
    install_path = os.path.abspath(install_path)

    if not os.path.isdir(build_path):
        os.makedirs(build_path)

    tools.git.prepare_submodule(src_path)
    old_cwd = os.getcwd()
    os.chdir(build_path)

    joined_options = " ".join([f"--prefix={install_path}"] + options)

    # if os.path.exists(os.path.join(src_path, "autogen.sh")):
    #     tools.os.run(
    #         f"{os.path.join(src_path, "autogen.sh")} {joined_options}", "generate"
    #     )

    tools.os.run(f"{os.path.join(src_path, "configure")} {joined_options}", "configure")
    tools.os.run(f"make -j{multiprocessing.cpu_count()}", "compile")
    tools.os.run(f"make install", "install")

    os.chdir(old_cwd)
