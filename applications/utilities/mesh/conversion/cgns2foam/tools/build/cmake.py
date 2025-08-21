import os
import multiprocessing
import tools.os
import tools.git


def do_install(
    src_path, build_path, install_path, generator=None, wrapper="", debug=True, **kwargs
):
    src_path = os.path.abspath(src_path)
    build_path = os.path.abspath(build_path)
    install_path = os.path.abspath(install_path)

    if not os.path.isdir(build_path):
        os.makedirs(build_path)

    tools.git.prepare_submodule(src_path)

    options = []
    if generator:
        options.append(f"-G{generator}")

    for key, value in kwargs.items():
        options.append(f"-D{key}={value}")

    options.append(f"-DCMAKE_INSTALL_PREFIX={install_path}")
    if debug:
        options.append("-DCMAKE_DEBUG_POSTFIX=d")

    joined_options = " ".join(options)

    tools.os.run(f"{wrapper} cmake {src_path} {joined_options}", "configuration", build_path)

    if debug:
        tools.os.run(
            f"cmake --build . --target install --config Debug --parallel {multiprocessing.cpu_count()}",
            "debug install",
            build_path
        )
    tools.os.run(
        f"cmake --build . --target install --config Release --parallel {multiprocessing.cpu_count()}",
        "release install",
        build_path
    )


def package_exists(name, compiler_id="GNU"):
    return (
        os.system(
            f"cmake --find-package -DNAME={name} -DCOMPILER_ID={compiler_id} -DLANGUAGE=CXX -DMODE=EXIST"
        )
        == 0
    )
