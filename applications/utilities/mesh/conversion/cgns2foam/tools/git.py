import configparser
import os
import tools.os
import uuid


def is_git_dir(path):
    if os.path.isdir(path):
        return tools.os.run(f"git rev-parse", "git rev-parse", path, None) == 0
    return False


def init_submodules(path, update=True):
    tools.os.run(f"git submodule init", "init-submodule", path)
    if update:
        update_submodules(path)


def update_submodules(path):
    tools.os.run(f"git submodule update --init --recursive", "update-submodule", path)


def checkout(path, tag):
    if is_git_dir(path):
        tools.os.run(f"git reset --hard", "reset", path)
        tools.os.run(f"git clean -fdx", "clean", path)
        tools.os.run(f"git checkout {tag}", "checkout", path)
        # tools.os.run(f"git pull", "pull", path)


def prepare_submodule(path):
    if is_git_dir(path):
        # tools.os.run(f"git reset --hard", "reset", path)
        # tools.os.run(f"git clean -fdx", "clean", path)
        # tools.os.run(f"git pull", "pull", path)

        init_submodules(path)


def submodules(path):
    if is_git_dir(path):
        # Create a ConfigParser instance
        config = configparser.ConfigParser()
        # Preserve case sensitivity in keys
        config.optionxform = str

        # Read the .gitmodules file
        config.read(os.path.join(path, ".gitmodules"))

        # Parse the submodule details
        submodules = []
        for section in config.sections():
            if section.startswith("submodule"):
                submodule_name = section.split('"')[1]  # Extract the submodule name
                d = {key: config[section][key] for key in config[section]}
                d["name"] = submodule_name
                submodules.append(d)

        return submodules

    return []


def bundle(source_path, bundle_dir, name=None):
    source_path = os.path.abspath(source_path)
    bundle_dir = os.path.abspath(bundle_dir)

    if is_git_dir(source_path):
        init_submodules(source_path)
        old_cwd = os.getcwd()
        os.chdir(source_path)
        if name is None:
            name = os.path.basename(source_path)
        uuid_name = uuid.uuid4()
        tools.os.run(f"git bundle create {uuid_name} --all", "bundle", source_path)
        if not os.path.isdir(bundle_dir):
            os.mkdir(bundle_dir)
        os.rename(f"{uuid_name}", os.path.join(bundle_dir, name))
        os.chdir(old_cwd)
        return

    print(f"not a git repo: {source_path}")


def unbundle(bundle_path, outdir):
    bundle_path = os.path.abspath(bundle_path)
    outdir = os.path.abspath(outdir)
    print(f"unbundling {bundle_path} to {outdir}")

    if is_git_dir(outdir):
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        tools.os.run(f"git clone {bundle_path} .", "unbundle", outdir)
        return

    print(f"not a git repo: {bundle_path}")


def bundle_submodules(search_path, bundles_dir, name_prefix=""):
    if is_git_dir(search_path):
        for sm in submodules(search_path):
            sm_name = sm["name"]
            prefix = f"{name_prefix}.{sm_name}"
            source_path = os.path.join(search_path, sm["path"])
            bundle(source_path, bundles_dir, prefix)
            bundle_submodules(source_path, bundles_dir, prefix)
        return
    print(f"not a git repo: {search_path}")


def unbundle_submodules(search_path, bundles_dir, name_prefix=""):
    if is_git_dir(search_path):
        for sm in submodules(search_path):
            sm_path = os.path.join(search_path, sm["path"])
            sm_name = sm["name"]
            prefix = f"{name_prefix}.{sm_name}"
            bundle_path = os.path.join(bundles_dir, prefix)
            if len(os.listdir(sm_path)) == 0:
                unbundle(bundle_path, sm_path)
            unbundle_submodules(sm_path, bundles_dir, prefix)
        return
    print(f"not a git repo: {search_path}")


def unbundle_submodules1(path, bundles):
    if is_git_dir(path):
        for sm in submodules(path):
            submodule_dir = os.path.abspath(os.path.join(path, sm["path"]))
            if len(os.listdir(submodule_dir)) == 0:
                submodule_name = os.path.basename(submodule_dir)
                bundle_path = os.path.abspath(os.path.join(bundles, submodule_name))
                unbundle(bundle_path, os.path.join(submodule_dir, ".."))
            else:
                print(f"not an empty directory: {submodule_dir}")
        return
    print(f"not a git repo: {path}")
