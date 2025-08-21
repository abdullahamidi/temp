import os
import shutil

SHELL_EXT = ".bat" if os.name == "nt" else ".sh"
EXE_EXT = ".exe" if os.name == "nt" else ""
LIST_SEP = ";" if os.name == "nt" else ":"
PATH_SEP = "\\" if os.name == "nt" else "/"


def run(cmd, description, workdir=None, expected=0):
    owd = os.getcwd()
    if workdir is not None:
        os.chdir(workdir)
    print(f"[execute] {cmd} - {os.getcwd()}")
    
    retval = os.system(cmd)
    os.chdir(owd)
    if expected is not None and retval != expected:
        print(f"{description} failed")
        exit(-1)

    if expected is None:
        return retval
    
    return True

def check(cmd, fail_msg):
    if shutil.which("vcvars64.bat") is None:
        print(fail_msg)
        exit(-1)

    return True


def add_env(key, value):
    old = os.environ.get(key)
    if old is not None:
        os.environ[key] = f"{value}{LIST_SEP}{old}"
    else:
        os.environ[key] = value
