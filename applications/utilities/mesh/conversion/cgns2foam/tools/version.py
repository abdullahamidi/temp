#!/usr/bin/env python3

import sys
import semver
import subprocess


def git(*args):
    return subprocess.check_output(["git"] + list(args))


def get_version():
    try:
        latest = (
            git("describe", "--tags").decode().strip().split("-")[0].strip("v") + ".0"
        )
        return {
            "latest": latest,
            "patch": semver.bump_patch(latest),
            "minor": semver.bump_minor(latest),
            "major": semver.bump_major(latest),
        }
    except subprocess.CalledProcessError:
        # No tags in the repository
        return {
            "latest": "1.0.0",
            "patch": "1.0.0",
            "minor": "1.0.0",
            "major": "1.0.0",
        }


def main():
    print(get_version()["patch"])

    return 0


if __name__ == "__main__":
    sys.exit(main())
