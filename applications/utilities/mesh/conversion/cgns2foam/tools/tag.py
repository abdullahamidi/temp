#!/usr/bin/env python3

import argparse
import os
import re
import sys
import version


def tag_repo(tag):
    url = os.environ["CI_REPOSITORY_URL"]

    # Transforms the repository URL to the SSH URL
    # Example input: https://gitlab-ci-token:xxxxxxxxxxxxxxxxxxxx@gitlab.com/threedotslabs/ci-examples.git
    # Example output: git@gitlab.com:threedotslabs/ci-examples.git
    push_url = re.sub(r".+@([^/]+)/", r"git@\1:", url)

    version.git("remote", "set-url", "--push", "origin", push_url)
    version.git("tag", tag)
    version.git("push", "origin", tag)


def main():
    parser = argparse.ArgumentParser(
        prog="version", description="Generates new version"
    )

    parser.add_argument("bump-type", type="patch" | "minor" | "major")

    args = parser.parse_args()

    tag_repo(version.get_version()[args.bump_type])

    return 0


if __name__ == "__main__":
    sys.exit(main())
