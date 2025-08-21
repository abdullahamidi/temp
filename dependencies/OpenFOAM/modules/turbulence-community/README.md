# OpenFOAM&reg; Turbulence Technical Committee Repository

## Overview

Welcome to the OpenFOAM&reg; [Turbulence Technical Committee](https://wiki.openfoam.com/Turbulence_Technical_Committee) repository - a hub of collaboration and innovation for OpenFOAM&reg; users.

This repository is included with [OpenFOAM](https://www.openfoam.com/)&reg;
as a collection of nested [submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules).

It hosts academic and industrial contributions in the realm of turbulence
modeling to cultivate a space for knowledge exchange, fostering collaboration
that propels the field and OpenFOAM&reg; forward.

Your ideas and insights are not just welcome – they are the building blocks of a
community dedicated to advancing turbulence modeling in unison. Please consider
to join us.

## How to Contribute

### Two-Tier System

The repository follows the so-called the *two-tier system*, employing
the tags `Tier-0` and `Tier-1` in order to effortlessly welcome
community contributions.

The contributions designated as `Tier-0` are only promoted within this
repository, while the `Tier-1` contributions are integrated as a submodule.

The requirements for the *two-tier system* are as follows:

| Tier   | Requirements |
| ---    | ---          |
| Tier-0 | 1. Include a `README.md` file as a user guide.<br>2. Up-to-date contact details of the maintainer. |
| Tier-1 | 1. Tier-0 requirements. <br>2. Error-free compilation with the latest OpenFOAM&reg; version. <br>3. Minimum one reproducible tutorial.|

### Current Contributions

#### Tier-0

| Contribution  | Description |
| ---           | ---          |
| [postChannelFlow](https://github.com/timofeymukha/postChannelFlow) | An improved version of the built-in `postChannel` utility. |
| [runTimeChannelBudgets](https://github.com/janneshopman/runTimeChannelBudgets) | A suite of function objects and utilities to compute the turbulent kinetic energy budget in channel flows. |
| [EllipticBlending](https://github.com/MAHTEP/EllipticBlending) | The k-epsilon Lag Elliptic Blending turbulence model. |
| [SpalartAllmarasRCsend](https://gitlab.com/mAlletto/openfoamtutorials/-/tree/master/SpalartAllmarasRCsend) | Spalart-Allmaras turbulence model with the Spalart-Shur Curvature Correction. |

#### Tier-1

| Contribution  | Description |
| ---           | ---          |
| [libWallModelledLES](https://github.com/timofeymukha/libWallModelledLES/) | A library for wall-modelled large eddy simulation. |

### Add Your Contributions

You can host your code in any Git-based source code repository
hosting service, e.g. Bitbucket or GitHub.

To add your contributions, simply follow the steps below by substituting
`<>`-enclosed texts (or adopting any other suitable method):

1. [Fork this repository, and clone the new repository if necessary](https://docs.gitlab.com/ee/user/project/repository/forking_workflow.html):

```bash
git clone https://gitlab.com/<user-name>/turbulence-community.git turbulence-community
```

2. [Add your repository as a submodule in a feature branch](https://docs.gitlab.com/ee/gitlab-basics/feature_branch_workflow.html):

```bash
cd turbulence-community/
```
```bash
git switch -c <feature-branch-name>
```
```bash
git submodule add <remote-repo-url> <submodule-name>
```
```bash
git commit -m "ENH: <submodule-name>: add a new submodule"
```

3. [Push your feature branch to GitLab](https://docs.gitlab.com/ee/user/project/push_options.html):

```bash
git push --set-upstream origin <feature-branch-name>
```

4. [Open a merge request targeting this repository after your internal tests](https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html).

## How to Use

1. **Source the OpenFOAM&reg; Repository:**

Assuming you have successfully acquired OpenFOAM&reg; by adhering to one of
the [Build Guides](https://develop.openfoam.com/Development/openfoam/-/wikis/building),
the first step is to source the OpenFOAM&reg; environment as follows:

```bash
source <installation path>/etc/bashrc
```

e.g. if installed under the `~/openfoam/OpenFOAM-v2312` directory:

```bash
source ~/openfoam/OpenFOAM-v2312/etc/bashrc
```

2. **Register and Update the Nested Submodules:**

On the first use, you need to register the submodules, and then update them.
You can execute both steps for all the nested submodules as follows while
you are at `$WM_PROJECT_DIR`:

```bash
cd $WM_PROJECT_DIR
```
```bash
git submodule update --init --recursive modules/turbulence-community
```

Executing this single-line command prepares the turbulence community repository
for compilation. Note that other [submodules](https://develop.openfoam.com/Development/openfoam/-/wikis/modules)
available in OpenFOAM&reg; are excluded by explicitly specifying this
submodule's name at the end of the command.

3. **Compile the Nested Submodules:**

Configure and build OpenFOAM&reg; without submodules by
following one of the [Guides](https://develop.openfoam.com/Development/openfoam/-/wikis/home).

If OpenFOAM&reg; has already been built, you can also navigate into the `modules`
directory, and initiate the compilation using all available processors with the
following commands:

```bash
cd $WM_PROJECT_DIR/modules/
```
```bash
./Allwmake -j
```

4. **Remove the Nested Submodules:**

If you need to remove the submodules or wish to restart the process,
you can simply carry out the task as follows:

```bash
cd $WM_PROJECT_DIR
```
```bash
git submodule deinit -f modules/turbulence-community
```

This command deregisters the submodule and clears the
`modules/turbulence-community` directory.

## License

This repository is licensed under the [GNU General Public License v3.0](LICENSE.md).
By contributing to this repository, you agree that your contributions will be
licensed under the same license.

## Code of Conduct

We encourage a positive and collaborative environment.
Please adhere to our [Code of Conduct](CODE_OF_CONDUCT.md)
when participating in this repository.

## Contact

[Please create an issue](https://gitlab.com/openfoam/community/tc-turbulence/turbulence-community/-/issues/new)
if you have any suggestions, questions or contributions. Please be aware that the
activities of the Turbulence Technical Committee are voluntary; consequently, delays
in responses should be expected.

If you would like to join the Turbulence Technical Committee,
[please get in touch](https://wiki.openfoam.com/Turbulence_Technical_Committee).

## Acknowledgments

We appreciate the efforts of the entire habitat of OpenFOAM&reg; users in
[contributing](CONTRIBUTORS.md) to the improvement and evolution of this
powerful CFD tool.

Thank you for your interest and contribution to the OpenFOAM&reg; community!


<!----------------------------------------------------------------------------->