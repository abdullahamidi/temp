
# TAIopenFOAMDev

This repository contains all the developments which have been done by the *Computational Fluid Mechanics department
in the Turkish Aerospace* using OpenFOAM as a development platform since 2019.

For now, it mainly contains the Turkish Aerospace Multiphysics Solver **TAMS-AERO** which is dedicated for the 
aerodynamic applications.

![TAIopenFOAMDev](TAMS_AERO_Logo.png)

The main dependencies of **TAMS-AERO** are:

- OpenFOAM (ESI)
- HiSA
- TIOGA
- foam-extend (Ported Code)

## Pre-installation requirements
Before building this project, make sure you have the following dependencies installed:

### Tools
- CMake (>= 3.22.1)
- A C++ compiler (GCC >= 12.3.0 or Clang >= 18.1.3)
### Libraries
- OpenMPI (>= 4.1.2)

## How to install on Ubuntu 22.04/24.04

    git clone https://git.tai.com.tr/t24367/taopenfoamdev.git
    cd taopenfoamdev
    git checkout <tag_name> # (if you need to compile the code at earlier tag, otherwise skip this step)
    git submodule update --init --recursive
    ./Allwmake

After compilation completes, source the tamsAero solver on your terminal:

    source dependencies/OpenFOAM/etc/bashrc

Now, you can check your installation version using `tamsAero -version`

## How to install off-line (e.g., on TAI HPC)

- Download the repository from https://git.tai.com.tr/tams/taopenfoamdev as compressed archive file (zip, tar.gz etc.)
- Transfer the downloaded archive to the HPC and extract its files (e.g., for tar):
<pre>
<code>tar -zxvf taopenfoamdev-master.tar.gz
cd taopenfoamdev-master
./AllcloneBundles
./Allwmake    
</code>
</pre>    
- Source and check tamsAero version as mentioned above

## How to setup development environment
- VSCode extensions
  - Clangd
  - Clang-Format

You need to also install Linux/Windows packages of these tools if they are not already installed on your system. For example, on Ubuntu, you can install Clangd and Clang-Format using the following commands: 

<pre>
<code>sudo apt-get update
sudo apt-get install clangd clang-format
</code>
</pre>

After cloning the repository, install the necessary VSCode extensions for C++ development. These extensions provide features like code completion, formatting, and navigation which enhance the development experience. To install these extensions, open the Extensions view in VSCode (Ctrl+Shift+X), search for "Clangd" and "Clang-Format", and click on the Install button next to each extension. 

Also you need to configure .clangd file in your project root directory with your OpenFOAM source path which is absolute path for `dependencies/OpenFOAM`. Clangd doesn't support environment variables, so you need to provide the full path of OpenFOAM installation directory in .clangd file.