# Start with NVIDIA's CUDA Ubuntu base image
FROM continuumio/miniconda3
RUN conda update -n base -c defaults conda

# Install Python 3.8
RUN conda install python=3.8

RUN conda install -c conda-forge -c salilab -y git openmm cudatoolkit=10.0 mdtraj gromacs boost=1.73 dssp
RUN git clone https://github.com/maccallumlab/martini_openmm.git \
    && cd martini_openmm \
    && python setup.py install \
    && cd ..
RUN pip install vermouth==0.9.3
# Install build tools for DSSP, libcif++, and libmcfp
#RUN apt-get update && apt-get install -y build-essential cmake

## Install zlib, eigen3, and boost development files
#RUN apt-get install -y zlib1g-dev libeigen3-dev libboost-all-dev
#
## Clone, build, and install libcif++
#RUN git clone https://github.com/PDB-REDO/libcifpp.git --recurse-submodules libcifpp
#WORKDIR libcifpp
#RUN cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_BUILD_TYPE=Release
#RUN cmake --build build
#RUN cmake --install build
#WORKDIR /
#
## Clone, build, and install libmcfp
#RUN git clone https://github.com/mhekkel/libmcfp.git
#WORKDIR libmcfp
#RUN mkdir build
#WORKDIR build
#RUN cmake ..
#RUN cmake --build .
#RUN cmake --install .
#WORKDIR /
#
## Clone and build DSSP
#RUN git clone https://github.com/PDB-REDO/dssp.git
#WORKDIR dssp
#RUN mkdir build
#RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#RUN cmake --build build
#RUN cmake --install build
#WORKDIR /
#
## Clean up
#RUN apt-get autoremove -y && apt-get clean
