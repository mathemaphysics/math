FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -y update
RUN apt-get -y install gcc g++ gfortran gdb cmake
RUN apt-get -y install locate
RUN apt-get -y install git
RUN apt-get -y install libgromacs-dev

ENTRYPOINT [ "/bin/bash" ]
