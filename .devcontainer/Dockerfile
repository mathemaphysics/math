FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update -y
RUN apt-get install -y gcc g++ gdb cmake
RUN apt-get install -y python3 python3-opencv
RUN apt-get install -y locate
RUN apt-get install -y git

ENTRYPOINT [ "/bin/bash" ]
