FROM gitpod/workspace-full-vnc

RUN sudo apt-get update && sudo sudo apt install -y \
  xorg-dev \
  libglu1-mesa-dev \
  python3-dev \
  libsdl2-dev \
  libc++-7-dev \
  libc++abi-7-dev \
  ninja-build \
  libxi-dev \
  libtbb-dev \
  libosmesa6-dev \
  libudev-dev \
  autoconf \
  libtool
