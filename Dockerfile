FROM gcc:latest
RUN apt update && apt upgrade -y && apt install -y cmake libc6-dbg gdb valgrind doxygen groff xxd libcairo2-dev libeigen3-dev