FROM gcc:latest
RUN apt update && apt upgrade -y && apt install -y groff xxd libcairo2-dev libeigen3-dev