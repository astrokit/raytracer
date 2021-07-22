FROM ubuntu:latest

RUN apt update; apt install -y gcc g++ cmake
COPY . /app-data
RUN cd /app-data/build; mkdir cmake-build; cd cmake-build; \
    cmake ../cmake/; make task1

WORKDIR /app-data/data/task1
ENTRYPOINT ["/app-data/build/cmake-build/bin/task1"]
