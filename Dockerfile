FROM python:3.12-bookworm

LABEL org.opencontainers.image.source=https://github.com/welch-lab/steamer

RUN apt update && apt install bedtools libhdf5-dev -y

RUN git clone https://github.com/welch-lab/steamer.git

RUN cd steamer && pip install ./
