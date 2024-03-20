FROM python:3.12-bookworm

RUN apt update && apt install bedtools -y

RUN git clone https://github.com/welch-lab/steamer.git

RUN cd steamer && pip install ./
