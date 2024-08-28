FROM python:3.12-bookworm

LABEL org.opencontainers.image.source=https://github.com/welch-lab/steamer

RUN apt-get update && apt-get install bedtools libhdf5-dev tabix -y

COPY requirements.txt /src/

WORKDIR /src/

RUN pip install --no-cache-dir -r requirements.txt

COPY pyproject.toml README.md src /src/

RUN pip install ./
