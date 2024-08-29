FROM python:3.12-bookworm

LABEL org.opencontainers.image.source=https://github.com/welch-lab/steamer

RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
    | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
    | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list

RUN apt-get update && apt-get upgrade -y && apt-get install bedtools libhdf5-dev tabix google-cloud-cli -y

COPY requirements.txt /src/

WORKDIR /src/

RUN pip install --no-cache-dir -r requirements.txt

COPY pyproject.toml README.md src /src/

RUN pip install ./
