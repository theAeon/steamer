FROM python:3.12-bookworm AS pysam-gcs

WORKDIR /src

RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
    | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
    | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list

ENV HTSLIB_CONFIGURE_OPTIONS "--enable-gcs"

RUN apt-get update && apt-get upgrade -y && apt-get install google-cloud-cli -y

RUN pip wheel git+https://github.com/theAeon/pysam.git@gcloud_auth

FROM python:3.12-bookworm AS builder

WORKDIR /src

RUN pip wheel "ncls==0.0.68" "pandas==1.5.3" "sorted-nearest==0.0.39"

FROM python:3.12-bookworm

RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
    | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
    | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list

RUN apt-get update && apt-get upgrade -y && apt-get install bedtools libhdf5-103-1 tabix google-cloud-cli -y

COPY --from=builder /src/pandas-1.5.3-cp312-cp312-linux_x86_64.whl /src/pandas-1.5.3-cp312-cp312-linux_x86_64.whl

COPY --from=builder /src/ncls-0.0.68-cp312-cp312-linux_x86_64.whl /src/ncls-0.0.68-cp312-cp312-linux_x86_64.whl

COPY --from=builder /src/sorted_nearest-0.0.39-cp312-cp312-linux_x86_64.whl /src/sorted_nearest-0.0.39-cp312-cp312-linux_x86_64.whl

COPY --from=pysam-gcs /src/pysam-0.22.1-cp312-cp312-linux_x86_64.whl /src/pysam-0.22.1-cp312-cp312-linux_x86_64.whl

COPY requirements.txt /src/

WORKDIR /src

RUN pip install --no-cache-dir -r requirements.txt --root-user-action ignore

COPY pyproject.toml README.md src /src/

RUN pip install -e ./ --root-user-action ignore --no-deps
