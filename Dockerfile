FROM python:3.12-bookworm

LABEL org.opencontainers.image.source=https://github.com/welch-lab/steamer

RUN <<EOF
set -e
apt update && apt install bedtools libhdf5-dev -y
git clone https://github.com/welch-lab/steamer.git
cd steamer && pip install -r requirements.txt
EOF

RUN pip install ./steamer
