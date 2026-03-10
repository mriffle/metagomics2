# Stage 1: Build frontend
FROM node:20-slim AS frontend-builder

WORKDIR /frontend
COPY frontend/package.json frontend/package-lock.json* ./
RUN npm install
COPY frontend/ ./
RUN npm run build

# Stage 2: Python application
FROM python:3.11-slim

# Version argument - can be overridden at build time
ARG VERSION=0.1.0

LABEL maintainer="Metagomics Team"
LABEL description="Metagomics 2 - Metaproteomics annotation and aggregation tool"
LABEL version="${VERSION}"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install DIAMOND
ARG DIAMOND_VERSION=2.1.21
RUN wget -q https://github.com/bbuchfink/diamond/releases/download/v${DIAMOND_VERSION}/diamond-linux64.tar.gz \
    && tar -xzf diamond-linux64.tar.gz \
    && mv diamond /usr/local/bin/ \
    && rm diamond-linux64.tar.gz

# Create app directory
WORKDIR /app

# Copy project files
COPY pyproject.toml README.md LICENSE ./
COPY src/ ./src/

# Install Python dependencies
RUN pip install --no-cache-dir -e .

# Copy built frontend
COPY --from=frontend-builder /frontend/dist ./frontend/dist

# Download and install reference data
ARG GO_VERSION=2026-01-23

# Download Gene Ontology OBO file
RUN mkdir -p /app/reference/go && \
    wget -q http://purl.obolibrary.org/obo/go/releases/${GO_VERSION}/go.obo \
    -O /app/reference/go/go.obo && \
    echo "${GO_VERSION}" > /app/reference/go/VERSION && \
    echo "source=http://purl.obolibrary.org/obo/go/releases/${GO_VERSION}/go.obo" >> /app/reference/go/VERSION

COPY docker/scripts/fetch_ncbi_taxonomy.sh /usr/local/bin/fetch_ncbi_taxonomy.sh
RUN chmod +x /usr/local/bin/fetch_ncbi_taxonomy.sh

# Download NCBI Taxonomy dump
RUN /usr/local/bin/fetch_ncbi_taxonomy.sh /app/reference/taxonomy

# Create data directory for persistent storage
RUN mkdir -p /data/jobs /data/databases /data/reference

# Copy default config files (users mount their own at /config)
COPY config/databases.example.json /config/databases.json
COPY config/server.example.json /config/server.json

# Copy entrypoint script
COPY docker-entrypoint.sh /app/docker-entrypoint.sh
RUN chmod +x /app/docker-entrypoint.sh

# Set environment variables
ENV METAGOMICS_DATA_DIR=/data
ENV METAGOMICS_CONFIG_DIR=/config
ENV PYTHONUNBUFFERED=1
ENV METAGOMICS_VERSION=${VERSION}
ENV DIAMOND_VERSION=${DIAMOND_VERSION}

# Expose API port
EXPOSE 8000

# Start both worker and server
ENTRYPOINT ["/app/docker-entrypoint.sh"]
