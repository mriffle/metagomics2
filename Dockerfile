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
ARG NCBI_TAXONOMY_DATE=2026-02-97

# Download Gene Ontology OBO file
RUN mkdir -p /app/reference/go && \
    wget -q http://purl.obolibrary.org/obo/go/releases/${GO_VERSION}/go.obo \
    -O /app/reference/go/go.obo && \
    echo "${GO_VERSION}" > /app/reference/go/VERSION && \
    echo "source=http://purl.obolibrary.org/obo/go/releases/${GO_VERSION}/go.obo" >> /app/reference/go/VERSION

# Download NCBI Taxonomy dump
RUN mkdir -p /app/reference/taxonomy && \
    wget -q https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz && \
    tar -xzf taxdump.tar.gz -C /app/reference/taxonomy && \
    rm taxdump.tar.gz && \
    echo "${NCBI_TAXONOMY_DATE}" > /app/reference/taxonomy/VERSION && \
    echo "source=https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" >> /app/reference/taxonomy/VERSION

# Create data directory for persistent storage
RUN mkdir -p /data/jobs /data/databases /data/reference

# Copy entrypoint script
COPY docker-entrypoint.sh /app/docker-entrypoint.sh
RUN chmod +x /app/docker-entrypoint.sh

# Set environment variables
ENV METAGOMICS_DATA_DIR=/data
ENV PYTHONUNBUFFERED=1
ENV METAGOMICS_VERSION=${VERSION}
ENV DIAMOND_VERSION=${DIAMOND_VERSION}

# Expose API port
EXPOSE 8000

# Start both worker and server
ENTRYPOINT ["/app/docker-entrypoint.sh"]
