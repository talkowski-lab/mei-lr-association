## Environment for assembly_coverage.wdl: samtools (depth).
FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends samtools \
    && rm -rf /var/lib/apt/lists/*
COPY scripts/ ./scripts/
CMD ["bash"]
