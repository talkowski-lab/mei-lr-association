## Environment for extract_insertion_methylation.py
## pysam ships manylinux wheels, so no build toolchain is needed.
FROM python:3.11-slim
RUN pip install --no-cache-dir pysam pandas
COPY scripts/ ./scripts/
CMD ["bash"]
