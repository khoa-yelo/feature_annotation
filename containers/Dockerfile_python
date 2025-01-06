#gete the official Python 3.12 slim image
FROM python:3.12-slim

# Set the working directory in the container
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

# Download and install BLAST+ tools
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-x64-linux.tar.gz \
    && tar -xzvf ncbi-blast-2.11.0+-x64-linux.tar.gz \
    && mv ncbi-blast-2.11.0+ /usr/local/blast \
    && ln -s /usr/local/blast/bin/* /usr/local/bin/

# Install Biopython
RUN pip install --no-cache-dir biopython

# Clone the SeqUtils repository and install it
RUN git clone https://github.com/khoa-yelo/SeqUtils.git

# Change working directory to SeqUtils
WORKDIR /app/SeqUtils

# Install dependencies for SeqUtils
RUN pip install --no-cache-dir -r requirement.txt

# Install SeqUtils as an editable package
RUN pip install --no-cache-dir -e .

# Add custom build directories to PATH
ENV PATH="/oak/stanford/groups/horence/khoa/scratch/build/splash2.6.1:/oak/stanford/groups/horence/khoa/scratch/build:/oak/stanford/groups/horence/khoa/scratch/build/hmmer-3.4/src:${P
ATH}"

# Add the PFAM HMM file path as an environment variable
ENV PFAM_HMM="/oak/stanford/groups/horence/khoa/scratch/data/pfam/Pfam-A.hmm"

# By default, run an interactive Python shell
CMD ["python3"]