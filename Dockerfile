FROM continuumio/miniconda3:latest

LABEL maintainer="Brown Beckley <brownbeckley94@gmail.com>"
LABEL description="Kleboscope - Comprehensive Klebsiella pneumoniae Genomic Typing Platform"

# Install system dependencies (procps for resource monitoring, jq for JSON parsing)
RUN apt-get update && apt-get install -y \
    procps \
    jq \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt/kleboscope

# Copy entire project
COPY . /opt/kleboscope/

# Create the Conda environment from environment.yml
RUN conda env create -f environment.yml && \
    conda clean -afy

# Make the environment the default for RUN commands
SHELL ["conda", "run", "-n", "kleboscope", "/bin/bash", "-c"]

# Run abricate database setup (one-time)
RUN abricate --setupdb

FROM continuumio/miniconda3:latest

# Run abricate database setup (one-time)
RUN abricate --setupdb

# Run AMR database setup (downloads latest AMRfinderPlus database)
RUN cd /opt/kleboscope/kleboscope/modules/kleb_amr_module && \
     python klebo_amrfinder.py --update-db   

# Set entrypoint to kleboscope command
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "kleboscope", "kleboscope"]
CMD ["-h"]
