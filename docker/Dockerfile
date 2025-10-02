FROM julia:1.11-bookworm

ENV JULIA_DEPOT_PATH=/opt/julia-depot

ENV DEBIAN_FRONTEND=noninteractive

RUN rm /var/lib/dpkg/info/libc-bin.* \
    && apt-get clean \
    && apt-get update \
    && apt-get install -y \
        git \
        cmake \
        wget \
        procps \
        clang \
        libblas-dev \
        liblapack-dev \
        libopenblas-dev \
        tabix \
        unzip \
        openjdk-17-jdk \
        zlib1g-dev \
        r-base \
        r-base-dev

# METAL
RUN wget https://github.com/statgen/METAL/archive/refs/tags/2020-05-05.zip \
    && unzip 2020-05-05.zip \
    && cd METAL-2020-05-05 \
    && mkdir build && cd build \
    && cmake -DCMAKE_BUILD_TYPE=Release .. \
    && make \
    && make test \
    && make install \
    && cp bin/metal /usr/local/bin/ \
    && chmod +x /usr/local/bin/metal

# Install SuSiE & Rfast
RUN Rscript -e 'install.packages(c("susieR", "Rfast"), repos = "https://cloud.r-project.org/")'

# Install PLINK2
RUN wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250129.zip \
    && unzip plink2_linux_x86_64_20250129.zip \
    && mv plink2 /usr/local/bin \
    && chmod +x /usr/local/bin/plink2

# Install conda
RUN wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" \
    && bash Miniforge3.sh -b -p /opt/miniforge3 \
    && /opt/miniforge3/bin/conda create -y -n regenie_env -c conda-forge -c bioconda regenie==4.1 \
    && /opt/miniforge3/bin/conda clean -a -y \
    && chmod -R 777 /opt/miniforge3/

# Update Path
ENV PATH=/opt/miniforge3/bin/:${PATH}

# Create Julia Sysimage

COPY . /opt/PopGen

RUN export R_HOME=/usr/lib/R \
    && julia --project=/opt/PopGen --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()' \
    && julia --project=/opt/PopGen --startup-file=no --threads=auto -e 'using PackageCompiler; create_sysimage(["PopGen"], sysimage_path="/opt/PopGen/sysimage.so", include_transitive_dependencies=false, cpu_target=PackageCompiler.default_app_cpu_target())'


