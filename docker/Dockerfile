FROM julia:1.11-bookworm

ENV JULIA_DEPOT_PATH=/opt/julia-depot

ENV DEBIAN_FRONTEND=noninteractive

RUN rm /var/lib/dpkg/info/libc-bin.* \
    && apt-get clean \
    && apt-get update \
    && apt-get install -y \
        git \
        libcurl4-openssl-dev \
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

# Install SAIGE
## I checkout a specific commit for a reproducible state: it is inspired from: https://github.com/saigegit/SAIGE/blob/main/docker/Dockerfile

RUN cd / \
    && git clone https://github.com/saigegit/SAIGE \
    && mv SAIGE app && cd app \
    && git checkout 3add05fc3f80309f8e4fab980e9a410f87cacc8f \
    && curl -fsSL https://pixi.sh/install.sh | sh \
    && mv /root/.pixi/bin/pixi /bin && pixi install \
    && rm -rf /root/.cache \
    && pixi run Rscript -e 'install.packages("lintools", repos="https://cloud.r-project.org")' \
    && curl -L https://github.com/chrchang/plink-ng/archive/refs/tags/v2.0.0-a.6.16.tar.gz | tar -zx \
    && mv plink-ng-2.0.0-a.6.16 plink-ng \
    && pixi run x86_64-conda-linux-gnu-cc -std=c++14 -fPIC -O3 -I.pixi/envs/default/include -L.pixi/envs/default/lib -o plink2_includes.a plink-ng/2.0/include/*.cc -shared -lz -lzstd -lpthread -lm -ldeflate \
    && mv plink2_includes.a .pixi/envs/default/lib \
    && pixi add --platform linux-64 tbb=2020.1 \
    && pixi run R CMD INSTALL . \
    && mv extdata/step1_fitNULLGLMM.R extdata/step2_SPAtests.R extdata/step3_LDmat.R extdata/createSparseGRM.R /usr/local/bin/ \
    && chmod a+x /usr/local/bin/step1_fitNULLGLMM.R /usr/local/bin/step2_SPAtests.R /usr/local/bin/step3_LDmat.R /usr/local/bin/createSparseGRM.R \
    && createSparseGRM.R --help \
    && step1_fitNULLGLMM.R --help \
    && step2_SPAtests.R --help \
    && step3_LDmat.R --help


# Install BCFTOOLS

RUN wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 \
    && tar -xvjf bcftools-1.22.tar.bz2 \
    && cd bcftools-1.22 \
    && ./configure --prefix=/usr/local/ \
    && make \
    && make install

# Install METAL
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

# Install PLINK
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip \
    && unzip plink_linux_x86_64_20250819.zip \
    && mv plink /usr/local/bin \
    && chmod +x /usr/local/bin/plink

# Install PLINK2
RUN wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250129.zip \
    && unzip plink2_linux_x86_64_20250129.zip \
    && mv plink2 /usr/local/bin \
    && chmod +x /usr/local/bin/plink2

# Install REGENIE
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


