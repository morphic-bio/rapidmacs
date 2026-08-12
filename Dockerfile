# RapidMACS container.
#
#   docker build -t rapidmacs .                              # slim, rapidmacs only
#   docker build -t rapidmacs:reproduce --target reproduce . # adds MACS3 3.0.3
#
#   docker run --rm -v "$PWD:/data" rapidmacs \
#     --preset shiftTAG -f FRAG -t /data/fragments.tsv.gz -g hs -n atac \
#     -p 1e-5 --bdgpeakcall-min-len 200 --bdgpeakcall-max-gap 30 \
#     --macs3-frag-narrowpeak /data/atac.narrowPeak \
#     --macs3-frag-summits /data/atac.summits.bed
#
# The image links Debian's htslib rather than the pinned submodule. htslib is
# used for BAM/BEDPE reading only and does not affect peak output.
#
# The default (last) stage is the slim runtime image. `reproduce` is opt-in via
# --target, so a plain `docker build` never produces the large one.

# ---------------------------------------------------------------- builder ---
FROM debian:bookworm-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential \
      libhts-dev \
      zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY Makefile ./
COPY include/ include/
COPY src/ src/
COPY cli/ cli/
RUN make -j"$(nproc)" \
 && strip bin/rapidmacs

# -------------------------------------------------------------- reproduce ---
# Adds released MACS3 3.0.3 and the paper scripts so the published comparison
# can be run inside the container. Opt in with --target reproduce.
FROM debian:bookworm-slim AS reproduce

RUN apt-get update && apt-get install -y --no-install-recommends \
      libhts3 zlib1g libgomp1 \
      python3 python3-venv \
      curl ca-certificates patch time procps util-linux gzip \
 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /src/bin/rapidmacs /usr/local/bin/rapidmacs
COPY LICENSE NOTICE /usr/share/doc/rapidmacs/
COPY paper/ /opt/rapidmacs/paper/
COPY docs/ /opt/rapidmacs/docs/

# Released MACS3 3.0.3 from the SHA-256-pinned sdist, at the path the paper
# scripts default to. Build tooling is installed, used, and removed in one
# layer so it does not persist in the image.
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential python3-dev \
 && DATA_ROOT=/opt/rapidmacs/.macs3-reference VARIANT=released \
      /opt/rapidmacs/paper/scripts/setup_macs3_reference.sh \
 && rm -rf /opt/rapidmacs/.macs3-reference/macs3-reference-src \
 && apt-get purge -y build-essential python3-dev \
 && apt-get autoremove -y \
 && rm -rf /var/lib/apt/lists/*

ENV MACS3=/opt/rapidmacs/.macs3-reference/macs3-3.0.3-venv/bin/macs3 \
    RAPIDMACS=/usr/local/bin/rapidmacs
WORKDIR /opt/rapidmacs
CMD ["/bin/bash"]

# ---------------------------------------------------------------- runtime ---
# Last stage, so this is what a plain `docker build` produces.
FROM debian:bookworm-slim AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
      libhts3 \
      zlib1g \
      libgomp1 \
 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /src/bin/rapidmacs /usr/local/bin/rapidmacs
COPY LICENSE NOTICE /usr/share/doc/rapidmacs/

# Peak calling is chromosome-parallel; give the container the cores you want it
# to use and pass --peak-caller-threads to match.
WORKDIR /data
ENTRYPOINT ["/usr/local/bin/rapidmacs"]
CMD ["--help"]

LABEL org.opencontainers.image.title="RapidMACS" \
      org.opencontainers.image.description="Narrow peak caller with byte-identical output to MACS3 v3.0.3" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.source="https://github.com/morphic-bio/rapidmacs"
