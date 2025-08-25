FROM debian:bookworm-slim
# Michael J. Foster
# github.com/mjfos2r

# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    perl \
    bash \
    gnuplot \
    && rm -rf /var/lib/apt/lists/*

# Setup all of our directories
RUN mkdir -p /data/ref /alignments /output /opt/mummer
# get mummer4 ...pls work...
RUN wget -P /opt/mummer -O /opt/mummer/mummer.tar.gz https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
RUN tar -xzf /opt/mummer/mummer.tar.gz -C /opt/mummer --strip-components=1 --no-same-owner \
    && cd /opt/mummer \
    && ./configure LDFLAGS=-static \
    && make \
    && make install

# check installations # love this little snippet :)
RUN command -v mummer &>/dev/null || exit 1
RUN command -v nucmer &>/dev/null || exit 1
RUN command -v promer &>/dev/null || exit 1

# Change context to /opt
WORKDIR /opt

COPY [ "entrypoint.sh", "simple_ava.sh", "/opt/" ]

# Set entrypoint
ENTRYPOINT [ "/bin/bash", "entrypoint.sh" ]
