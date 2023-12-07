FROM rachelicr/pisca-branch-master

RUN apt install -y gcc git r-base && \
    pip install phyfum

ENTRYPOINT phyfum
