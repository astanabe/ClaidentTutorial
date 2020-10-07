ls RawSequences/*.sha256 | xargs -L 1 -P 8 -I {} sh -c 'sha256sum -c {}'
