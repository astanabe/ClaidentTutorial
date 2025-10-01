git filter-repo --path-glob 'NonoverlappedPairedEnd_*' --path-glob 'OverlappedPairedEnd_*' --path-glob 'PairedEnd_*' --path-glob 'SingleEnd_*' --path-glob '*_RawSequences_*' --invert-paths --force
git remote add origin https://github.com/astanabe/ClaidentTutorial.git
#edit .git/config
git push --all --force origin
