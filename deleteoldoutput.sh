git filter-branch --force --index-filter 'git rm -r --cached --ignore-unmatch NonoverlappedPairedEnd_* OverlappedPairedEnd_* PairedEnd_* SingleEnd_* *_RawSequences_*' -- --all
git push --all --force origin
