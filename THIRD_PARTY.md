# Third-party code and adaptations

## Circle_finder

- Repository: https://github.com/pk7zuva/Circle_finder
- Commit: `3eb333db2ea6277dde36cbf640be9afeb710c717`
- File: `circle_finder-pipeline-bwa-mem-samblaster.sh`
- Original source and MIT license are retained in `third_party/Circle_finder/`.
- `cf_lines()` extracts the original step 1–2 commands for mapping and step 3–11 commands for calling. The original final wildcard `rm` is never executed.
- Each original command is run with Bash strict error handling and pipefail. `grep` status 1 (no matches) is normalized to success; status >1 remains an error.
- Unique sample working directories and fixed local argument names avoid cross-sample collisions. Scientific AWK/filter expressions are not rewritten.
- These safety changes and current dependency versions still require output comparison on the cluster before large-scale use.

## CReSIL

- Repository: https://github.com/visanuwan/cresil
- Commit: `78356affa6adca6a41ea297c3d7d93f5ecb5ae6a`
- `envs/cresil-upstream.yml` is copied from that revision. Installation fetches upstream with its own MIT license.
- preflight verifies Git blob hashes of installed `trim.py` and `identify_wgls.py` against the pinned revision.
- `scripts/cresil_wgls_checked.py` invokes the installed WGLS module without changing graph/detection logic. It replaces unchecked `subprocess.call` with checked execution, using Bash pipefail for shell pipelines, and normalizes three explicitly enumerated no-candidate exits to a documented empty result.
- Zero trimmed regions, unknown exit messages, subprocess errors and exceptions remain failures; they are not silent biological negatives.
- For the WGLS FASTA indexing used by upstream, a FASTA copy is created while full FASTQ qualities remain available for evidence extraction.
- CReSIL reported regions can include acyclic graphs; `ctc=True` is used for candidate counts, not as physical-circle proof.

All third-party authors retain their copyrights. Cite the original methods when using results.
