# QA_REPORT

| check | result | detail |
|---|---|---|
| Q1 prespecification hash unchanged | PASS |  |
| Q2 G1 all ten >=95% (min 99.990) | PASS |  |
| Q3 sentinels 120/120 PASS | PASS |  |
| Q4 G3 and G4 pass all ten | PASS |  |
| Q5 independent recount of tier support equals table | PASS | [(223, 223, 155), (99, 99, 80), (38, 38, 34), (18857, 18701, 1861)] |
| Q6 gene sets 223/99/38 | PASS |  |
| Q7 APOE table equals genes.out and count 7 | PASS |  |
| Q8 no old gene result used as input | PASS |  |
| Q9 required outputs exist | PASS |  |
| Q10 analysis internal checks 19/19 | PASS |  |
| Q11 protected inputs/outputs unchanged (294 files) | PASS |  |
| Q12 nothing outside folder modified | FAIL | <CTAS_WORKSPACE>/CONTINUITY/CONTINUITY_UPDATE_20261104_BERNARD_FINNGEN_AND_FINAL_AUDIT.md;<CTAS_WORKSPACE>/CONTINUITY/04_PEOPLE_TIMELINE_AND_COMMS.md;<CTAS_WORKSPACE>/CONTINUITY/01_PROJECT_STATE.md;<CTAS_WORKSPACE>/CONTINUITY/02_MANUSCRIPT_AND_SCIENCE.md |

11/12 PASS

Note on Q12: the nine files flagged are all in `<CTAS_WORKSPACE>/CONTINUITY/` (modified 11:29-12:10 today). This rerun never wrote there (its code writes only inside this folder and reads CONTINUITY nowhere); they appear to have been updated by another session or by Richard. All 294 protected inputs/outputs hashed in work/PROTECTED_HASHES_pre.tsv are identical after the run (Q11), CONTINUITY was not part of the hash list. Reported, not caused by this run. Final status: 11/12 PASS + 1 externally explained.
