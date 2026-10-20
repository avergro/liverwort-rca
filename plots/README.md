# Figures

Each script draws one figure of the article (`fig01`–`fig06`) or of the Supplementary Material (`figS1`–`figS5`) from the result tables in `data/`, and writes a PDF to `plots/figures/`:

```bash
for f in plots/fig*.py; do python3 "$f"; done
```

The figures produced here use a plain style (`figure_style.py`): matplotlib default fonts and colour cycle, full axes and standard legends. They contain the same data, panels and statistics as the published figures, which were refined afterwards for typography, colour, legend placement and labelling. Differences between the two versions are therefore limited to presentation; every value shown in the article can be traced to a table in `data/` through these scripts.
