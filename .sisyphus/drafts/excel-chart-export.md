# Draft: Excel Export Chart Sheets

## Requirements (confirmed)
- Add per-gene "{Gene}_Chart" sheets to Excel export
- Each chart sheet: data table (col C, row 2) + embedded openpyxl BarChart (col K)
- Y-axis: rich text with RED gene name: "Relative mRNA expression level of GENE/β-Actin"
- SEM error bars (plus-only, custom values)
- White-filled bars, black outline, no legend, gap width 219, overlap -27
- Chart size: width=15, height=7.5
- HK gene name needed for Y-axis → add param to export_to_excel

## Technical Decisions
- Engine approach: User prefers hybrid (xlsxwriter → reopen with openpyxl) as "cleaner"
- openpyxl already in requirements.txt ✓
- xlsxwriter only used in export_to_excel ✓

## Research Findings
- processed_data: Dict[gene: DataFrame], ONE ROW PER CONDITION (not per replicate)
- Fold_Change and SEM are already per-condition aggregates → no aggregation needed for chart
- Fold_Change == Relative_Expression (duplicated column)
- export_to_excel has ZERO test coverage
- 81 tests total: 78 pass, 3 fail (PPT-related, unrelated)
- No tests import/mock xlsxwriter → engine switch safe

## Call Sites (both identical pattern)
1. Line 4719: direct download button
2. Line ~5055: ZIP bundle writestr()
- Both need hk_gene parameter added

## Open Questions
- Hybrid approach vs full openpyxl switch → which to recommend?
- Test strategy for new chart code?
- Condition ordering in chart bars?
- HK gene in Y-axis: use actual hk_gene name or hardcode "β-Actin"?
- COMPLETED: librarian openpyxl research (see below)

## openpyxl API Research (from librarian)

### Bar Chart Styling
- `series.graphicalProperties.solidFill = "FFFFFF"` for white fill
- `series.graphicalProperties.line.solidFill = "000000"` for black outline
- `chart.type = "col"`, `chart.grouping = "clustered"`

### Error Bars (PLUS-ONLY with custom SEM)
- Use `ErrorBars(plus=nds, errDir='y', errBarType="plus", errValType="cust")`
- Build custom values: `NumVal(i, None, v=x)` → `NumData(pt=numvals)` → `NumDataSource(numLit=nd)`
- Reference: https://github.com/uskysd/openpyxl-errorbar/blob/master/errorbar.py

### Rich Text Y-axis Title
- Use `openpyxl.chart.text.RichText` (NOT cell.rich_text)
- `Paragraph(r=[RegularTextRun(rPr=CharacterProperties(sz=800, solidFill="000000"), t="text")])`
- Font size in 1/100 points: 800 = 8pt
- Mix runs with different solidFill colors in same Paragraph

### Hybrid Approach Gotchas
- File locking: xlsxwriter must be fully closed before openpyxl load
- Performance: two serialization passes
- With BytesIO: must seek(0) before load_workbook, and load_workbook needs a file-like object

### Chart Anchoring
- Simple: `ws.add_chart(chart, "K2")` anchors top-left to K2
- Chart size: `chart.width = 15`, `chart.height = 7.5` (in cm)

### Gap Width / Overlap
- `chart.style` for preset styles
- `series.graphicalProperties` for individual bar styling
- Gap width: likely need XML-level attribute or chart.gapWidth property

## Scope Boundaries
- INCLUDE: export_to_excel function, call sites, chart generation
- EXCLUDE: existing sheets, GraphGenerator class, Plotly graphs, test fixes
