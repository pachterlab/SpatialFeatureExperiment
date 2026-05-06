# Internal functions also used in Voyager and Wayfarer

Not meant for the user, but exporting to be used internally in Voyager.
But one day I may clean these up and remove the internal note for people
building on top of SFE.

## Usage

``` r
.value2df(value, use_geometry, feature = NULL)

.check_features(x, features, colGeometryName = NULL, swap_rownames = NULL)

.warn_symbol_duplicate(x, symbols, swap_rownames = "symbol")

.symbol2id(x, features, swap_rownames)

.check_sample_id(x, sample_id, one = TRUE, mustWork = TRUE)

.rm_empty_geometries(g, MARGIN)

.check_rg(type, x, sample_id)

.ext_(x)

.check_tx_file(
  file,
  spatialCoordsNames,
  gene_col,
  phred_col,
  min_phred,
  flip,
  BPPARAM = SerialParam(),
  save_memory = FALSE
)

.check_tx_file_mem(
  file,
  spatialCoordsNames,
  gene_col,
  phred_col,
  min_phred,
  flip,
  BPPARAM = SerialParam()
)
```

## Value

Internal
