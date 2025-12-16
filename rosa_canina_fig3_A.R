##############################################
#  CIRCULAR GENOME PLOT WITH DUMMY DATA
#  Requires: circlize, ComplexHeatmap
##############################################

# Install packages if needed
if (!requireNamespace("circlize", quietly=TRUE)) install.packages("circlize")
if (!requireNamespace("ComplexHeatmap", quietly=TRUE)) install.packages("ComplexHeatmap")

library(circlize)
# library(ComplexHeatmap)

set.seed(123)

##############################################
# 1. CREATE DUMMY CHROMOSOME SIZES
##############################################

chrom_sizes <- data.frame(
  chr = paste0("Chr", 1:7),
  size = sample(90:150, 7, replace=TRUE) * 1e6   # 90–150 Mb
)

chrom_sizes

##############################################
# 2. INITIALIZE CIRCLIZE WITH CHROMOSOME SECTORS
##############################################

circos.clear()
circos.par(
  "track.height" = 0.08,
  gap.after = c(rep(2, 6), 10)   # Larger gap after last chromosome
)

circos.initialize(
  factors = chrom_sizes$chr,
  xlim = cbind(rep(0, 7), chrom_sizes$size)
)

##############################################
# 3. ADD OUTER TRACKS (S1–S4 HISTOGRAMS)
##############################################

make_dummy_hist <- function(n = 300) {
  abs(rnorm(n)) * runif(1, 0.3, 1.2)
}

S1 <- make_dummy_hist()
S2 <- make_dummy_hist()
S3 <- make_dummy_hist()
S4 <- make_dummy_hist()

cols_outer <- c("purple4", "skyblue", "orchid", "orange1")
outer_data <- list(S1 = S1, S2 = S2, S3 = S3, S4 = S4)

for (i in seq_along(outer_data)) {
  circos.trackPlotRegion(
    ylim = c(0, max(outer_data[[i]])),
    track.height = 0.06,
    bg.border = NA,
    panel.fun = function(x, y) {
      chr = CELL_META$sector.index
      n = length(outer_data[[i]])
      xs = seq(0, chrom_sizes$size[chrom_sizes$chr == chr], length.out=n)
      circos.lines(xs, outer_data[[i]], col = cols_outer[i], lwd = 1)
    }
  )
}

##############################################
# 4. INNER TRACKS (REPEAT FAMILIES)
##############################################

repeat_families <- list(
  CENH3 = make_dummy_hist(),
  CANR4 = make_dummy_hist(),
  ATHILA = make_dummy_hist(),
  Gypsy = make_dummy_hist(),
  Bianca = make_dummy_hist(),
  Copia = make_dummy_hist(),
  DNA_transposon = make_dummy_hist(),
  LINE = make_dummy_hist()
)

cols_inner <- c(
  "red", "magenta", "orange", "purple",
  "blue", "darkgreen", "chartreuse3", "goldenrod"
)

for (i in seq_along(repeat_families)) {
  circos.trackPlotRegion(
    ylim = c(0, max(repeat_families[[i]])),
    track.height = 0.05,
    bg.border = NA,
    panel.fun = function(x, y) {
      chr = CELL_META$sector.index
      n = length(repeat_families[[i]])
      xs = seq(0, chrom_sizes$size[chrom_sizes$chr == chr], length.out=n)
      circos.lines(xs, repeat_families[[i]], col = cols_inner[i], lwd = 1)
    }
  )
}

##############################################
# 5. ADD CHROMOSOME LABELS
##############################################

for (chr in chrom_sizes$chr) {
  circos.text(
    circlize::get.cell.meta.data("xcenter", sector.index=chr, track.index=1),
    circlize::get.cell.meta.data("ycenter", sector.index=chr, track.index=1) + 0.5,
    labels = chr, cex = 0.9
  )
}

##############################################
# 6. ADD LEGEND
##############################################

legend(
  "topright",
  legend = c(names(outer_data), names(repeat_families)),
  col = c(cols_outer, cols_inner),
  lwd = 3, cex = 0.75, bty = "n", title = "Tracks"
)
