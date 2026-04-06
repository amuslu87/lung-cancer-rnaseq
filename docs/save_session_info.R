# Run this once after completing all analysis to capture package versions
sink("docs/sessionInfo.txt")
sessionInfo()
sink()
cat("Saved: docs/sessionInfo.txt\n")
