#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

infile <- "nonmonophyly_affecting_basal_w5000.any_ingroup_nonmono_summary.tsv"
out_pdf <- "nonmonophyly_affecting_basal_w5000.any_ingroup_nonmono_basal_status.pdf"
out_png <- "nonmonophyly_affecting_basal_w5000.any_ingroup_nonmono_basal_status.png"

summary <- read.delim(infile, check.names = FALSE)

plot_data <- data.frame(
  Region = rep(summary$Region, each = 3),
  Category = rep(
    c("South-Asia-only basal", "Non-South-Asian basal", "Ambiguous"),
    times = nrow(summary)
  ),
  Windows = c(
    rbind(
      summary$South_Asia_only_basal_windows,
      summary$Non_South_Asian_basal_windows,
      summary$Ambiguous_windows
    )
  ),
  Total = rep(summary$At_least_one_ingroup_strict_nonmono_windows, each = 3)
)

plot_data$Fraction <- plot_data$Windows / plot_data$Total
plot_data$Label <- sprintf("%.1f%%", 100 * plot_data$Fraction)
plot_data$Region <- factor(plot_data$Region, levels = c("Autosome", "X"))
plot_data$Category <- factor(
  plot_data$Category,
  levels = c("Ambiguous", "Non-South-Asian basal", "South-Asia-only basal")
)

p <- ggplot(plot_data, aes(x = Region, y = Fraction, fill = Category)) +
  geom_col(width = 0.56, color = "white", linewidth = 0.35) +
  geom_text(
    aes(label = Label),
    position = position_stack(vjust = 0.5),
    size = 2.7,
    color = "black"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "South-Asia-only basal" = "#2B7A78",
      "Non-South-Asian basal" = "#C65D32",
      "Ambiguous" = "#8A8F98"
    )
  ) +
  labs(
    x = NULL,
    y = "Proportion of non-monophyletic local-tree windows",
    fill = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 8),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(margin = margin(r = 8)),
    plot.margin = margin(8, 10, 8, 8)
  )

ggsave(out_pdf, p, width = 4.6, height = 3.4, units = "in", device = cairo_pdf)
ggsave(out_png, p, width = 4.6, height = 3.4, units = "in", dpi = 600)

message("[DONE] ", out_pdf)
message("[DONE] ", out_png)
