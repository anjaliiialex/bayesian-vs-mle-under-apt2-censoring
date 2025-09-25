# plot_dagum_figures.R
library(ggplot2)
library(gridExtra)

# dagum functions (consistent with your implementation)
dagum_pdf <- function(x, a, b, p) {
  (a * p / x) * (x/b)^(a * p) * (1 + (x/b)^a)^(-(p + 1))
}
dagum_cdf <- function(x, a, b, p) {
  (1 + (b/x)^a)^(-p)
}
dagum_surv <- function(x, a, b, p) {
  1 - dagum_cdf(x, a, b, p)
}
dagum_hazard <- function(x, a, b, p) {
  f <- dagum_pdf(x, a, b, p)
  S <- dagum_surv(x, a, b, p)
  f / S
}

x <- seq(1e-3, 10, length.out = 2000)

# Choose representative parameter sets
b <- 1.0
param_sets_pdf <- list(
  c(a=0.8, b=b, p=0.8),
  c(a=1.5, b=b, p=1.2),
  c(a=2.5, b=b, p=1.2)
)

df_pdf <- do.call(rbind, lapply(param_sets_pdf, function(ps) {
  data.frame(x = x, y = dagum_pdf(x, ps['a'], ps['b'], ps['p']),
             a = ps['a'], p = ps['p'])
}))
df_pdf$label <- paste0("a=", df_pdf$a, ", p=", df_pdf$p)

p_pdf <- ggplot(df_pdf, aes(x=x, y=y, color=label)) +
  geom_line(size=1) +
  xlim(0, 10) + ylim(0, NA) +
  labs(title="Dagum PDF for representative parameters", x="x", y="f(x)") +
  theme_minimal()

# Survival plots: vary p
param_sets_surv <- list(
  c(a=1.5, b=b, p=0.5),
  c(a=1.5, b=b, p=1.2),
  c(a=1.5, b=b, p=2.5)
)
df_surv <- do.call(rbind, lapply(param_sets_surv, function(ps) {
  data.frame(x = x, S = dagum_surv(x, ps['a'], ps['b'], ps['p']),
             a = ps['a'], p = ps['p'])
}))
df_surv$label <- paste0("a=", df_surv$a, ", p=", df_surv$p)

p_surv <- ggplot(df_surv, aes(x=x, y=S, color=label)) +
  geom_line(size=1) +
  xlim(0, 10) + ylim(0, 1) +
  labs(title="Dagum Survival Function for varying p (a fixed)", x="x", y="S(x)") +
  theme_minimal()

# Hazard plots: show different hazard shapes
param_sets_hazard <- list(
  c(a=0.8, b=b, p=1.2),
  c(a=1.5, b=b, p=1.2),
  c(a=3.0, b=b, p=1.2)
)
df_hz <- do.call(rbind, lapply(param_sets_hazard, function(ps) {
  hz <- dagum_hazard(x, ps['a'], ps['b'], ps['p'])
  data.frame(x = x, h = hz, a = ps['a'], p = ps['p'])
}))
df_hz$label <- paste0("a=", df_hz$a, ", p=", df_hz$p)

p_hz <- ggplot(df_hz, aes(x=x, y=h, color=label)) +
  geom_line(size=1) +
  xlim(0, 10) + ylim(0, NA) +
  labs(title="Dagum Hazard Function for representative a values", x="x", y="h(x)") +
  theme_minimal()

# Save as PNGs (high resolution)
ggsave("dagum_pdf_param_sets.png", p_pdf, width = 9, height = 5, dpi = 300)
ggsave("dagum_survival_param_sets.png", p_surv, width = 9, height = 5, dpi = 300)
ggsave("dagum_hazard_param_sets.png", p_hz, width = 9, height = 5, dpi = 300)

cat("Plots generated: dagum_pdf_param_sets.png, dagum_survival_param_sets.png, dagum_hazard_param_sets.png\n")
