# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import pandas as pd

# 1. Load the generated CSV file
df = pd.read_csv("D253_displacement_over_time.csv")

# 2. Set up the plot layout
plt.figure(figsize=(8, 4.5), dpi=300)

# Plot distance over time
plt.plot(
    df["time_ns"],
    df["displacement_A"],
    color="#1f77b4",
    linewidth=1.2,
    alpha=0.85,
    label="D253 C\u00C5 Distance",
)

# Optional: Add a rolling average curve to show long-term structural trends
if len(df) > 50:
    df["rolling_avg"] = df["displacement_A"].rolling(window=50).mean()
    plt.plot(
        df["time_ns"],
        df["rolling_avg"],
        color="#d62728",
        linewidth=2.0,
        label="50-frame Moving Average",
    )

# 3. Labeling and Formatting
plt.title(
    "D253 Structural Displacement Relative to Inactive State (4N6H)",
    fontsize=12,
    pad=12,
    fontweight="bold",
)
plt.xlabel("Time (ns)", fontsize=11)
plt.ylabel("Displacement (\u00C5)", fontsize=11)
plt.xlim(df["time_ns"].min(), df["time_ns"].max())

# Enhance gridlines and layout
plt.grid(True, linestyle="--", alpha=0.5)
plt.legend(frameon=True, loc="upper right")
plt.tight_layout()

# 4. Save figure
plt.savefig("D253_displacement_plot.png", dpi=300)
print("Plot saved successfully as 'D253_displacement_plot.png'")