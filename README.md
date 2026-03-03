# TRACLUS

An R package implementing the TRACLUS (TRAjectory CLUStering) algorithm from Lee, Han and Whang (2007). Developed as a course project for *Data Analytics -- Unsupervised Learning* at TU Dresden.

## Installation

The package requires R (>= 4.0) and a working C++ compiler (for the Rcpp components).

```r
# Install from a local clone of the repository
devtools::install("path/to/TRACLUS")

# Or directly from GitHub
devtools::install_github("MartinHoblisch/TRACLUS")
```

## Quick start

```r
library(TRACLUS)

# Run the full pipeline on the built-in toy dataset
result <- traclus(make_toy_trajectories(), eps = 20, min_lns = 3)

# Inspect cluster assignments
table(result$segments$cluster_id, useNA = "ifany")

# Plot the result
plot_traclus_result(result$segments, result$representatives)
```

For a detailed walkthrough of the algorithm and an application to the HURDAT2 Atlantic hurricane dataset, see the package vignette:

```r
vignette("TRACLUS-introduction", package = "TRACLUS")
```

## Algorithm

TRACLUS clusters trajectories in three phases:

1. **Partitioning.** Each trajectory is split into characteristic line segments using the Minimum Description Length (MDL) principle. The greedy scan balances encoding cost of the partition against the deviation cost of the raw points.

2. **Clustering.** The segments are grouped using a density-based algorithm adapted from DBSCAN. Two segments are neighbours if their weighted distance (perpendicular, parallel, angle) is within `eps`. Segments with at least `min_lns` neighbours qualify as core segments. Clusters with fewer than `min_lns` distinct source trajectories are discarded.

3. **Representative trajectories.** For each cluster, a sweep line moves along the average direction. Wherever at least `min_lns` segments are simultaneously active, the average transverse position is recorded as a waypoint.

## Repository structure

```
TRACLUS/
├── R/
│   ├── traclus.R                 # Main pipeline function traclus()
│   ├── traclus_partition.R       # MDL partitioning
│   ├── traclus_distances.R       # Distance functions (perpendicular, parallel, angle)
│   ├── traclus_clustering.R      # Density-based segment clustering
│   ├── traclus_representative.R  # Sweep-line representative trajectories
│   ├── plot_result.R             # plot_traclus_result() visualisation
│   ├── utils_toy.R              # Toy dataset and basic trajectory plotting
│   ├── utils_hurdat.R           # HURDAT2 parser and equirectangular projection
│   └── RcppExports.R            # Auto-generated Rcpp bindings
├── src/
│   ├── traclus_rcpp.cpp          # C++ distance computation and neighbourhood search
│   └── RcppExports.cpp           # Auto-generated Rcpp bindings
├── inst/
│   └── extdata/
│       └── hurdat2-1851-2024-040425.txt  # HURDAT2 best-track data (NOAA NHC)
├── tests/
│   ├── testthat.R                      # Test runner
│   └── testthat/
│       ├── test-traclus.R              # Full pipeline tests
│       ├── test-traclus_distances.R    # Distance function tests
│       ├── test-traclus_partition.R    # MDL partitioning tests
│       ├── test-traclus_clustering.R   # Clustering tests
│       └── test-traclus_representative.R  # Representative trajectory tests
├── man/                          # Auto-generated documentation (.Rd files)
├── vignettes/
│   ├── TRACLUS-introduction.Rmd        # Detailed vignette with toy and HURDAT2 examples
│   └── fig_distance_components.png     # Distance function figure (Lee et al. 2007)
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
├── LICENSE.md
├── .Rbuildignore
├── .gitignore
├── TRACLUS.Rproj
└── README.md
```

## Key functions

| Function | Purpose |
|---|---|
| `traclus()` | Full pipeline: partition, cluster, compute representatives |
| `partition_trajectories()` | MDL-based partitioning of all trajectories |
| `partition_trajectory()` | MDL-based partitioning of a single trajectory |
| `cluster_segments()` | Density-based segment clustering |
| `compute_all_representatives()` | Sweep-line representative trajectories for all clusters |
| `compute_representative()` | Representative for a single cluster |
| `dist_segments()` | Weighted TRACLUS distance between two segments |
| `dist_perpendicular()`, `dist_parallel()`, `dist_angle()` | Individual distance components |
| `estimate_eps()` | Entropy-based heuristic for eps selection |
| `plot_traclus_result()` | Visualise clustering results (two display modes) |
| `plot_trajectories()` | Simple trajectory plot for exploration |
| `read_hurdat2()` | Parse HURDAT2 data and project to Cartesian km |
| `make_toy_trajectories()` | Built-in toy dataset for testing |

## Parameters

- **`eps`**: Distance threshold for neighbourhood queries. Use `estimate_eps()` as a starting point, then experiment.
- **`min_lns`**: Minimum neighbourhood size for core segments. Also the minimum trajectory cardinality per cluster and the density threshold for the sweep-line procedure.
- **`w_perp`, `w_par`, `w_angle`**: Weights for the three distance components (default: all 1).
- **`gamma`**: Minimum spacing between representative trajectory waypoints.

## Performance

The pairwise distance computation in the clustering phase is implemented in C++ via Rcpp. This is the computational bottleneck of the algorithm (quadratic in the number of segments) and the C++ implementation avoids n-squared R-level function call overhead.

## References

Lee, J.-G., Han, J. and Whang, K.-Y. (2007). Trajectory clustering: a partition-and-group framework. *Proceedings of the 2007 ACM SIGMOD International Conference on Management of Data*, 593-604. https://doi.org/10.1145/1247480.1247546

National Hurricane Center (2025). Atlantic hurricane database (HURDAT2). National Oceanic and Atmospheric Administration. https://www.nhc.noaa.gov/data/
