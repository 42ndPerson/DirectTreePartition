import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import sys
from typing import Optional, Protocol, cast

class _LinregressResult(Protocol):  # minimal protocol for our usage
    slope: float
    intercept: float
    rvalue: float

# --- Helper Utilities for Safe Regression ---
def _filter_positive_finite(n_series: pd.Series, t_series: pd.Series):
    """Return filtered n, T keeping only rows with n>0, T>0 and both finite."""
    mask = (
        (n_series > 0)
        & (t_series > 0)
        & np.isfinite(n_series)
        & np.isfinite(t_series)
    )
    return n_series[mask], t_series[mask]

def _safe_linregress(x, y) -> Optional[_LinregressResult]:
    """Perform linregress with guards.

    Returns None if regression is not meaningful (length < 2 or constant).
    """
    x = np.asarray(x)
    y = np.asarray(y)
    if x.size < 2 or y.size < 2:
        return None
    # Constant check (variance == 0)
    if np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return None
    try:
        reg_any = linregress(x, y)
    except Exception:
        return None
    reg = cast(_LinregressResult, reg_any)
    # Some SciPy versions can yield nan if degenerate; detect & discard
    try:
        values = [float(getattr(reg, attr)) for attr in ("slope", "intercept", "rvalue")]
    except Exception:
        return None
    invalid = any(np.isnan(v) for v in values)
    return None if invalid else reg

def define_complexity_models():
    """
    Returns a dictionary of common complexity functions to test.
    The key is the name (e.g., 'O(n)'), and the value is a
    lambda function that takes 'n' and returns f(n).
    """
    models = {
        'O(1)': lambda n: np.ones_like(np.asarray(n, dtype=float), dtype=float),
        'O(log n)': lambda n: np.log(np.asarray(n, dtype=float)),
        'O(sqrt n)': lambda n: np.sqrt(np.asarray(n, dtype=float)),
        'O(n)': lambda n: np.asarray(n, dtype=float),
        'O(n log n)': lambda n: np.asarray(n, dtype=float) * np.log(np.asarray(n, dtype=float)),
        'O(n^1.5)': lambda n: np.asarray(n, dtype=float)**1.5,
        'O(n^2)': lambda n: np.asarray(n, dtype=float)**2,
        'O(n^2 log n)': lambda n: (np.asarray(n, dtype=float)**2) * np.log(np.asarray(n, dtype=float)),
        'O(n^2.5)': lambda n: np.asarray(n, dtype=float)**2.5,
        'O(n^3)': lambda n: np.asarray(n, dtype=float)**3,
    }
    return models

def load_data(csv_file: str) -> pd.DataFrame:
    """
    Loads the benchmark CSV file and prepares it for analysis.
    
    The CSV is expected to have a MultiIndex header [version, n_size]
    and an index of 'run'.
    """
    # Load the CSV, specifying the two-level header and the index column
    df = pd.read_csv(csv_file, header=[0, 1], index_col=0)
    
    # Get descriptive statistics (mean, std, etc.)
    # We transpose (.T) so that (version, n_size) becomes the index
    stats_df = df.describe().T
    
    # Reset the index to turn (version, n_size) into columns
    stats_df = stats_df.reset_index()
    
    # Ensure n_size is a numeric type for calculations
    stats_df['n_size'] = pd.to_numeric(stats_df['n_size'])
    
    return stats_df

def plot_analysis(version, n, T, log_n, log_T, log_log_reg, models, model_results, best_fit_name, version_color=None):
    """
    Generates two plots for a given algorithm version:
    1. A log-log plot to visualize the O(n^k) fit.
    2. A linear-scale plot showing the data against the best-fitting model.
    """
    plt.figure(figsize=(12, 6))

    # Determine layout based on availability of log-log regression
    has_loglog = log_log_reg is not None
    if has_loglog:
        # --- Plot 1: Log-Log Plot ---
        plt.subplot(1, 2, 1)
        # Scatter plot of the actual data (log-log scale)
        plt.scatter(log_n, log_T, label='Mean Time (Data)', alpha=0.8, color=version_color)
        # Plot the O(n^k) linear regression line
        fit_T_log = log_log_reg.intercept + log_log_reg.slope * log_n
        plt.plot(log_n, fit_T_log, color='red',
                 label=f'Fit: O(n^{log_log_reg.slope:.2f})\n(R² = {log_log_reg.rvalue**2:.4f})')
        plt.xlabel('log(n_size)')
        plt.ylabel('log(mean_time)')
        plt.title(f'Log-Log Complexity for {version}')
        plt.legend()
        plt.grid(True, which="both", ls="--", alpha=0.5)

        # --- Plot 2: Best Fit Model (Linear Scale) ---
        plt.subplot(1, 2, 2)
    else:
        # Single plot if no log-log regression could be computed
        plt.subplot(1, 1, 1)
    
    # Get the data for the best-fitting model
    best_func = models[best_fit_name]
    best_reg = model_results[best_fit_name]['reg_object']
    
    # X-axis is the complexity function, f(n)
    f_n = best_func(n)
    
    # Scatter plot of time vs. f(n)
    plt.scatter(f_n, T, label='Mean Time (Data)', alpha=0.8, color=version_color)
    
    # Plot the linear regression line for the best fit
    if best_reg is not None:
        fit_T_linear = best_reg.intercept + best_reg.slope * f_n
        plt.plot(f_n, fit_T_linear, color='green',
                 label=f'Fit: {best_fit_name}\n(R² = {best_reg.rvalue**2:.4f})')

    plt.xlabel(f'f(n) = {best_fit_name}')
    plt.ylabel('Mean Time (seconds)')
    plt.title(f'Best Fit Model for {version}')
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)
    
    plt.tight_layout()
    plt.suptitle(f'Complexity Analysis for {version}', fontsize=16, y=1.03)
    plt.show()

def analyze_versions(data: pd.DataFrame, models: dict):
    """
    Performs and prints the complexity analysis for each algorithm version
    found in the data.
    """
    versions = data['version'].unique()
    # Prepare color map for versions
    cmap = plt.get_cmap('tab10')
    version_to_color = {v: cmap(i % 10) for i, v in enumerate(versions)}

    # Prepare summary rows
    summary_rows = []
    
    for version in versions:
        print(f"\n--- 📊 Analysis for {version} ---")
        
        # Isolate data for this version and sort by n_size
        df_version = data[data['version'] == version].sort_values('n_size')
        
        n = df_version['n_size']
        T = df_version['mean']  # Mean time

        # --- 0. Clean data (remove non-positive / non-finite entries) ---
        n, T = _filter_positive_finite(n, T)
        if n.empty or T.empty:
            print("⚠️ Skipping: no positive finite data points after cleaning.")
            print("-" * (26 + len(version)))
            continue
        
        # --- 1. Log-Log Analysis (to find k in O(n^k)) ---
        # Ensure all values remain positive for log transform
        log_n = np.log(n)
        log_T = np.log(T)

        log_log_reg = _safe_linregress(log_n, log_T)

        print(f"📈 **Log-Log O(n^k) Fit:**")
        if log_log_reg is None:
            print("   Regression not meaningful (insufficient or constant data).")
        else:
            print(f"   Best fit is **O(n^{log_log_reg.slope:.3f})**")
            print(f"   R-squared: {log_log_reg.rvalue**2:.5f}")

    # --- 2. Specific Model Fitting (R-squared comparison) ---
        model_results = {}
        print(f"\n📋 **Specific Model Fits (Time vs. f(n)):**")
        
        for name, func in models.items():
            # Compute complexity transform f(n)
            try:
                f_n = func(n)
            except Exception:
                f_n = np.array([])
            reg = _safe_linregress(f_n, T)
            if reg is None:
                r2_display = float('nan')
                model_results[name] = {
                    'r_squared': np.nan,
                    'reg_object': None,
                }
            else:
                r2_display = reg.rvalue**2
                model_results[name] = {
                    'r_squared': r2_display,
                    'reg_object': reg,
                }
            print(f"   {name:<12}: R-squared = {r2_display:.5f}")

        # Find and report the best-fitting model
        # Select best fit among models with finite R²
        finite_models = [
            (name, res['r_squared'])
            for name, res in model_results.items()
            if np.isfinite(res['r_squared'])
        ]
        best_fit_name_selected = None
        if finite_models:
            best_fit_name_selected = max(finite_models, key=lambda t: t[1])[0]
            print(
                f"\n   ==> **Best Fit:** **{best_fit_name_selected}** "
                f"(R² = {model_results[best_fit_name_selected]['r_squared']:.5f})"
            )
            print("-" * (26 + len(version)))
            # --- 3. Plotting ---
            plot_analysis(
                version,
                n,
                T,
                log_n,
                log_T,
                log_log_reg,
                models,
                model_results,
                best_fit_name_selected,
                version_color=version_to_color.get(version),
            )
        else:
            print("\n   ==> No meaningful model fit (all regressions invalid).")
            print("-" * (26 + len(version)))

        # --- 4. Collect summary data for this version ---
        row = {
            'version': version,
            'k_exponent': np.nan if log_log_reg is None else float(log_log_reg.slope),
            'k_r2': np.nan if log_log_reg is None else float(log_log_reg.rvalue**2),
        }
        for name in models.keys():
            r2v = model_results.get(name, {}).get('r_squared', np.nan)
            row[f'R2[{name}]'] = r2v
        if finite_models and best_fit_name_selected is not None:
            row['best_model'] = best_fit_name_selected
            row['best_model_r2'] = float(model_results[best_fit_name_selected]['r_squared'])
        else:
            row['best_model'] = None
            row['best_model_r2'] = np.nan
        summary_rows.append(row)

    # --- 5. Print summary table across versions ---
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        # Order columns: version, k_exponent, k_r2, best_model, best_model_r2, then model R2s
        base_cols = ['version', 'k_exponent', 'k_r2', 'best_model', 'best_model_r2']
        model_cols = [c for c in summary_df.columns if c.startswith('R2[')]
        summary_df = summary_df[base_cols + model_cols]
        print("\n====== Summary Across Versions ======")
        print(summary_df.to_string(index=False, float_format=lambda x: f"{x:.4f}" if np.isfinite(x) else "nan"))
        # Optional: save to CSV for later inspection
        try:
            summary_df.to_csv('complexity_summary.csv', index=False)
            print("Saved summary to 'complexity_summary.csv'.")
        except Exception:
            pass

def main():
    """
    Main function to run the analysis.
    """
    csv_file = "benchmarking_mini_results.csv"
    models = define_complexity_models()
    
    try:
        data = load_data(csv_file)
    except FileNotFoundError:
        print(f"Error: '{csv_file}' not found.")
        print("Please run your benchmarking script first to generate the file.")
        sys.exit(1)
    except Exception as e:
        print(f"Error loading or processing '{csv_file}': {e}")
        print("Ensure the file format is correct (MultiIndex header).")
        sys.exit(1)
        
    analyze_versions(data, models)

if __name__ == "__main__":
    main()