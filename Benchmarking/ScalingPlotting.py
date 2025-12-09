import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os

def plot_results(csv_file="benchmarking_results.csv"):
    """
    Reads the benchmark CSV and plots the results as a line chart
    with confidence intervals. Dynamically adapts to column structure.
    """
    if not os.path.exists(csv_file):
        # Fallback or check for other common names
        if csv_file == "benchmarking_results.csv" and os.path.exists("benchmarking_results_rgebp.csv"):
            print(f"Default file '{csv_file}' not found, switching to 'benchmarking_results_rgebp.csv'")
            csv_file = "benchmarking_results_rgebp.csv"
        else:
            print(f"Error: The file {csv_file} was not found.")
            print("Please run the benchmarking script first to generate the data.")
            return

    try:
        # Read the CSV.
        # header=[0, 1] reads the MultiIndex columns.
        # index_col=0 sets the *first column* (which is the 'run' index) as the index.
        df = pd.read_csv(csv_file, header=[0, 1], index_col=0)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        print("Ensure the CSV file is not empty and has the expected format.")
        return

    # Detect parameter name (e.g., 'n_size' or 'shape')
    # df.columns.names should be ['version', 'parameter_name']
    if len(df.columns.names) > 1 and df.columns.names[1]:
        param_name = str(df.columns.names[1])
    else:
        param_name = 'parameter'

    print(f"Detected parameter: {param_name}")

    # Extract unique parameter values in order to preserve sorting for categorical data
    # (stack() sorts levels, which ruins custom orders like shapes)
    ordered_params = []
    seen = set()
    for _, p in df.columns:
        if p not in seen:
            ordered_params.append(p)
            seen.add(p)

    # --- Data Reshaping ---
    # Stack the multi-level columns into the index
    long_df = df.stack(level=[0, 1])
    
    # Reset the index to turn 'run', 'version', and param_name into columns
    long_df = long_df.reset_index()
    
    # Rename the columns for clarity
    long_df.columns = ['run', 'version', param_name, 'duration']
    
    # Try to convert to numeric if possible
    try:
        long_df[param_name] = pd.to_numeric(long_df[param_name])
        is_numeric = True
    except ValueError:
        is_numeric = False

    # If not numeric, enforce categorical order
    if not is_numeric:
        long_df[param_name] = pd.Categorical(
            long_df[param_name], 
            categories=ordered_params, 
            ordered=True
        )

    print("--- Data prepared for plotting (first 5 rows) ---")
    print(long_df.head())
    print("\nPlotting...")

    # Determine tag from filename
    base_name = os.path.splitext(os.path.basename(csv_file))[0]
    tag = ""
    title_tag = ""
    if base_name.startswith("benchmarking_results_"):
        raw_tag = base_name.replace("benchmarking_results_", "")
        tag = raw_tag + "_"
        title_tag = f" ({raw_tag.upper()})"
    elif base_name != "benchmarking_results":
        tag = base_name + "_"
        title_tag = f" ({base_name.upper()})"

    # --- Plotting ---
    import re
    sns.set_theme(style="whitegrid", font_scale=1.1)
    plt.figure(figsize=(12, 7))

    # Sort versions naturally to ensure legend is in logical order (1, 2, 10 not 1, 10, 2)
    def natural_key(string_):
        return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', str(string_))]
    
    unique_versions = sorted(long_df['version'].unique(), key=natural_key)

    # Create the line plot
    ax = sns.lineplot(
        data=long_df,
        x=param_name,
        y='duration',
        hue='version',
        style='version', # Also vary line style for accessibility
        hue_order=unique_versions,
        style_order=unique_versions,
        markers=True,
        dashes=False,
        estimator='mean',
        errorbar=('pi', 50) 
    )

    # --- Customize Plot ---
    if not is_numeric:
        # Rotate x labels if they are strings (shapes)
        plt.xticks(rotation=45)
    else:
        ax.set_xscale('linear')
    
    ax.set_yscale('linear')

    # Set titles and labels
    formatted_param_name = param_name.replace("_", " ").title()
    if formatted_param_name == "N Size":
        formatted_param_name = "Graph Vertices"

    plt.title(f'Algorithm Performance by {formatted_param_name}{title_tag}', fontsize=16, weight='bold')
    plt.xlabel(formatted_param_name, fontsize=12)
    plt.ylabel('Mean Duration (seconds)', fontsize=12)

    # Improve legend
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1), title='Algorithm Version')

    # Final adjustments
    plt.tight_layout(rect=(0, 0, 0.8, 1)) # Adjust layout to make room for legend
    
    # Save the figure
    output_filename = f"benchmark_{tag}{param_name}_performance_chart.png"
    plt.savefig(output_filename)
    print(f"Chart saved as {output_filename}")

    # Display the plot if interactive
    if sys.stdout.isatty():
        plt.show()


if __name__ == "__main__":
    # Allow passing filename as argument
    if len(sys.argv) > 1:
        plot_results(sys.argv[1])
    else:
        plot_results()