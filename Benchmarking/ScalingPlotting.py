import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_results(csv_file="benchmarking_results.csv"):
    """
    Reads the benchmark CSV and plots the results as a line chart
    with confidence intervals.
    """
    try:
        # Read the CSV.
        # header=[0, 1] reads the MultiIndex columns.
        # index_col=0 sets the *first column* (which is the 'run' index) as the index.
        df = pd.read_csv(csv_file, header=[0, 1], index_col=0)
    except FileNotFoundError:
        print(f"Error: The file {csv_file} was not found.")
        print("Please run the benchmarking script first to generate the data.")
        return
    except Exception as e:
        print(f"Error reading CSV: {e}")
        print("Ensure the CSV file is not empty and has the expected format.")
        return

    # --- Data Reshaping ---
    # The data is "wide": columns are ('version', 'n_size')
    # We need to "melt" or "stack" it into a "long" format
    # for plotting with seaborn.
    #
    # Stack the multi-level columns into the index
    long_df = df.stack(level=[0, 1])
    
    # Reset the index to turn 'run', 'version', and 'n_size' into columns
    long_df = long_df.reset_index()
    
    # Rename the columns for clarity
    long_df.columns = ['run', 'version', 'n_size', 'duration']
    
    # 'n_size' was part of the column name, so it's a string
    # Convert it to a numeric type for plotting
    long_df['n_size'] = pd.to_numeric(long_df['n_size'])

    print("--- Data prepared for plotting (first 5 rows) ---")
    print(long_df.head())
    print("\nPlotting...")

    # --- Plotting ---
    sns.set_theme(style="whitegrid", font_scale=1.1)
    plt.figure(figsize=(12, 7))

    # Create the line plot
    # - x='n_size': n_size on the x-axis
    # - y='duration': duration on the y-axis
    # - hue='version': Create a separate line for each algorithm version
    # - marker='o': Add a marker at each data point
    # - errorbar=('pi', 50): This is the key part!
    #   'pi' = percentile interval
    #   50 = show the 50% interval centered on the median (i.e., 25th to 75th percentile)
    # - estimator='mean': The line itself will be the mean of all runs
    ax = sns.lineplot(
        data=long_df,
        x='n_size',
        y='duration',
        hue='version',
        style='version', # Also vary line style for accessibility
        markers=True,
        dashes=False,
        estimator='mean',
        errorbar=('pi', 50) 
    )

    # --- Customize Plot ---
    # Set log scales
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    
    # Linear optionally
    ax.set_xscale('linear')
    ax.set_yscale('linear')

    # Set titles and labels
    plt.title('Algorithm Performance by Graph Size', fontsize=16, weight='bold')
    plt.xlabel('Graph Size (n_size)', fontsize=12)
    plt.ylabel('Mean Duration (seconds)', fontsize=12)

    # Improve legend
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1), title='Algorithm Version')

    # Final adjustments
    plt.tight_layout(rect=(0, 0, 0.8, 1)) # Adjust layout to make room for legend
    
    # Save the figure
    output_filename = "benchmark_performance_chart.png"
    plt.savefig(output_filename)
    print(f"Chart saved as {output_filename}")

    # Display the plot
    plt.show()


if __name__ == "__main__":
    plot_results()