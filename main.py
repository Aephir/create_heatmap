import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from math import log10 , floor


def round_it(val, type, significant=2):
    try:
        if type == "pEC50":
            rounded_val = round(val, significant - int(floor(log10(abs(val)))))
        else:
            rounded_val = round(val, significant - int(floor(log10(abs(val)))) - 1)
        if 10 <= rounded_val:
            rounded_val = "{:.1f}".format(rounded_val)

        # Use this for 2 decimals when number >10
        # if 10 <= rounded_val < 11:
        #     format_str = "{:." + str(significant) + "f}"
        #     return format_str.format(rounded_val)
        # elif rounded_val >= 11:
        #     return str(int(rounded_val))

        elif rounded_val < 10:
            format_str = "{:." + str(significant) + "f}"
            rounded_val = format_str.format(rounded_val)
        else:
            rounded_val = str(rounded_val)
    except ValueError:  # Catch non-numeric values
        if type == "EC50":
            rounded_val = ">10,000"
        elif type == "pEC50":
            rounded_val = "<5"
        else:  # type == "SEM":
            rounded_val = ""
    return str(rounded_val)

def text_color_for_bg(bg_color):
    if bg_color == '#000000':  # Assuming 'none' is the color passed for NaN cells
        bg_color = "#D3D3D3"

    # Convert the hex color to RGB values
    bg_color = bg_color[1:]  # Remove the '#' at the start
    r, g, b = tuple(int(bg_color[i:i + 2], 16) for i in (0, 2, 4))

    # Calculate brightness
    brightness = (0.299 * r) + (0.587 * g) + (0.114 * b)
    # print(brightness)

    # Return white for dark backgrounds, black otherwise
    return "white" if brightness < 128 else "black"


def read_and_prepare_data(file_path):
    # Read the CSV
    data = pd.read_csv(file_path, sep='\t', header=None)

    # Extract compound names
    compounds = data[0].tolist()

    # Extract data lists
    pEC50_numeric = data.iloc[:, 1::3].replace('NaN', np.nan).astype(float).values
    SEM_list = data.iloc[:, 2::3].replace('NaN', np.nan).astype(float).values
    n_array = data.iloc[:, 3::3].replace('NaN', np.nan).values
    n_list = np.where(np.isnan(n_array), np.nan, n_array.astype(int))

    # Create array for heatmap annotations
    data_info = []
    for pEC50_row, SEM_row, n_row in zip(pEC50_numeric, SEM_list, n_list):
        row_info = []
        for pEC50, SEM, n in zip(pEC50_row, SEM_row, n_row):
            # Convert pEC50 to EC50
            ec50_value = pow(10, -1 * pEC50) * 1000000000
            ec50_str = round_it(ec50_value, "EC50")
            pEC50_str = round_it(pEC50, "pEC50")
            SEM_str = round_it(SEM, "SEM")
            try:
                n = int(n)
            except ValueError:
                pass

            if SEM_str:  # If SEM_str is not empty
                row_info.append(f'{ec50_str}\n{pEC50_str}±{SEM_str} (n={n})')
            else:
                row_info.append(f'{ec50_str}\n{pEC50_str}')
        data_info.append(row_info)

    # Convert arrays to DataFrames
    df_color = pd.DataFrame(pEC50_numeric, columns=["SST1", "SST2", "SST3", "SST4", "SST5"], index=compounds)
    df_info = pd.DataFrame(data_info, columns=["SST1", "SST2", "SST3", "SST4", "SST5"], index=compounds)

    return df_color, df_info

def plot_heatmap(df_color, df_info):
    # Calculate figure height based on number of rows
    n_rows = df_color.shape[0]
    row_height = 16/24
    fig_height = n_rows * row_height

    # Custom gradient colormap
    colors = ["#363636", "#4D91BE", "#F8CD3D"]
    cmap = LinearSegmentedColormap.from_list('custom', colors, N=256)

    sns.set(font="Arial")
    fig, ax = plt.subplots(figsize=(10, fig_height))
    sns.heatmap(df_color, cmap=cmap, cbar=False, ax=ax, linewidths=2, linecolor='black')  # Notice cbar=False

    # Extract position details of the heatmap
    pos = ax.get_position()
    cbar_width = 0.02
    gap = 0.01

    # Define colorbar axes based on main heatmap position
    cbar_ax = fig.add_axes([pos.x1 + gap, pos.y0, cbar_width, pos.height])

    # Create a colorbar
    norm = plt.Normalize(df_color.min().min(), df_color.max().max())
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, cax=cbar_ax,
                        ticks=np.arange(np.floor(df_color.min().min()), np.ceil(df_color.max().max()) + 1))
    cbar.outline.set_edgecolor('black')
    cbar.outline.set_linewidth(2)

    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(2)
        spine.set_edgecolor('black')

    ax.set_xticklabels(df_color.columns, rotation=0, horizontalalignment='center')
    ax.set_yticklabels(df_color.index, rotation=0, verticalalignment='center')

    # Convert df_color values to RGB colors using the cmap
    norm_data = (df_color - df_color.min().min()) / (df_color.max().max() - df_color.min().min())
    colors_rgb = cmap(norm_data.to_numpy().flatten())

    # Convert RGB colors to hex format
    colors_hex = ['#%02x%02x%02x' % (int(r * 255), int(g * 255), int(b * 255)) for r, g, b, _ in colors_rgb]

    color_iter = iter(colors_hex)
    for i in range(df_info.shape[0]):
        for j in range(df_info.shape[1]):
            # Get the background color for the cell from colors_hex
            cell_bg_color = next(color_iter)

            # Determine the appropriate text color based on the cell's background color
            text_color = text_color_for_bg(cell_bg_color)

            # Get the text to display for the cell
            text = df_info.iloc[i, j]

            # Add the text annotation to the cell with the determined text color
            ax.text(j + 0.5, i + 0.5, text, ha='center', va='center', color=text_color)

    return fig, ax


def save_and_show(fig, output_path):
    fig.savefig(output_path, dpi=300, transparent=True, bbox_inches="tight")
    plt.show()


def main():
    files = [
        # "data/fj1_analogs.csv",
        # "data/constomatin_screen.csv",
        # ... [Add more files as required]
        "data/example.csv"
    ]

    for file in files:
        df_color, df_info = read_and_prepare_data(file)
        fig, ax = plot_heatmap(df_color, df_info)
        output_path = file.replace('.csv', '_heatmap.png')
        save_and_show(fig, output_path)


if __name__ == '__main__':
    main()
