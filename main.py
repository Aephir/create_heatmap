import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from math import log10, floor

files = [
    # "data/example.csv",
    "data/fj1_analogs.csv"
]


def create_heatmap(file_list, significant_digits=2, pec50_empty="<5", ec50_empty=">10,000") -> None:
    for file in file_list:
        df_color, df_info, pEC50_numeric = read_and_prepare_data(file, significant_digits, pec50_empty, ec50_empty)
        fig, ax = plot_heatmap(df_color, df_info, pEC50_numeric)
        output_path = file.replace(".csv", "_heatmap.png")
        save_and_show(fig, output_path)


def round_to_sf(number, sf):
    """
    Round the given number to specified significant figures (sf).

    Parameters:
        number (float): The number to be rounded.
        sf (int): The number of significant figures.

    Returns:
        float: Rounded number.
    """
    if number == 0:
        return 0
    else:
        order_of_magnitude = floor(log10(abs(number)))
        scale_factor = 10 ** order_of_magnitude
        return round(number / scale_factor, sf - 1) * scale_factor


def round_it(val, type_of_value, significant, pec50_empty, ec50_empty):
    try:
        if type_of_value == "pEC50":
            rounded_val = round(val, significant - int(floor(log10(abs(val)))))
            # Always format with one decimal place
            return f"{rounded_val:.2f}"
        else:  # EC50
            # Adjust the significant digits for EC50
            if 0 < val < 1:
                # The minimum number of decimal places is 2 when val < 1.
                dec_places = max(2, significant - int(floor(log10(abs(val)))) - 1)
            elif 1 <= val < 10:
                # For 1 <= val < 10, always show 1 decimal place
                dec_places = 1
            else:
                # For val >= 10, use round_to_sf
                rounded_val = round_to_sf(val, significant)
                # Ensure the type is float before calling is_integer
                return str(int(rounded_val)) if float(rounded_val).is_integer() else str(rounded_val)

            rounded_val = round(val, dec_places)

        # Formatting
        if type_of_value == "EC50":
            # Ensure that the rounded_val has the correct number of decimal places by formatting it as a string.
            rounded_val_str = f"{rounded_val:.{dec_places}f}"
            # When rounded_val is 10 or greater, do not show decimal places
            if rounded_val >= 10:
                rounded_val_str = f"{int(rounded_val)}"
            return rounded_val_str

        return str(rounded_val)

    except ValueError:  # Catch non-numeric values
        if type_of_value == "EC50":
            return ec50_empty
        elif type_of_value == "pEC50":
            return pec50_empty
        else:  # type_of_value == "SEM":
            return ""


def text_color_for_bg(bg_color):
    if bg_color == "#000000":  # Assuming 'none' is the color passed for NaN cells
        bg_color = "#D3D3D3"

    # Convert the hex color to RGB values
    bg_color = bg_color[1:]  # Remove the '#' at the start
    r, g, b = tuple(int(bg_color[i : i + 2], 16) for i in (0, 2, 4))

    # Calculate brightness
    brightness = (0.299 * r) + (0.587 * g) + (0.114 * b)

    # Return white for dark backgrounds, black otherwise
    return "white" if brightness < 128 else "black"


def read_and_prepare_data(file_path, significant_digits, pec50_empty, ec50_empty):
    # Read the CSV
    data = pd.read_csv(file_path, sep="\t", header=None)

    # Extract compound names
    compounds = data[0].tolist()

    # Extract data lists
    pEC50_numeric = data.iloc[:, 1::3].replace("NaN", np.nan).astype(float).values
    SEM_list = data.iloc[:, 2::3].replace("NaN", np.nan).astype(float).values
    n_array = data.iloc[:, 3::3].replace("NaN", np.nan).values
    n_list = np.where(np.isnan(n_array), np.nan, n_array.astype(int))

    # Create array for heatmap annotations
    data_info = []
    for pEC50_row, SEM_row, n_row in zip(pEC50_numeric, SEM_list, n_list):
        row_info = []
        for pEC50, SEM, n in zip(pEC50_row, SEM_row, n_row):
            # Convert pEC50 to EC50
            ec50_value = pow(10, -1 * pEC50) * 1000000000
            ec50_str = round_it(ec50_value, "EC50", significant_digits, pec50_empty, ec50_empty)
            pEC50_str = round_it(pEC50, "pEC50", significant_digits, pec50_empty, ec50_empty)
            SEM_str = round_it(SEM, "SEM", significant_digits, pec50_empty, ec50_empty)
            try:
                n = int(n)
            except ValueError:
                pass

            if SEM_str:  # If SEM_str is not empty
                row_info.append(f"{ec50_str}\n{pEC50_str}±{SEM_str} (n={n})")
            else:
                row_info.append(f"{ec50_str}\n{pEC50_str}")
        data_info.append(row_info)

    # Convert arrays to DataFrames
    df_color = pd.DataFrame(
        pEC50_numeric, columns=["SST1", "SST2", "SST3", "SST4", "SST5"], index=compounds
    )
    df_info = pd.DataFrame(
        data_info, columns=["SST1", "SST2", "SST3", "SST4", "SST5"], index=compounds
    )

    return df_color, df_info, pEC50_numeric


def plot_heatmap(df_color, df_info, pEC50_numeric):
    # Calculate figure height based on number of rows
    n_rows = df_color.shape[0]
    row_height = 16 / 24
    fig_height = n_rows * row_height

    # Custom gradient colormap
    colors = ["#363636", "#4D91BE", "#F8CD3D"]
    cmap = LinearSegmentedColormap.from_list("custom", colors, N=256)

    sns.set(font="Arial")
    fig, ax = plt.subplots(figsize=(10, fig_height))
    sns.heatmap(
        df_color, cmap=cmap, cbar=False, ax=ax, linewidths=2, linecolor="black"
    )  # Notice cbar=False

    # Extract position details of the heatmap
    pos = ax.get_position()
    cbar_width = 0.02
    gap = 0.01

    # Define colorbar axes based on main heatmap position
    cbar_ax = fig.add_axes([pos.x1 + gap, pos.y0, cbar_width, pos.height])

    # Create a colorbar
    norm = plt.Normalize(df_color.min().min(), df_color.max().max())
    cbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=ax,
        cax=cbar_ax,
        ticks=np.arange(
            np.floor(df_color.min().min()), np.ceil(df_color.max().max()) + 1
        ),
    )
    cbar.outline.set_edgecolor("black")
    cbar.outline.set_linewidth(2)

    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(2)
        spine.set_edgecolor("black")

    ax.set_xticklabels(df_color.columns, rotation=0, horizontalalignment="center")
    ax.set_yticklabels(df_color.index, rotation=0, verticalalignment="center")

    # Convert df_color values to RGB colors using the cmap
    norm_data = (df_color - df_color.min().min()) / (
        df_color.max().max() - df_color.min().min()
    )
    colors_rgb = cmap(norm_data.to_numpy().flatten())

    # Convert RGB colors to hex format
    colors_hex = [
        "#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255))
        for r, g, b, _ in colors_rgb
    ]

    # Get the minimum and maximum pEC50 values used for the heatmap
    data_min = np.nanmin(pEC50_numeric)
    data_max = np.nanmax(pEC50_numeric)

    # Modify the heatmap with custom text and color.
    modify_dict = custom_text(data_min, data_max, cmap)


    # Create a mask based on modify_dict to skip annotations
    mask = np.zeros(df_info.shape, dtype=bool)
    for cell_coords in modify_dict:
        row_idx, col_idx = cell_coords
        mask[row_idx, col_idx] = True

    color_iter = iter(colors_hex)
    for i in range(df_info.shape[0]):
        for j in range(df_info.shape[1]):
            # If the cell is flagged in the mask, skip over it
            if mask[i, j]:
                continue
            # Get the background color for the cell from colors_hex
            cell_bg_color = next(color_iter)

            # Determine the appropriate text color based on the cell's background color
            text_color = text_color_for_bg(cell_bg_color)

            # Get the text to display for the cell
            text = df_info.iloc[i, j]

            # Add the text annotation to the cell with the determined text color
            ax.text(j + 0.5, i + 0.5, text, ha="center", va="center", color=text_color)

    modify_cells(ax, modify_dict)

    return fig, ax


def modify_cells(ax, modify_dict):
    for cell_coords, modification in modify_dict.items():
        row_idx, col_idx = cell_coords  # Unpack the tuple to get row and column indices

        # Set the cell background color
        ax.add_patch(plt.Rectangle((col_idx, row_idx), 1, 1, fc=modification["color"], edgecolor=None))

        # Modify the text if specified
        if "text" in modification:
            ax.text(col_idx + 0.5, row_idx + 0.5, modification["text"],
                    ha="center", va="center", color=text_color_for_bg(modification["color"]))

    # Redraw grid lines
    for edge, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(2)
        spine.set_edgecolor("black")

    # Redraw the grid lines
    for i in range(len(ax.get_xticklabels())):
        ax.axvline(x=i, color='black', linewidth=2)
    for i in range(len(ax.get_yticklabels())):
        ax.axhline(y=i, color='black', linewidth=2)


def custom_text(data_min, data_max, cmap) -> dict[tuple, dict[str, str]]:
    """
    Specify text and color for given cells.
    Cells are referenced by a tuple with Y, X coordinates (first index is row, second is column)
    In the 'modifications' dictionary, use:
        (0, 1): {"color": "#FF0000", "text": "NewText1"}
    to manually add both text and color (in HEX format), or
        (1, 3): {"color": get_color_from_value(4, data_min, data_max, cmap), "text": "NewText2"}
    To give a number value (4 in this case) to automatically color by the color to that number in the heatmap.
    """

    def get_color_from_value(value, vmin, vmax, colormap):
        if isinstance(value, (int, float)):  # Check if the value is a number
            # Normalize the value to the range of the data
            normalized = (value - vmin) / (vmax - vmin)
            normalized = max(min(normalized, 1), 0)  # Clamp between 0 and 1
            # Get the corresponding color from the colormap
            return colormap(normalized)
        else:
            return value  # Return the value unchanged if it's not a number

    # Modify specific cells based on modify_dict
    modifications = {  # Sample dictionary for testing
        # (0, 1): {"color": "#FF0000", "text": "NewText1"},
        # (1, 3): {"color": get_color_from_value(4, data_min, data_max, cmap), "text": "NewText2"}
    }

    # Convert RGB colors to hex format if they are not already in hex
    for cell, modification in modifications.items():
        color = modification["color"]
        if isinstance(color, tuple):  # If the color is an RGB tuple
            modification["color"] = '#{:02x}{:02x}{:02x}'.format(int(color[0] * 255), int(color[1] * 255), int(color[2] * 255))

    return modifications


def save_and_show(fig, output_path):
    fig.savefig(output_path, dpi=300, transparent=True, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    create_heatmap(files, )
