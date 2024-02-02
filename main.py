"""
Author:      @Aephir (Walden Bjørn-Yoshimoto)
Repository:  https://github.com/Aephir/create_heatmap
Version:     2.0.6

Generates heatmap images from a csv file with data in the format:
    Any number of columns with 3 sub-columns in each (plus one column for row names).
    The 3 sub-columns should have:
        1) pEC50 value as a float or integer.
        2) Error range (± CI95, SEM, etc.) as a float or integer.
        3) Number of independent repeats (n) as an integer.
    Rows as needed (e.g., different compounds or conditions) with row names in the first column. No column headers.
Parameters that can be adjusted are in the "adjustable_parameters" dictionary.
Current options are:
    1) "files": A list of the files cvs you want to create heatmaps from.
    2) "significant_digits": The number of significant digits to use (NOT significant decimals)
        - Caveat: For pEC50 values, "significant_digits" is number of decimals, since this is a log scale.
    3) "column_names": names of column headers
    4) "pec50_empty" the text to use for pEC50 if no value in cell
    5) "ec50_empty" the text to use for EC50 if no value in cell
    6) "modifications": Allow custom color and text for specified cells. Use case could be:
        - If you cannot dissolve a compound enough to use the standard "pec50_empty" and "ec50_empty".
        - If you can determine approximate values, and want to indicate this.
        The tuple (e.g., (0, 1)) determines the location in the heatmap to change (index, column).
        Remember to start counting at 0. For example:
            (0, 1) means index=0 (or row=0), column=1, so first rw, second column.
            (23, 4) means index=23, column=4, so 24th row, fifth column.
        You supply a numeric value that will act as if that value was in your dataset.
    7) "font_settings": settings for the font type, size, and weight (e.g., "normal", "bold") to use
    8) "show_plot": Show the plot after the script is done (as opposed to just saving the file)
    9) "save_as_type": Which types of files should be saved. Current options are "png" and "pdf".
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from math import log10, floor
from typing import Any, Union
from numpy import ndarray, dtype
from pandas import DataFrame
from const import (
    SUBSTITUTE_SCREEN,
    SUBSTITUTE_ANALOGS
)

datasets = {
    'screen': ['data/consomatin_screen_ci95.csv', SUBSTITUTE_SCREEN],
    'analogs': ['data/fj1_analogs_ci95.csv', SUBSTITUTE_ANALOGS]
}
user_input = input('Which dataset do you want to use?\n\n[1] Screen\n[2] Analogs\n\n')
if int(user_input) == 1:
    dataset = 'screen'
elif int(user_input) == 2:
    dataset = 'analogs'
else:
    print(f'"{user_input}" is not a valid selection. Terminating.')
    quit()

# Adjustable parameters
parameters = {
    "files": [
        datasets[dataset][0]
    ],
    "significant_digits": 2,
    "column_names": ["SST1", "SST2", "SST3", "SST4", "SST5"],
    "pec50_empty": "<5",
    "ec50_empty": ">10,000",
    "modifications": datasets[dataset][1],
        # {  # Uncomment line(s) below and add values if/as needed
        # # (0, 1): {"color": "#FF0000", "text": "Super\nModified"},  # This version is not currently working!
        # (1, 1): {"color": 6, "text": ">1,000\n<6"},  # Fj1 @ SST1
        # (23, 0): {"color": 5.3, "text": ">1,000\n<6"},  # TT-232 @ SST1
        # (23, 3): {"color": 5.8, "text": ">1,000\n<6"},  # TT-232 @ SST4
        # },
    "font_settings": {'font': 'Arial', 'size': 12, 'weight': 'normal'},
    "show_plot": True,
    "save_as_type": ["pdf", "png"]  # available options are currently "pdf" and "png".
}


class HeatmapGenerator:
    def __init__(self, adjustable_parameters: dict[Any]):
        self.files: list[str] = adjustable_parameters["files"]
        self.significant_digits: int = adjustable_parameters["significant_digits"]
        self.pec50_empty: str = adjustable_parameters["pec50_empty"]
        self.ec50_empty: str = adjustable_parameters["ec50_empty"]
        self.cmap: LinearSegmentedColormap = self.create_colormap()
        self.font_settings: dict[str, Any] = adjustable_parameters["font_settings"]
        self.modifications = adjustable_parameters["modifications"]
        self.rows: int = 24  # placeholder - will be updated in ``read_and_prepare_data
        self.columns: int = 5  # placeholder - will be updated in ``read_and_prepare_data
        self.column_names: list[str] = adjustable_parameters["column_names"]
        self.show = adjustable_parameters["show_plot"]
        self.save_as_type = adjustable_parameters["save_as_type"]

    def create_heatmaps(self) -> None:
        """Main method to create heatmap"""
        for file in self.files:
            df_color, df_info, pec50_numeric = self.read_and_prepare_data(file)
            df_color_modified, df_info_modified = self.replace_values(df_color, df_info)
            fig, ax = self.plot_heatmap(df_color_modified, df_info_modified, pec50_numeric)
            output_path = file.replace(".csv", "_heatmap")
            self.save_and_show(fig, output_path)

    def round_to_sf(self, number: Union[int, float]) -> float:
        """
        Rounds the number to show the number of significant digits.
        Importantly, NOT the number of significant decimals, but significant digits.
        """
        if number == 0:
            return 0
        else:
            return round(number, self.significant_digits - int(floor(log10(abs(number)))) - 1)

    def round_it(self, val, type_of_value):
        """

        :param val:
        :param type_of_value:
        :return:
        """
        try:
            if type_of_value == "pEC50":
                rounded_val = round(val, self.significant_digits - int(floor(log10(abs(val)))))
                # Always format with two decimal place
                return f"{rounded_val:.2f}"
            else:  # EC50
                # Here we only change EC50 rounding to tackle floating point errors
                # Use round_to_sf to get a preliminary rounded value with the correct number of significant figures
                rounded_val = self.round_to_sf(val)

                # Now determine how many decimal places should be shown
                if val < 1:
                    # Values less than 1 should have at least 2 decimal places
                    dec_places = max(2, self.significant_digits - int(floor(log10(abs(val)))) - 1)
                elif val < 10:
                    # Values between 1 and 10 should have at least 1 decimal place
                    dec_places = 1
                else:
                    # Values greater than or equal to 10 should have decimal places based on the significant figures
                    dec_places = self.significant_digits - int(floor(log10(rounded_val))) - 1

                # For val >= 10, do not show decimal places if the rounded value is an integer
                if rounded_val >= 10 and float(rounded_val).is_integer():
                    return str(int(rounded_val))

                # Formatting the rounded value with the determined number of decimal places
                rounded_val_str = f"{rounded_val:.{dec_places}f}" if dec_places > 0 else f"{int(rounded_val)}"
                return rounded_val_str

        except ValueError:  # Catch non-numeric values
            if type_of_value == "EC50":
                return self.ec50_empty
            elif type_of_value == "pEC50":
                return self.pec50_empty
            else:  # type_of_value == "SEM":
                return ""

    @staticmethod
    def text_color_for_bg(bg_color):
        if bg_color == "#000000":  # Assuming 'none' is the color passed for NaN cells
            bg_color = "#D3D3D3"

        # Convert the hex color to RGB values
        r, g, b = tuple(int(bg_color[1:][i: i + 2], 16) for i in (0, 2, 4))

        # Calculate brightness
        brightness = (0.299 * r) + (0.587 * g) + (0.114 * b)

        # Return white for dark backgrounds, black otherwise
        return "white" if brightness < 128 else "black"

    def read_and_prepare_data(self, file_path: str):
        # Read the CSV
        data: DataFrame = pd.read_csv(file_path, sep="\t", header=None)

        rows, columns = data.shape

        self.rows = int(rows)
        self.columns = int((columns - 1) / 3)

        # Extract row names
        compounds: list[str] = data[0].tolist()

        # Extract data lists
        pec50_numeric: list[Union[int, float]] = data.iloc[:, 1::3].replace("NaN", np.nan).astype(float).values
        sem_list: list[Union[int, float]] = data.iloc[:, 2::3].replace("NaN", np.nan).astype(float).values
        n_array = data.iloc[:, 3::3].replace("NaN", np.nan).values
        n_list: ndarray[Any, dtype[Any]] = np.where(np.isnan(n_array), np.nan, n_array.astype(int))

        # Create array for heatmap annotations
        data_info = []
        for pEC50_row, SEM_row, n_row in zip(pec50_numeric, sem_list, n_list):
            row_info = []
            for pEC50, SEM, n in zip(pEC50_row, SEM_row, n_row):
                # Convert pEC50 to EC50
                ec50_value: Union[int, float] = pow(10, -1 * pEC50) * 1000000000
                ec50_str = self.round_it(ec50_value, "EC50")
                pec50_str = self.round_it(pEC50, "pEC50")
                sem_str = self.round_it(SEM, "SEM")
                try:
                    n = int(n)
                except ValueError:
                    pass

                if sem_str:  # If sem_str is not empty
                    row_info.append(f"{ec50_str}\n{pec50_str}±{sem_str} (n={n})")
                else:
                    row_info.append(f"{ec50_str}\n{pec50_str}")
            data_info.append(row_info)

        # Convert arrays to DataFrames
        df_color = pd.DataFrame(
            pec50_numeric, columns=self.column_names, index=compounds
        )
        df_info = pd.DataFrame(
            data_info, columns=self.column_names, index=compounds
        )

        return df_color, df_info, pec50_numeric

    def replace_values(self, df_color, df_info):
        # self.modifications

        for cell_coords, mod_details in self.modifications.items():
            df_color.iat[cell_coords] = mod_details['color']
            df_info.iat[cell_coords] = mod_details['text']
        return df_color, df_info

    @staticmethod
    def create_colormap():
        colors = ["#363636", "#4D91BE", "#F8CD3D"]
        return LinearSegmentedColormap.from_list("custom", colors, N=256)

    @staticmethod
    def calculate_figure_dimensions(n_rows):
        row_height = 16 / 24
        return n_rows * row_height

    @staticmethod
    def draw_heatmap(df_color, ax, cmap):
        sns.heatmap(
            df_color, cmap=cmap, cbar=False, ax=ax, linewidths=2, linecolor="black"
        )

    def add_colorbar(self, fig, ax, df_color):
        # Extract position details of the heatmap
        pos = ax.get_position()
        cbar_width = 0.02
        gap = 0.01

        # Define colorbar axes based on main heatmap position
        cbar_ax = fig.add_axes([pos.x1 + gap, pos.y0, cbar_width, pos.height])

        # Create a colorbar
        norm = plt.Normalize(df_color.min().min(), df_color.max().max())
        cbar = fig.colorbar(
            plt.cm.ScalarMappable(norm=norm, cmap=self.cmap),
            ax=ax,
            cax=cbar_ax,
            ticks=np.arange(
                np.floor(df_color.min().min()), np.ceil(df_color.max().max()) + 1
            ),
        )
        cbar.outline.set_edgecolor("black")
        cbar.outline.set_linewidth(2)
        return cbar

    @staticmethod
    def customize_spines(ax):
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(2)
            spine.set_edgecolor("black")

    @staticmethod
    def set_axis_labels(ax, df_color):
        ax.set_xticklabels(df_color.columns, rotation=0, horizontalalignment="center")
        ax.set_yticklabels(df_color.index, rotation=0, verticalalignment="center")

    def calculate_colors_for_annotations(self, df_color):
        # Normalize data for color mapping
        norm_data = (df_color - df_color.min().min()) / (df_color.max().max() - df_color.min().min())
        colors_rgb = self.cmap(norm_data.to_numpy().flatten())

        # Convert RGB colors to hex format
        colors_hex = [
            "#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255))
            for r, g, b, _ in colors_rgb
        ]

        colors_hex_array = np.array(colors_hex).reshape((self.rows, self.columns))

        return colors_hex_array

    def add_text_annotations(self, ax, df_info, mask, colors_hex, cells_to_skip):
        # Assume colors_hex is a list of lists with the same structure as df_info
        for i in range(df_info.shape[0]):
            for j in range(df_info.shape[1]):
                if mask[i, j] or (i, j) in cells_to_skip:
                    # If the cell is flagged in the mask or is in the list of cells to skip, skip over it
                    continue

                # Get the background color for the cell from colors_hex
                cell_bg_color = colors_hex[i][j]  # Note: this assumes colors_hex is a list of lists
                # Determine the appropriate text color based on the cell's background color
                text_color = self.text_color_for_bg(cell_bg_color)
                # Get the text to display for the cell
                text = df_info.iloc[i, j]
                # Add the text annotation to the cell with the determined text color
                ax.text(j + 0.5, i + 0.5, text, ha="center", va="center", color=text_color)

    def modify_cells(self, ax, modify_dict, colors_hex):
        for cell_coords, mod_details in modify_dict.items():
            row_idx, col_idx = cell_coords
            cell_text = mod_details.get('text', '')
            text_color = mod_details.get('text_color', 'black')
            fontweight = mod_details.get('fontweight', self.font_settings['weight'])
            fontsize = mod_details.get('fontsize', self.font_settings['size'])

            # Retrieve the cell's background color from the colors_hex array
            background_color = colors_hex[row_idx, col_idx]

            # Determine the appropriate text color if it wasn't explicitly provided
            if 'text_color' not in mod_details:
                text_color = self.text_color_for_bg(background_color)

            # Update the text annotation with custom properties
            ax.text(col_idx + 0.5, row_idx + 0.5, cell_text,
                    ha="center", va="center", color=text_color,
                    fontweight=fontweight, fontsize=fontsize)

    def get_masks_and_modifications(self, df_info, pec50_numeric):
        # Initialize the mask with False values, i.e., do not skip any cells by default
        mask = np.zeros(df_info.shape, dtype=bool)

        # Get the minimum and maximum pEC50 values used for the heatmap
        data_min = np.nanmin(pec50_numeric)
        data_max = np.nanmax(pec50_numeric)

        # Generate the modifications dictionary using a custom method
        # For the sake of this example, let's assume the custom_text method
        # generates a dictionary where the keys are (row_index, column_index) tuples
        # and the values are the modification details for each cell.
        modify_dict = self.custom_text(data_min, data_max)  # , pec50_numeric

        # Update the mask based on the keys in the modify_dict
        for cell_coords in modify_dict.keys():
            row_idx, col_idx = cell_coords
            mask[row_idx, col_idx] = True

        return mask, modify_dict

    def plot_heatmap(self, df_color, df_info, pec50_numeric):
        cmap = self.create_colormap()
        fig_height = self.calculate_figure_dimensions(df_color.shape[0])
        sns.set(font="Arial")
        fig, ax = plt.subplots(figsize=(10, fig_height))

        self.draw_heatmap(df_color, ax, cmap)
        self.add_colorbar(fig, ax, df_color)
        self.customize_spines(ax)
        self.set_axis_labels(ax, df_color)

        colors_hex = self.calculate_colors_for_annotations(df_color)
        mask, modify_dict = self.get_masks_and_modifications(df_info, pec50_numeric)
        cells_to_skip = self.cells_to_skip()
        self.add_text_annotations(ax, df_info, mask, colors_hex, cells_to_skip)

        plt.rcParams.update({'font.size': self.font_settings['size']})
        plt.rcParams.update({'font.family': self.font_settings['font']})

        self.modify_cells(ax, modify_dict, colors_hex)

        return fig, ax

    def cells_to_skip(self):
        return [key for key in self.modifications.keys()]

    def custom_text(self, data_min, data_max) -> dict[tuple, dict[str, str]]:
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

        modifications = {
            cell_coords: {
                "color": get_color_from_value(mod_details["color"], data_min, data_max, self.cmap)
                if isinstance(mod_details["color"], (int, float)) else mod_details["color"],
                "text": mod_details["text"]
            }
            for cell_coords, mod_details in self.modifications.items()
        }

        # Convert RGB colors to hex format if they are not already in hex
        for cell, mod_details in modifications.items():  # Changed 'modifications' to 'mod_details' here
            color = mod_details["color"]
            if isinstance(color, tuple):  # If the color is an RGB tuple
                mod_details["color"] = '#{:02x}{:02x}{:02x}'.format(
                    int(color[0] * 255),
                    int(color[1] * 255),
                    int(color[2] * 255)
                )

        return modifications

    def save_and_show(self, fig, output_path):
        saved_files = []
        if "png" in self.save_as_type:
            plt.savefig(output_path + ".png", transparent=True, bbox_inches="tight")
            saved_files.append(output_path + ".png")
        if "pdf" in self.save_as_type:
            fig.savefig(output_path + ".pdf", format="pdf", transparent=True, bbox_inches="tight")
            saved_files.append(output_path + ".pdf")
        if self.show:
            plt.show()
        saved_files_str = ""
        for i in saved_files:
            saved_files_str += "\n" + i
        end_message = f"The script has completed. The following files have been saved:\n" + saved_files_str + "."
        print(end_message)


heatmap = HeatmapGenerator(parameters)
heatmap.create_heatmaps()
