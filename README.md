# HeatmapGenerator

This is a script to create heatmaps for publications.

This is an example for how it can look generated from the date in "data/example.csv".
![Heatmap Example](https://github.com/Aephir/create_heatmap/blob/main/data/example_heatmap.png?raw=true)

In the current version (v2.0.0), the following has been implemented:

Generating heatmap images from a csv file with data divided in sets of 3 sub-columns for each (plus one column for row names).

The 3 sub-columns should have:

1) pEC50 value as a float or integer.
2) Error range (± CI95, SEM, etc.) as a float or integer.
3) Number of independent repeats (n) as an integer.

There can be as many rows as needed (e.g., different compounds or conditions) with row names in the first column, and there should be no column headers. 

## Adjustable parameters 
Parameters that can be adjusted are in the "adjustable_parameters" dictionary. Current options are:

### files
A list of the files cvs you want to create heatmaps from.

### significant_digits
The number of significant digits to use (NOT significant decimals)

__Caveat:__ For pEC50 values, "significant_digits" is number of decimals.

### column_names
Names of column headers

### pec50_empty
The text to use for pEC50 if no value in cell

### ec50_empty
The text to use for EC50 if no value in cell

### modifications
Allow custom color and text for specified cells. Use case could be:
- If you cannot dissolve a compound enough to use the standard "pec50_empty" and "ec50_empty".
- If you can determine approximate values, and want to indicate this.

The tuple (e.g., (0, 1)) determines the location in the heatmap to change (index, column).

_Remember to start counting at 0!_.

For example:
- (0, 1) means index=0 (or row=0), column=1, so first rw, second column.
- (23, 4) means index=23, column=4, so 24th row, fifth column.

You can give the "color" parameter either as a HEX value color code (e.g., "#000000" for black) or as a numeric value where the color will become the closest available in the current heatmap.

E.g., giving a value of 6 where the heatmap values/colors range is 5–9 will give the color corresponding to 6. A value of 4 in the same heatmap will ive the color corresponding to 5 (closest to 6 available)

### font_settings
Settings for the font type, size, and weight (e.g., "normal", "bold") to use

### show_plot
Show the plot after the script is done (as opposed to just saving the file)

### save_as_type
Which types o files should be saved. Current options are "png" and "pdf".


## Known limitations
- Background "color" for empty cells is "clear", and there is a white "plus" through the cells. While this does not matter i exporting a png and displaying on white background, it limits other uses.
