# HeatmapGenerator

This is a script to create heatmaps for publications.

This is an example for how it can look generated from the date in "data/example.csv".
![Heatmap Example](https://github.com/Aephir/create_heatmap/blob/main/data/example_heatmap.png?raw=true)

In the current version (v2.0.6), the following has been implemented:

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
Hint: Currently it prompts you for which experiment to calculate from. This will select the correct file (if keeping Excel file names as indicated in the script)

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

__*Remember to start counting at 0!*__.

For example:
- (0, 1) means index=0 (or row=0), column=1, so first rw, second column.
- (23, 4) means index=23, column=4, so 24th row, fifth column.

You supply a numeric value as the 'color' variable, that will act as if that value was in your dataset, 

You supply any string as the 'text' variable, that will be printed in the corresponding cell. Remember, 'backlash + n' `\n` is a "new line" character to have your text on two lines.  

### font_settings
Settings for the font type, size, and weight (e.g., "normal", "bold") to use

### show_plot
Show the plot after the script is done (as opposed to just saving the file)

### save_as_type
Which types o files should be saved. Current options are "png" and "pdf".

### Custom text and color for specific cells
This is set in the `const.py` file. Usefull for when you have incomplete curves and want to show that EC50 is above or approximate.
Hint: Currently it prompts you for which experiment to calculate from. This will select the correct variable to import from `const.py`, as long as you edit the variable to fit your data.

## Known limitations
- Background "color" for empty cells is "clear", and there is a white "plus" through the cells. While this does not matter if exporting a png and displaying on white background, it limits other uses.
