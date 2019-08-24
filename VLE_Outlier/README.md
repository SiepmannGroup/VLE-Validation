usage: python outlier.py input weight t_c

To run the sample input file:
python outlier.py sample_input.txt 74.12 523.4

input: the path to the input file containing VLE data
weight: the molecular weight of the molecule in g/mol
t_c: the critial temperature of the molecule in K

The data in the input file is parsed using the column names in the first line:

temperature: <code>T / temperature</code>

vapor density: <code>rhog / density</code>

pressure: <code>p / pressure</code>

compressibility factor: <code>z / compressibility</code>

