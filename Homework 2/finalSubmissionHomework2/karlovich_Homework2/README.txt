README.txt

The command to run the file is
python3 hw2_karlovich.py -q {Query File} -r {Reference File} -m [Optional Matches .txt file]

The query and reference files must put the query on the second line of the file as the
implementation here skips the first line of the file (which is often header information)

So an example query file would look like

1. Query ABCD, random other garbage
2. FSGHSGBUEFHAFJGSHGEUSGSGSEGUSGHSEGUSFHSAGUAFUAF
3. {blank}
...
where the first line is whatever information the second line is the query and after
the query there is nothing else.

The Python file hw2_karlovich.py should be placed in the folder where it can access
all of the query, reference and matches files.

The query and reference sequences are required but the matches parameter is optional
If you include a match file (with the signal -m) then the program will run
and anchored version of needleman_wunsch.  If no values is provided, the script
will run a non-anchored version of needleman_wunsch

At the top of the python script there is a global variable OUTPUT_FILE which
is by default set to "output.txt" but can be changed by the user.

Multi-processing mode: Lines that need to be edited/looked over: 35-55, 125, 279-285
In the __main__ function there is a line of commented code, uncomment that code and
comment the three lines of code below it to switch into bulk multi-processing mode
When you do this you must edit the header of needleman_wunsch so that it's second
variable "seq2" has a default value equal to the reference you'd like to run another
sequence against.  Currently the multi_processing function will take in seq1, permute
the sequence NUM_OF_ITERATIONS times and then return a list of the score of the
10,000 runs.

In the code comments there is an example of the function header that
needleman_wunsch() needs to look like in order to run Multi-processing.  When you are done
it's not a bad idea to change the needleman_wunsch header back into requiring two variables.

The [default] output of the program is formatted as such:

Score: {Score}
Path Back: {path back consisting of "d", "v", and "h", for diagonal, vertical and horizontal respectively}

-------------------------
Sequence 1: {the value on line 2 of the query file}
Sequence 2: {the value on line 2 of the reference file}
--- Aligned Sequences ---
Seq1: {sequence 1 but with gaps that are aligned to seq2}
Seq2: {sequence 2 but with gaps that are aligned to seq1}
