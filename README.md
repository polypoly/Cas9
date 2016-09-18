## SeqAlign2
### Shiny app for multiple sequence alignment.

You can just run the SeqAlign2.R from R console or R studio.

You can use this tool to find if our CRISPR/Cas9 transgenic lines contain targeted mutations, or use it to check if your plasmids contain the sequence you want.

For CRISPR/Cas9 transgenic lines, we suggest designing your sequence primers to have the same directionality as the target sequence. The reason would take some time to explain. Please just do so. 

The user interface is self-explanatory. You would need to input the target sequence. Just copy and paste. You can select if the target sequence needs to be reversed and complemented. Then you can use "Choose files" to select files from a proper directory. For this step, you can also drag multiple sequence files that end with the .seq directly onto the "Choose files" button. Next, you can select which range of DNA sequences you would like to use. The reason to do this is to remove low-quality nucleotides that usually are located at the two ends. The next step is to give the output file a name. The default name is SeqInfo. Please give a meanful name for your output file. You can see if an output that whether the target sequence exists in the DNA sequences.

You can download the pdf and CSV files.  
  * The CSV file contains: 
    1. The name of the input file 
    2. Length of DNA sequences used in this analysis
    3. Whether the target sequence is in the DNA sequences
    4. DNA sequences employed in the study.
  * The pdf file contains the MSA plot.

## Example:
For target sequence, use CGCGTCTTGTCGAACGAAGCCT.
For example DNA sequences, find them in the test.zip folder
The output should be identical to the SeqInfo.pdf and SeqInfo.R files in this repository. 

Have fun!
