# SeqAlign2
Shiny app for multiple sequence alignment.

You can just run this R code from R console or R studio.

You can use this tool to find if our CRISPR/Cas9 transgenic lines contain targeted mutations, or use it to check if your plasmids contain the sequence you want.

For CRISPR/Cas9 transgenic lines, we suggest to design your sequence primers to have the same directionality as the target sequence. This would take sometime to explain. Please just do so. 

The user interface is self-explanatory. You would need to input the target sequence. Just copy and paste. You can select if the target sequence needs to be reversed and complemented. Then you can use "Choose files" to select files from proper directory. For this step, you can also drag multiple sequence files that ends with the .seq directly onto the "Choose files" botton. Next you can select which range of sequences you would like to use. The reason to do this is to remove low quality nucleotides that usually are located at the two ends. You can see if an output that whether the target sequence exists in the sequence.

You can download the pdf and csv files.  
  * The csv file contains: 
    1. The name of the input file 
    2. Length of DNA sequences used in this analysis
    3. Whether the target sequence is in the DNA sequences
    4. DNA sequences used in the analysis.
  * The pdf file contains the msa plot.
