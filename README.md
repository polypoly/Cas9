## SeqAlign2
### Shiny app for multiple sequence alignment using the MSA package from bioconductor.

You could run the SeqAlign2.R from R studio and view it from your default web browser.

You could use this tool to find if your CRISPR/Cas9 transgenic lines contain targeted mutations, or use it to check if your plasmids contain the sequence you want.

For CRISPR/Cas9 transgenic lines, we suggest designing your sequencing primers to have the same directionality as the target sequence. The reason is due to the biallelic nature of Cas9 mutants. CRISPR/Cas9-generated mutants are often small deletions or insertions (INDELs). Two different INDELs (or one wild-type and one INDEL for the two alleles) for the two alleles would lead to confusion for Sanger sequencing. By using sequencing primers with the same directionality as the target sequence, you can see the INDEL events happened within the target sequences. By using sequencing primers with the opposite directionality as the target sequence, you would not be sure where the INDEL events occur.

The biallelic sequences could also be observed in sequencing chromatograms (usually saved in .ab1 format).  

The user interface is self-explanatory. You would need to input the target sequence. Just copy and paste from MS word or excel sheet. You can select if the target sequence needs to be reversed and complemented. Then you can use "Choose files" to select files from a proper directory. For this step, you can also drag multiple sequence files that end with the .seq directly onto the "Choose files" button. Next, you can select which range of DNA sequences you would like to use. The reason to do this is to remove low-quality nucleotides that usually are located at the two ends. The next step is to give the output file a name. The default name is SeqInfo. Please give a meanful name for your output file. In the Shiny main panel, you can see whether the target sequence exists in the DNA sequences. In fact, the main panel info is exactly the same as the CSV file that you can download.

## Output files:
You can download the pdf and CSV files.  
  * The CSV file contains: 
    1. The name of the input file
    2. Length of DNA sequences used in this analysis
    3. Whether the target sequence is in the DNA sequences
    4. DNA sequences employed in the study
  * The pdf file contains the MSA plot.

## Example:
For target sequence, use CGCGTCTTGTCGAACGAAGCCT.  
For example DNA sequences, find them in the test.zip folder.  
The output should be identical to the SeqInfo.pdf and SeqInfo.R files in this repository. 

## References:
Bodenhofer U, Bonatesta E, Horejs-Kainrath C and Hochreiter S (2015). “msa: an R package for multiple sequence alignment.” Bioinformatics, 31(24), pp. 3997–3999.

## Version info:
RStudio Version 0.99.903
MSA Version 1.4.5


Have fun!
