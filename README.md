# Summary:
    Compare primary alignment scores from several SAM/BAM files. For each read (read pair if paired-end [PE]), output score from the respective input (blank if absent) as well as which one has the best score (multiple if ties). Three ways of calculating the score: 
    *) AS: original score reported in the tag 'AS';
    *) heuristic_raw: raw mapped length - penalty * #. of mismatches;
    *) heuristic_nonoverlap: nonoverlapping mapped length - penalty * #. of mismatches.

# Usage:
    perl $0 --inFiles <input1.bam,input2.bam,...> --labels <species1,species2,...> [--outfile output.txt]  [--endType PE] [--nameSorted] [--ncores 1] [--scheme AS] [--AS-per-mate] [--NM-per-mate] [--penalty 2] [--debug info] [--version]

# Options:
    --inFiles, -i   input files, comma separated, can be SAM or BAM format;
    --labels        labels of input files, comma separated;
    --outFile, -o   output file; if omitted, write to STDOUT; otherwise, if ending with '.gz', will be GZ compressed;
    --endType       choose from SE or PE (default: PE);
    --nameSorted    whether the input files are ALL name sorted. If not, they will be name sorted and saved to
                    temporary files with '.nameSorted.bam' suffix;
    --ncores        number of cores to use for name sorting (default: 1);
    --scheme        score scheme, choose from AS, heuristic_raw and heuristic_nonoverlap (default: AS);
    --AS-per-mate   boolean switch, if on, a read pair's AS equals to the sum of two mates' AS. This should be 
                    turned off if STAR, otherwise if Bowtie2 or BWA inputs (default: off); 
    --NM-per-mate   boolean switch, if on, instead of 'nM', a read pair's #. of mismatches equals to the sum of 
                    two mates' 'NM' (default: off)--this requires presence of tag 'NM'; note that STAR by default
                    does NOT output 'NM' (default: off); 
    --penalty       linear panelty of # of mismatches in heuristic scores (default: 2);
    --debug         level of logging, choose from fatal, error, warn, info, debug, trace (default: info);
    --version, -v   output version number.

# Output:
    The output has n+3 columns where n is the number of input files/labels:
        *) read ID;
        *) score from 1st input;
        ...
        *) score from n-th input;
        *) best score;
        *) source label(s) of the best hit(s);

# Example:
    perl sam_best_hits.pl --inFiles human.sam,mouse.sam,bacteria.sam --labels human,mouse,bacteria 2>err
    ReadID  human   mouse   bacteria        BestScore       BestLabel
    NB501328:197:HMK3KBGX7:1:11101:1498:13562       49      41              49      human
    NB501328:197:HMK3KBGX7:1:11101:1744:6598        39      47              47      mouse
    NB501328:197:HMK3KBGX7:1:11101:2468:19747       105     44              105     human
    NB501328:197:HMK3KBGX7:1:11101:2835:13039       35      35      31      35      human,mouse
    NB501328:197:HMK3KBGX7:1:11101:2927:13134               78      44      78      mouse
    NB501328:197:HMK3KBGX7:1:11101:3669:18422               105     107     107     bacteria
