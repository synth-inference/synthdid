### Synth MATLAB Code (11/07/2006) written for MATLAB 7.0
### by Alberto Abadie, Alexis Diamond, and Jens Hainmueller (all Harvard University)
### Contact: jhainm@harvard.edu

The Synth MATLAB implements synthetic control methods for causal inference
in comparative case studies with aggregate data as developed in Abadie and
Gardeazabal (2003) and Abadie, Diamond, and Hainmueller (2006).

Files:
synth_code.m --> Main Script that runs the example and reproduces the main
                 results of the paper (the results may differ
                 slighly due to tolerance settings of optimization). 
                 Researchers are encouraged to adjust this code for their
                 needs; please give due credit.


loss_function.m --> loss function called by synth_code.m

MLAB_data.txt --> 39 by 39 data matrix to run the example:
                   First row contains state numbers:
                     Alabama 1; Arkansas 2; Colorado 4; Connecticut 5; Delaware 6;
                     Georgia 7; Idaho 8; Illinois 9; Indiana 10; Iowa 11; Kansas 12;
                     Kentucky 13; Louisiana 14; Maine 15; Minnesota 16; Mississippi 17;
                     Missouri 18; Montana 19; Nebraska 20; Nevada 21; New Hampshire 22;
                     New Mexico 23; North Carolina 24; North Dakota 25; Ohio 26; Oklahoma 27;
                     Pennsylvania 28; Rhode Island 29; South Carolina 30; South Dakota 31;
                     Tennessee 32; Texas 33; Utah 34; Vermont 35; Virginia 36; West Virginia 37;
                     Wisconsin 38; Wyoming 39; California 3.
                  Predictors are stored in rows 2 to 8:
                     row 2: income, row 3: retail price, row 4: percent_15-19; row 5: beer
                     consumption (all averaged 1980 to 1988);
                     row 6: smoking 1988, row 7: smoking 1980; row 8: smoking 1975;
                  Outcome Data (smoking consumption in packs per capita) is stored in rows 9 to 39
                     for the years 1970, 1971,...,2000


