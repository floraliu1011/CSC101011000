## Intro  
This is Flora's Repo for ROP299 project. 
  
The Repo contains following scripts:  

 1. `setup.sh`   
 	 * Autometically splits up large file and generates cluster job submission file  
 	 * Orders file directories as it builds, leading to more organized workflow. It builds following subdirectories:
 	 	* `/Data` contains the original data files
 	 	* `/Intermediate` contains the split file
 	 	* `/Output` contains the split result and logging files
 2. `/Parser`  
	* Contains scripts for parsing sequences implemented in Python3  
	* Separates IO from main functionalities, alow easier extension and modification
 3. `/Analyser`
	* Contains scripts for enumerating parsed elements implemented in R
	* `/Analyser/analyse.R` takes in 

 ##Utils  
 * For downloading the Repo, copy the following command into the command line  
```git clone https://github.com/floraliu1011/ROP299.git```
 **[git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) needs to be installed in order to run this command**  
 
* Move the Repo, as well as the sequencing result, to the scratch folder. Make sure that the forward read file and and reverse read file is located in the same directory as `setup.sh` and `/Script`. Run the following command to invoke the building process:
```sh setup.sh <path_to_forward_read> <path_to_reverse_read> [num_split]```

* Use the following command to submit job  
	```qsub submission```  
	
* Inspect the result under `/Output`  

* Open `Analyser/analyse.R` in [R Studio](https://www.rstudio.com) and fill in the path inside the comment block to obtain output data

##Extension and Modification

To change the behavior of the parser, change the funcitons in `Script/Parser/utils.py`. As long as the functionality of `parse` is preserved  
To modify the submission file, modify the `setup.sh`

 ##Potential Bugs & Attention

1. The submission files (`submission`) assembled by `setup.sh` might contain line break inside job commands.  
For example:  
```
(sh script.sh file_A     
file_B file_C)&
```  
might leads to early abortion.  
Please ensure the line breaks are removed inside each job command.
2. The file is preferably been split into a multiple of 16 to obtain better use of the computation power

 ##Further reading

* [Scinet User tutorial](https://wiki.scinet.utoronto.ca/wiki/images/5/54/SciNet_Tutorial.pdf)
* [Scinet Wiki](https://wiki.scinet.utoronto.ca/wiki/index.php/SciNet_User_Support_Library)
* [R package data.table](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html)
* 