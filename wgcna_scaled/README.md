# Mass WGCNA

Run WGCNA on 1000s of GSE id's and forget about it.

## What does this do?

Takes in a file containing GSE accessions from Gene Expression Omnibus(GEO), downloads it and runs WGCNA on them giving diagnostics on all the runs. Saves the modules for each GSE id and gives an option to add the modules to a GeneQuery module file.

## Usage

1. Clone the repo
2. Change directory 

    ```  
    $ cd masswgcna
    ```

3. Create a file containing the list of GSE Ids on which you want to perform WGCNA

    ```
    $ cat gse_list.txt
    ```
    
    GSE42261
    
    GSE425
    
    GSE43021
    
    GSE45430
    
    GSE5122
    
    GSE5654
    
    GSE57194
    
    GSE59808

4. Run the following command with the gse_list file path. Set the flag "-t" fas "ma" for Microarray datasets and "rna" for RNASeq datasets. 

    ```  
    $ ./mass_wgcna.sh -t ma -g gse_list.txt
    ```

5. After the run has completed, check the results of all GSE runs in the mass_run log file. This will tell you if a GSE id was processed successfully or it failed due to some reason.

    ```
    $ cat mass_run_2019.03.29-08.12.06.log
    ```
    
    2019-03-29 10:39:52 INFO:GSE104099:WGCNA successfully completed!!!
    
    2019-03-29 10:40:44 INFO:GSE106096:WGCNA successfully completed!!!
    
    2019-03-29 10:41:55 INFO:GSE12417:WGCNA successfully completed!!!
    
    2019-03-29 10:42:11 INFO:GSE12417:WGCNA successfully completed!!!
    
    2019-03-29 10:42:26 INFO:GSE12417:WGCNA successfully completed
    
    2019-03-29 10:42:46 ERROR:GSE13710:less than 20 samples present
    
    2019-03-29 10:42:47 ERROR:GSE13710:less than 20 samples present
    
    2019-03-29 10:44:25 INFO:GSE14468:WGCNA successfully completed!!!
    
    2019-03-29 10:46:45 INFO:GSE1729:WGCNA successfully completed!!!
  
6. Check why a particular dataset failed to run in the log directory.

    ```
    $ cd dataset_logs/
    ```
    
    ```
    $ cat GSE13710.log
    ```

7. All modules will be saved in a new directory called 'modules/' with a seperate module file for each GEO dataset.

8. All title info of GSE ids is saved in a "module-info.txt" file on the same path.

    ```
    $ cat module-info.txt
    ```


### Adding to a genequery module file
If you already have a genequery modules file and you want to append your modules to that file run the following 

```
$ ./add_modules.sh ~/genequery_rnaseq_integration/gqcmd-internal_new/hs.modules.gmt 
```

### Adding to a genequery GSE info file
If you already have a genequery GSE info file and you want to append your modules to that file run the following 

```
$ ./add_gse_info.sh ~/genequery_rnaseq_integration/gqcmd-internal_new/gse-info.txt 
```

### Configuring this on an instance
Run the script configure.sh with aws access key as the first argument and aws secret ket as the second argument to connect to the s3 bucket for this. This script will install R, dependencies and WGCNA. It will also set up aws cli with the credentials.

```
$ ./configure.sh <aws_access_key> <aws_secret_key>
```

### Caveats

A run may fail on a GSE id when processing of a platform annotation fails. This case is not that common(3 in 10).


