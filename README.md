# JAGWAS
Joint Analysis of multi-phenotype GWAS

<br />
Current version: 1.0.1

## Quick Installation 

Option 1: Use the binary executable file for Linux
* Download the binary executable JAGWAS from the main directory or run the following line of code:
```
git clone https://github.com/baihongguo/JAGWAS
cd JAGWAS
```
Then change the permission:
```
chmod a+x JAGWAS
```

Option 2: Build JAGWAS Library Dependencies  
   * C++11 compiler or later 
   * BLAS/LAPACK. For Intel processors, we recommend that JAGWAS be compiled with an optimized math routine library such as the Intel oneAPI Math Kernal Library to replace BLAS/LAPACK for optimal performance.  
   * Armadillo Library. 
<br />
<br />

## Dependencies
C/C++ Compiler
 * A compiler with C++11 (or later) support is required.
 
LAPACK and BLAS
 * The LAPACK (Linear Algebra PACKage) and BLAS (Basic Linear Algebra Subprograms) libraries are used for matrix operations in JAGWAS.

Armadillo Library
 * The Armadillo library is used for matrix operations, which is dependent on LAPACK and BLAS.


## Usage

### Running JAGWAS

1. [Command Line Options](#command-line-options)  
1. [Input File Format](#input-file-format)
1. [Output File Format](#output-file-format)
1. [Examples](#examples)

<br /> 
<br />

### Command Line Options

Once JAGWAS is installed, the executable can be used to run the program.  
```
Input Options:
--outputFilePath
  Path to the output file.

--cor_matrix
  A text file including the correlation matrix of the phenotypes, and the order of the phenotypes in the correlation matrix should match with the order of the summary statistics files in the input command. 

--nrow
  The number of rows of results each time JAGWAS processes and writes out, may affect memory usage, 1000-50000 is recommended. The maximum is the total number of SNPs.

--MAF
  Minimum threshold value [0, 0.5] to exclude variants based on the minor allele frequency.

--score_test
  If the summary statistics files are from a score test, 1 = True, 0 = False. 0 means the summary statistics files are from a Wald test.

--beta_se
  If the summary statistics files are in a beta/se format from the Wald test, 1 = True, 0 = False. When score_test=0, and beta_se=0, this means that the summary statistics files have the z-scores.

--logP
  Whether to output the p-values in log10 scale, thus can getting p-values smaller than the double precision limit, 1 = True, 0 = False.

--delim
  Delimiter separating values in the summary statistics file. Tab delimiter should be represented as \t, space delimiter as \0, comma delimiter as ,.

--fileNames
  Paths to the summary statistics files.
```

<br /> 
<br />

### Input File Format
The input summary files can be tab, single-space, or comma-separated.
* #### Score test File Format
The output from the GMMAT score test is set as the standard input format by JAGWAS:
```diff 
SNP	A1	A2	N	AF	SCORE	VAR	PVAL
1	AC	A	5127	0.389814452896431	34.4751108518404	1047.78879525852	0.286854635988476
2	G	A	5127	0.113190081919251	-18.3952465352991	446.694169234623	0.384102004874046
3	A	G	5127	0.0446306807099668	13.2671776127244	414.396331767185	0.514572577835737
4	C	T	5128	0.251488202028081	-32.4559032611082	1963.18654661393	0.463858002357953
```


* #### Wald test File Format
The output from the fastGWA Wald test is set as the standard input format by JAGWAS:
```diff 
CHR	SNP	POS	A1	A2	N	AF1	BETA	SE	P	INFO
1	rs1	665266	T	C	2288	0.981101	-1.29776	2.95353	0.660376	0.664293
1	rs2	714596	T	C	2288	0.96833	0.786756	2.06255	0.702871	0.84736
1	rs3	715265	C	T	2288	0.964588	0.66013	1.86351	0.723159	0.932847
1	rs4	715367	A	G	2288	0.964621	0.873248	1.86112	0.638923	0.93636
```

or

```diff 
CHR	SNP	POS	A1	A2	N	AF1	BETA	SE	Zscore
1	rs1	665266	T	C	2288	0.981101	-1.29776	2.95353	2.5	
1	rs2	714596	T	C	2288	0.96833	0.786756	2.06255	3.6	
1	rs3	715265	C	T	2288	0.964588	0.66013	1.86351	-2.5	
1	rs4	715367	A	G	2288	0.964621	0.873248	1.86112	3.6	
```

Only AF, BETA/SE, and SCORE/VAR or Zscore are a must-have. If any of the other columns are missing from your summary statistics files, simply fill them out using NA. 
Note that all the input summary statistics files must include the same sets of variants and also in the same order. 

<br /> 
<br />

### Output File Format  
```diff
CHR	SNP	POS	A1	A2	N	AF1	P
1	rs1	665266	T	C	2288	0.981101	0.875815
1	rs2	714596	T	C	2288	0.96833	0.724438
1	rs3	715265	C	T	2288	0.964588	0.855705
1	rs4	715367	A	G	2288	0.964621	0.863347
```


### Examples  
<br />

To run JAGWAS using the example data, execute JAGWAS with the following code.
```unix
./JAGWAS --outputFilePath outputFilePath/example_JAGWAS.txt --cor_matrix pathTo/correlation_matrix.txt --nrow 3 --MAF 0 --score_test 0 --beta_se 1 --logP 0 --delim \t --fileNames pathTo/discovery_QT0.txt pathTo/discovery_QT1.txt

```
The results should look like. 
```diff
CHR	SNP	POS	A1	A2	N	AF1	P
1	rs539032812	665266	T	C	22880	0.981101	0.875815
1	rs149887893	714596	T	C	22880	0.96833	0.724438
1	rs12184267	715265	C	T	22880	0.964588	0.855705
1	rs12184277	715367	A	G	22880	0.964621	0.863347
1	rs12184279	717485	C	A	22880	0.964754	0.818601
1	rs144155419	717587	G	A	22880	0.983974	0.569762
1	rs116801199	720381	G	T	22880	0.964323	0.811919
1	rs12565286	721290	G	C	22880	0.964186	0.793858
1	rs2977670	723891	G	C	22880	0.0389927	0.960416
```
<br />

## Recent Updates 
[Version 1.0.1] - May 26, 2024:
* Change the output variable log10(P) to log10_P to fix the issue that results$log10(P) not working in R since log10() is a built-in function.
* Add a delim argument to give more flexibility to the format of summary statistics files, user now can choose tab, single-space, or comma as the delimiter.  

[Version 1.0.0] - April 14, 2024:
* The original version of JAGWAS

## References
<p>If you use JAGWAS, please cite
<li>Guo, B. et al. Efficient multi-phenotype genome-wide analysis identifies genetic associations for unsupervised deep-learning-derived high-dimensional brain imaging phenotypes. <em>medRxiv</em> 2024.12.06.24318618 
[<b>Preprint</b>]. Dec 8, 2024. DOI: <a href="https://www.medrxiv.org/content/10.1101/2024.12.06.24318618v1">10.1101/2024.12.06.24318618</a>.</li></p>
