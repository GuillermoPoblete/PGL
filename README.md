## How to use PGL

1. Install R Language https://cloud.r-project.org/
2. Download the files from this repository, and place them in a drive with at least 3.1 GB of free disk space.
3. Copy your genes and/or ROIs lists to be processed in the same folder as the PGL script.
4. Run the Source() command so you can access the PGL's functions from the command line:
```
source('PGL.R')
```
5. Run the command 
```
Result<-GetROIs('ROIs.csv'), replacing the filename for your own ROIs list file name.
```
or 
```
Result<-GetGenes('Genes.csv'), replacing the filename for your your genes list file name.
```
You can also use the ROIs and genes lists uploaded in this repository to try PGL for the first time.

The first time you run one of these commands, a customized version of the Allen Human Brain Atlas will be download, so you will need an Internet connection and 3.1 GB of free disk space.

When the process is completed, you will find a CSV file containing the results in the same directory where you placed the PGL script and the lists to be processed. This new CSV file is named as follows: The original list name + 'PGL Result' + the date and time when the function ended + '.csv'
