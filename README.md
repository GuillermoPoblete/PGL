# PGL

How to use PGL

1- Install R Language https://cloud.r-project.org/
1- Download the files from this repository, and place it in a drive with at least 3.1 GB of free disk space.
2- Copy the genes and/or ROIs lists to be processed in the same folder as the PGL script.
2- Run the Source command so you can access to the PGL's functions from the comman line:
  source('PGLAPI.GitHub.R')
3- Run the command GetROIs('list.csv') or GetGenes('list.csv')
The first time you run one of these commands, a customized version of the Allen Human Brain Atlas will be download. This will require an Internet conection and 3.1 GB of free disk space.
4- When the process is completed, you will find a CSV file containing the results in the same directory where you placed the PGL script and the lists to be processed. This CSV file is named as follows: The list name + 'PGL Result' + the date and time when the function ended + '.csv'




