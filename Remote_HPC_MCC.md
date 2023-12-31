# Remote cluster (MCC)

## Connect
`ssh user_name@mcc-dtn.ccs.uky.edu` 

Once you connect, you will be working with two computers: your local computer and the remote server. When running commands, pay attention to which computer you are working with.

To disconnect from MCC, run `exit`

___

## Transferring files 

### Copy files from your computer to the cluster
This requires you to be in the directory on your local machine where the `local_file` is (i.e. you should be in a terminal shell on your local machine, not ssh'ed into the cluster).
`scp local_file user_name@server_address:/path/to/cluster/directory`

### Copy files from the cluster to your computer
Again, this is from a location in your local machine (not ssh'ed into the cluster).
`scp user_name@server_address:/path/to/cluster/remote_file /path/to/local/directory`<br>

If the destination is the current working directory, `scp user_name@server_address:/path/to/cluster/file.txt .` can be used instead of the full path.

**_Question:_** How do you modify the `scp` command to copy a directory?

### Another way to transfer files
The command `sftp` provides another way to transfer files, but this time in a more interactive way. It's basically just like ssh'ing, but you enter a special shell that allows you to `get` and `put` files. So if you're on your local machine in a directory with a file you want to transfer, you can go `sftp [username]@mcc.uky.edu` and will be prompted to enter your PW. Once you hit enter, you will be in your home directory but with the prompt `sftp> `. Now you can `cd` around your directory structure in MCC, `ls` and `pwd` to your heart's content, and find the location you want to put a file. Then you can simply type `put [filename]` (note, this is the file name on your local machine that you want to transfer to the cluster) and it will transfer. Likewise, to get a file off the cluster you can `get [filename]` and it will transfer to the local machine. 

___

## Modules
Some programs are already installed on LCC, however, they will not work until they are activated.

`module list` lists the activated programs onto your account.

`module avail` lists the programs installed on LCC. Use `module avail | grep "program_name"` to search for a specific program (note that most module names on are in lowercase). Alternatively, you can use `module spider "program_name"` to search for a specific program.

`module load module_name` activates the program

`module unload module name` deactivates/removes the program from your account

___

## Singularities
These are like modules. The full list of available programs is [here](https://ukyrcd.atlassian.net/wiki/spaces/UKYHPCDocs/pages/2920537/Software+list+for+singularity+containers+for+conda+packages+in+the+MCC+cluster)

Instead of loading them, however, add the singularity information (in the "Notes" area) to your batch script. 

___

## Slurm job submission & batch scripts
LCC uses Slurm to manage user demand and system resources. 

Commands, or "jobs", are submitted with batch scripts. Batch script files must end in `.sh` and must start with the following lines:

```
#!/bin/bash
#SBATCH --time=             # How long you want to use the resources
#SBATCH --nodes=            # How many processors; typically 1
#SBATCH --ntasks=           # How many cores on the processor
#SBATCH --partition=        # Processor type
#SBATCH --account=          # CPU hours are monitored and charged

Optional
#SBATCH --mail-type=END     # Notification email when job is done
#SBATCH --mail-user=        # Your email address
```
The format for --time= is day-hour:minute:second, e.g., 00-01:00:00 means 1 hour.

If you choose to have notification emails sent to you, note that they often end up in spam/junk folders.

**_Task:_** Make a new file named "batch_header.sh" and add the following lines. We will use this file as a template for making batch files.
```
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --time 01-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --account=coa_jdu282_brazil_bootcamp2023
#SBATCH --mail-type END
#SBATCH --mail-user <your email address>
```
___

Once you submit a job, you can check it's status with `squeue | grep lcc_user_name`. If all processors are in use, your job will wait until resources are available. 

Each job gets a job number and a corresponding `slurm-job_number.out` file. This file has the information that is normally printed to the screen as the program is running. 

To cancel a submitted job use `scancel job_number` (get the job number from `squeue | grep mcc_user_name`)

NB: Never run a job without submitting a batch script.

### run a job
Now let's edit that `batch_header.sh` content to contain an actual command that we can watch. Create a new job file named `counting.sh`, and paste the following into it:
```
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --time 20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=coa_jdu282_brazil_bootcamp2023
#SBATCH --mail-type ALL
#SBATCH --mail-user [enter your email here]

echo starting
sleep 10
echo "it's been 10 seconds"
sleep 20
echo "it's been 30 seconds"
sleep 30 
echo "it's been 60 seconds. Exciting!"
sleep 5
echo "it's been 65 seconds. Amazing."
sleep 6
echo "it's been 71 seconds. What are we deviating from 5 second intervals???"
```
Save the file, and then check that the job submission file looks with cat. <br>
Let's submit the job using the following command: `sbatch counting.sh`

Once you have the job submitted, check that it's running by using `squeue`. Try that command out by itself. Lots of stuff, right? That's everyone's jobs that are running on the cluster right now. We can subset that in two ways. First, we could just grep out our username:
```
squeue | grep "username"
```
Or you can use an option in squeue:
```
squeue -u "username"
```

Let's check on the status of the job. Standard out (stdout) for computers is the normal output of a command/execution. Standard error (stderr) is any error messages that arise from a command/execution. The default location for stdout when MCC is running a job is in a file called `slurm-[jobID].out`. So let's see what's in that file. From when you submitted your job, or from squeue, you can figure out what the job ID is for that job. Try using `cat` to see what's in that file. How long ago did you submit that file? Can you piece it together based on what's in the output file?

Another way to interact with the stdout from a job is to write that output to another file. Take the `counting.sh` file and at the end of each "echo" line, add the following `  >> counting_output`. Resubmit the job, and now see if you can follow the status of the job in real time using `cat` and the `counting_output` file.

Also, note that I modified this line in the job script `#SBATCH --mail-type` to be `ALL`. You should have emails in your inbox that document when these jobs started and ended. The exit code is important in this email, as it lets you know if the job finished with no issues (exit code 0) or with an error or other issue. See [here](https://hpc-discourse.usc.edu/t/exit-codes-and-their-meanings/414) for more info on exit codes. They can be very useful when troubleshooting errors!
