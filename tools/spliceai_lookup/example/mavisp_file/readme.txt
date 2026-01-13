## Usage 

### API activation

Activate the local API.

In order to open the local API open a new terminal and run the following commands:
```
source /data/user/marnaudi/spliceai_lookup/spliceai-env/bin/activate
conda activate /home/marnaudi/.conda/envs/bioenv/
```
run the following command to check if there are some running process:

```
lsof -i :8080 
``` 
you need to kill all the processes with your name in the port 8080

To kill them, check the PID code and run the following command for each process to kill
```
kill -9 PID
```
Check that no process is active with
```
lsof -i :8080
```
Activate the local API server with
```
cd /data/user/marnaudi/spliceai_lookup/SpliceAI-lookup/
bash start_local_server.sh

```

### splice_lookup.py 
```
module load python/3.10/modulefile
python splice_lookup.py -i spliceai_input_file.csv  -d 500 -t 9
```


You can run an example within the example folder
```
bash run.sh
```
