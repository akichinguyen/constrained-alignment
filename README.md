1. Download ksw2-master and mummer-3.9.4alpha. 
2. Put all the folder and files in folder ksw2-master (exp: ksw2-master/test_ref/test)<br />
2.1. To reproduce the results reported in the paper, run the command <br />
- ./test_hash 100
- python time_graph.py<br />
2.2. To create a new test:<br />
- delete all files in /test and /test/result 
- change directory in file run_mummer.sh to your mummer-3.9.4alpha directory 
- run the bash script run_mummer.sh
