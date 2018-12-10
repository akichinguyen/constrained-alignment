#! /bin/bash
echo Please enter number of references you want
read number_ref
g++ generate_data.cpp -o generate_data
./generate_data $number_ref
i=0
while [ $i -lt $number_ref ]
do
	/home/thunguyen/Downloads/mummer-3.9.4alpha/./mummer /home/thunguyen/Downloads/ksw2-master/test_ref/test/ref_$i.fa /home/thunguyen/Downloads/ksw2-master/test_ref/test/query.fa >> /home/thunguyen/Downloads/ksw2-master/test_ref/test/result/result_$i.txt
	((i++))
done
g++ test_hash.cpp ksw2_extz.c ksw2_gg.c -o test_hash
./test_hash $number_ref
python time_graph.py