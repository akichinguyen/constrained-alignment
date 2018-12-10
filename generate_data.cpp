#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
// #include "final.h"
using namespace std;

string generate_ref (string query, int max_len_insert_string, int max_insert_position){
	string nu[4] = {"A","C","G","T"};
	int number_insert_position = rand()%max_insert_position + 1;
	for(int i = 1; i<=number_insert_position; i++){
		int position = rand()%(query.length());
		int len_insert_string = rand()%max_len_insert_string +1 ;
		string s;
		for (int j = 1; j<=len_insert_string; j++){
			s += nu[rand()%4];
		}
		query = query.substr(0,position) + s + query.substr(position, query.length()-position);
	}
	return query;
}
void write_ref(int number_ref, string query, int max_len_insert_string, int max_insert_position){
	for(int i = 0; i<number_ref; i++){
		string ref = generate_ref(query,max_len_insert_string, max_insert_position);
		// cout <<ref<<endl;
		string file_dir = "./test_ref/test/ref_"+to_string(i)+".fa";
		ofstream file(file_dir);
		file << ">ref_"+to_string(i)<<endl;
		file << ref;
		file.close();
	}
}
void write_query(string query){
	string file_dir = "./test_ref/test/query.fa";
	ofstream file(file_dir);
	file << ">query"<<endl;
	file << query;
	file.close();
}
string read_query(string filename){
	ifstream file(filename);
	string line, id, genome;
	while(getline(file, line)){
		// cout << line << endl;
		if(line.empty()){
			continue;
		}
		if(line[0] == '>'){
			id = line.substr(1);
		}
		else{
			genome += line;
		}
	}
	file.close();
	return genome;
}
int main(int argc, char*argv[]){
	int number_ref=atoi(argv[1]);
	string query = read_query("./test/MT-human.fa");
	query = query.substr(0,300);
	write_query(query);
	write_ref(number_ref,query, 2, number_ref);
}

