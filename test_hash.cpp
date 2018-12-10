#include <iostream> 
#include <fstream>
#include <sstream>
#include <unordered_map> 
#include <cstring>
#include <iterator>
#include <vector>
#include <algorithm>
#include <list>
#include <regex>
#include <chrono>
#include <sstream>
#include "ksw2.h"
// #include "final.h"

using namespace std; 
using namespace std::chrono;

// cache to store alignment results
unordered_map<int, unordered_map<int, unordered_map<string, string>>> genome_cache;

// finding alignment results from the cache. Return the alignment result if it exits, otherwise, return empty string
string cache_finding (int off_q, int len_q, string sub_ref){
	if(genome_cache.find(off_q)!= genome_cache.end()){
		unordered_map<int, unordered_map<string, string>> sub_q = genome_cache[off_q];
		if(sub_q.find(len_q) != sub_q.end()){
			unordered_map<string, string> sub_r = sub_q[len_q];
			if(sub_r.find(sub_ref) != sub_r.end()){
				return sub_r[sub_ref];
			}
		}
	}
	return "";
}

// if the alignment does not exist in the cache, store it for later queries
void add_to_cache (int off_q, int len_q, string sub_ref, string alignment){
	if(genome_cache.find(off_q) == genome_cache.end()){
		unordered_map<string, string> res;
		res.insert(make_pair(sub_ref, alignment));
		unordered_map<int, unordered_map<string, string>> sub_r;
		sub_r.insert(make_pair(len_q, res));
		genome_cache.insert(make_pair(off_q,sub_r));
	}
	else{
		unordered_map<int, unordered_map<string, string>> sub_r = genome_cache[off_q];
		if(sub_r.find(len_q) == sub_r.end()){
			unordered_map<string, string> res;
			res.insert(make_pair(sub_ref, alignment));
			sub_r.insert(make_pair(len_q,res));
		}
		else{
			unordered_map<string, string> res = sub_r[len_q];
			if(res.find(sub_ref) == res.end()){
				res.insert(make_pair(sub_ref, alignment));
				sub_r.insert(make_pair(len_q,res));
			}
		}
	}
}

// perform global alignment using ksw2. Return the alignment result
string global_aligment(const char *tseq, const char *qseq, int sc_mch, int sc_mis, int gapo, int gape){
	void *km = 0;
	// int8_t q = 4, e = 2;
	int w = -1, rep = 1;
	ksw_extz_t ez;
	int i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis;
	int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
	int tl = strlen(tseq), ql = strlen(qseq);
	uint8_t *ts, *qs, c[256];
	memset(&ez, 0, sizeof(ksw_extz_t));
	memset(c, 4, 256);
	c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
	c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; 
	ts = (uint8_t*)malloc(tl);
	qs = (uint8_t*)malloc(ql);
	for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; 
	for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
	ez.score = ksw_gg(km, ql, qs, tl, ts, 5, mat, gapo, gape, w, &ez.m_cigar, &ez.n_cigar, &ez.cigar);
	string rs = "";
	for (i = 0; i < ez.n_cigar; ++i) 
		rs = rs + to_string(ez.cigar[i]>>4) + "MID"[ez.cigar[i]&0xf];
	free(ez.cigar); free(ts); free(qs);
	return rs;
}

// perform global and extension alignment. Return alignment result as a string
string extension_alignment(const char *tseq, const char *qseq, int sc_mch, int sc_mis, int gapo, int gape){
	int i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; 
	int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
	int tl = strlen(tseq), ql = strlen(qseq);
	uint8_t *ts, *qs, c[256];
	ksw_extz_t ez;

	memset(&ez, 0, sizeof(ksw_extz_t));
	memset(c, 4, 256);
	c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
	c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; 
	ts = (uint8_t*)malloc(tl);
	qs = (uint8_t*)malloc(ql);
	for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; 
	for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
	ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);
	string rs = "";
	for (i = 0; i < ez.n_cigar; ++i) 
		rs = rs + to_string(ez.cigar[i]>>4) + "MID"[ez.cigar[i]&0xf];
	free(ez.cigar); free(ts); free(qs);
	return rs;
}

// reverse the order of the alignment result
string reverse_align(string align){
	regex rg("M|I|D");
	vector<int> index;
	for(auto it = std::sregex_iterator(align.begin(), align.end(), rg); it != std::sregex_iterator(); ++it)
	{
		index.push_back(it->position());
	}
	string rs = "";
	for(int i = index.size()-2; i>=0; i--){
		rs += align.substr(index[i]+1, index[i+1]-index[i]);
	}
	rs += align.substr(0, index[0]+1);
	return rs;
}

//perform alignment on the whole query and reference by appending results of subproblems.
string alignment(string ref, string q, list<int>ref_index, list<int>q_index, list<int>match_len, int sc_mch, int sc_mis, int gapo, int gape){
	string result = "";

	//do alignment for the most left subsequences
	int begin_match_ref= ref_index.front();
	int begin_match_q = q_index.front();
	
	if(begin_match_q == 0 && begin_match_ref != 0 ){
		result += to_string(begin_match_ref) +"D";
	}
	if(begin_match_q != 0 && begin_match_ref == 0){
		result += to_string(begin_match_q) +"I";
	}
	if(begin_match_q != 0 && begin_match_ref != 0){
		string begin_ref = ref.substr(0,begin_match_ref);
		string begin_q = q.substr(0,begin_match_q);
		string align = cache_finding(0, begin_match_q, begin_ref);

		if( align == ""){
			reverse(begin_ref.begin(), begin_ref.end());
			reverse(begin_q.begin(), begin_q.end());
			string align_rs = extension_alignment(begin_ref.c_str(), begin_q.c_str(), sc_mch, sc_mis, gapo, gape);
			align_rs = reverse_align(align_rs);
			result += align_rs;
			add_to_cache(0, begin_match_q,begin_ref, align_rs);
		}
		else{
			result += align;
		}
		
	}

	//do global alignment for subsequences between exact matches
	std::vector<int> ref_vector{begin(ref_index), end(ref_index)};
	std::vector<int> q_vector{begin(q_index), end(q_index)};
	std::vector<int> len_vector{begin(match_len), end(match_len)};
	for(int i = 0; i< q_index.size()-1; ++i){
		int off_q = q_vector[i] + len_vector[i];
		int len_q = q_vector[i+1]-off_q;
		string sub_q = q.substr(off_q, len_q);
		string sub_ref = ref.substr(ref_vector[i]+len_vector[i], ref_vector[i+1]- ref_vector[i]-len_vector[i]);
		result += to_string(len_vector[i]) + "M";
		string find_cache_rs = cache_finding(off_q, len_q, sub_ref);
		if(find_cache_rs != ""){
			result += find_cache_rs;
		}
		else{
			string glob_align = global_aligment(sub_ref.c_str(), sub_q.c_str(), sc_mch, sc_mis, gapo, gape);
			add_to_cache(off_q, len_q, sub_ref, glob_align);
			result += glob_align;
		}

	}

	//do extension alignment for the most right subsequences
	int end_match_ref = ref_index.back();
	int end_match_q = q_index.back();
	int end_match_len = match_len.back(); 
	if(end_match_q + end_match_len == q.length() && end_match_ref + end_match_len < ref.length()){
		result += to_string(end_match_len)+"M";
		result += to_string(ref.length()- end_match_ref- end_match_len) +"D";
	}
	else if(end_match_q + end_match_len < q.length() && end_match_ref + end_match_len == ref.length()){
		result += to_string(end_match_len)+"M";
		result += to_string(q.length()- end_match_q- end_match_len) +"I"; 
	}
	else if (end_match_q + end_match_len == q.length() && end_match_ref + end_match_len == ref.length()){
		result += to_string(end_match_len)+"M";
	}
	else{
		result += to_string(end_match_len)+"M";
		end_match_ref = end_match_ref + end_match_len;
		end_match_q = end_match_q + end_match_len;
		int end_len_q = q.length() - end_match_q;
		string end_sub_ref = ref.substr(end_match_ref, ref.length()-end_match_ref);
		string end_align = cache_finding(end_match_q, end_len_q, end_sub_ref);
		if(end_align == ""){
			string end_sub_q = q.substr(end_match_q, end_len_q);
			string end_align_rs = extension_alignment(end_sub_ref.c_str(), end_sub_q.c_str(), sc_mch, sc_mis, gapo, gape);
			result += end_align_rs;
		}
		else{
			result += end_align;
		}
		// result += end_align;



	}
	return result;
}

//read fasta file
string read_file(string filename){
	ifstream file(filename);
	string line, id, genome;
	while(getline(file, line)){
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

//read files contain exact match results
vector<vector<string>> read_exact_match(string filename){
	ifstream file(filename);
	std::vector<vector<string>> match_pair;
	string line, id;
	while(getline(file, line)){
		// cout << line << endl;
		if(line.empty()){
			continue;
		}
		if(line[0] == '>'){
			id = line.substr(1);
		}
		else{
			line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1");
			istringstream buf(line);
    		istream_iterator<string> beg(buf), end;
   			vector<string> tokens(beg, end);
			match_pair.push_back(tokens);
		}
	}
	file.close();
	return match_pair;
}  

int main(int argc, char*argv[]) 
{ 	
	
	int number_ref = atoi(argv[1]);
	string time_compare = "./test_ref/time_compare.txt";
	ofstream file(time_compare);
	for(int i = 0; i<number_ref; ++i){
		string query = read_file("./test_ref/test/query.fa");
		string ref_dir = "./test_ref/test/ref_" + to_string(i)+".fa";
		string ref = read_file(ref_dir);
		list<int> ref_index ;
		list<int> q_index ;
		list<int> len ;
		string index_dir = "./test_ref/test/result/result_"+to_string(i)+".txt";
		vector<vector<string>> matches = read_exact_match(index_dir);
		if(matches.size() > 1){
			cout<<i<<endl;
			for(unsigned j = 1; j<matches.size(); j = j+2){
				// cout << matches.at(j)[0] <<" "<< matches.at(j)[1]<< " "<<matches.at(j)[2]<<endl;
				ref_index.push_back(stoi(matches.at(j)[0])-1);
				q_index.push_back(stoi(matches.at(j)[1])-1);
				len.push_back(stoi(matches.at(j)[2]));
			}

			cout<< ref.substr(210, 36)<<endl;
			cout<<query.substr(209,36)<<endl;
			cout <<ref.substr(246, 39)<<endl;
			cout << query.substr(245,39)<<endl;
			auto start = high_resolution_clock::now(); 
			cout<<'*'<<alignment(ref, query, ref_index, q_index, len, 1, -2, 2, 1).c_str()<<endl;
			auto stop = high_resolution_clock::now(); 
			auto duration = duration_cast<microseconds>(stop - start);
			// cout<<'*'<< duration.count()<<endl;

			auto start2 = high_resolution_clock::now(); 
			cout<<'-'<<global_aligment(ref.c_str(), query.c_str(), 1, -2, 2, 1).c_str()<<endl;
			auto stop2 = high_resolution_clock::now(); 
			auto duration2 = duration_cast<microseconds>(stop2 - start2);
			// cout<< '-'<<duration2.count()<<endl;

			
			file << duration.count()<<" "<< duration2.count()<<endl;
			
		}
		
	}
	file.close();	
		

	
	// cout << "all"<<query.substr(108,4) << endl;
	
	// cout <<"all ref "<< ref.substr(684,4)<<endl;
	// list<int> ref_index {9375, 9446, 9535, 9574, 9598, 9629, 9676, 9727, 9776};
	// list<int> q_index {47, 118, 207, 246, 270, 301, 348, 399, 448};
	// list<int> len {28, 28, 26, 23, 30, 33, 29, 48, 20};
	// cout<<alignment(ref, query, ref_index, q_index, len, 1, -2, 2, 1).c_str()<<endl;
	// string s1 = "atATGCTAcCGCat";
	// string s2 = "AT";
	// printf("%s",extension_alignment(s1.c_str(),s2.c_str(), 1, -2, 2, 1).c_str());
	
	// cout << q_index.size()<<","<<ref_index.size()<<","<<len.size()<< endl;

	// cout<<alignment(ref, query, ref_index, q_index, len, 1, -2, 2, 1).c_str()<<endl;
	
	return 0;

} 