#ifndef FINAL_ALIGN
#define FINAL_ALIGN

#include <stdio.h>
#include <stdlib.h>
#include <string>
std::string generate_ref (std::string query, int max_len_insert_string, int max_insert_position);
void write_ref(int number_ref, std::string query, int max_len_insert_string, int max_insert_position);
std::string read_query(std::string filename);
void write_query(std::string query);
#endif