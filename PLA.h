/* PROGRAM FUNCTION HEADER FILE */

//Using include guard
#ifndef PLA_H
#define PLA_H

//Including the required libraries and header files
#include <string>
#include <vector>

//Including the required entities/identifiers
using std::vector;
using std::string;


//Function prototypes
void obtaining_complements(vector<vector<int>> &mintermvector, int number_of_functions, int max_min, vector<vector<int>> &mintermvectorcomplement);
void function_minimization(vector<vector<int>> &mintermvector, int func, int minterm, int numstring, int number_of_functions, 
                            int number_of_variables, vector<vector<string>> &minimizedbinary);
void removing_duplicates(vector<vector<string>> &vectorfunction, int number_of_functions, int func, int minterm);  
void removing_common_terms(vector<vector<string>> &bigvector,vector<vector<string>> &smallvector, int number_of_functions);
void seperating_strings(vector<vector<vector<string>>> &individualstring, vector<vector<string>> &stringfunction, int number_of_functions);
unsigned int count_ones(unsigned int number);
template <class T>
string to_string(T s);
vector<string> variable_letters(int number_of_variables);
string convert_to_letters(string a, int number_of_variables);
string decimal_to_binary(int number, int max_digits);
string combiningstrings(string a, string b, int p, int number_of_variables);
string combiningrecordstrings(string a, string b, int count);
void clearing_2dvector(vector<vector<string>> &vector2d);
void clearing_2dvector_int(vector<vector<int>> &vector2d);
void clearing_3dvector(vector<vector<vector<string>>> &vector3d);
void sorting_vector(vector<vector<int>> &tosort, vector<vector<int>> &sorted, int func, int minterm, int number_of_functions, int number_of_variables);
void performing_tabulation(int number_of_functions, int number_of_variables, 
                            vector<vector<string>> &tablerecord, vector<vector<string>> &historyrecord, 
                            vector<vector<string>> &tablerecordleftover, vector<vector<string>> &historyrecordleftover, 
                            vector<vector<string>> &binary, vector<vector<string>> &stringminterm);
void obtaining_essential_terms(int number_of_functions, vector<vector<int>> &sorted, vector<vector<vector<string>>> &stringhistory, 
                                vector<vector<string>> &stringminterm, vector<vector<string>> &essentialprimeimplicants,
                                vector<vector<string>> &essentialprimehistory, vector<vector<string>> &tablerecord, vector<vector<string>> &historyrecord);
void obtaining_notcovered_minterms(int number_of_functions, vector<vector<int>> &sorted, vector<vector<vector<string>>> &stringhistoryessential, 
                                    vector<vector<string>> &stringminterm, vector<vector<string>> &notcovered);
void minimizing_nonessential_minterms(int number_of_functions, vector<vector<vector<string>>> &stringhistorynonessential, 
                                        vector<vector<string>> &notcovered, vector<vector<int>> &countmeter, vector<vector<string>> &nonessential, 
                                        vector<vector<string>> &finalnonessential, int func, int minterm, int numstring);
vector<string> function_letters(int number_of_functions);
void making_combinations(vector<string> &funcstring, int n, int num, vector<vector<int>> &combinations, vector<vector<string>> &funcname);
void combination_function(vector<string> &funcstring, vector<string> &data, int start, int end, int index, int num, int dataindex[], 
                            vector<vector<int>> &combinations, int &o, vector<vector<string>> &funcname);  
void best_function_combination(vector<vector<string>> &minimizedbinary, vector<vector<string>> &minimizedbinarycomplement, 
                                int number_of_functions, int func, int minterm, vector<string> &finalterms, vector<string> &finalfunctions, 
                                vector<vector<string>> &functionrecord, int number_of_variables);
void smallest_combination(vector<vector<string>> &allcombinations, vector<vector<string>> &funcname, int number_of_functions, 
                            vector<string> &finalterms, vector<string> &finalfunctions, vector<vector<string>> &functionrecord, 
                            vector<vector<string>> &allfunctions, vector<string> &funcstring);
void printing_PLA(int number_of_functions, int number_of_variables, vector<string> &finalfunctions, vector<string> &finalterms, 
                    vector<vector<string>> &functionrecord, vector<string> &finaltermsletters);
void inputs_and_mainfunctioncalls(int func, int minterm, int numstring);

#endif
