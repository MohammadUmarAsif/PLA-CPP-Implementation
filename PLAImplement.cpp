/* PROGRAM FUNCTION IMPLEMENTATION FILE */

//Including the required libraries and header files
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "PLA.h"

//Including the required entities/identifiers
using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::string;
using std::sort;
using std::ostringstream;
using std::setw;

/* Function Implementations */

//Function to obtain the complement of the entered minterm function
void obtaining_complements(vector<vector<int>> &mintermvector, int number_of_functions, int max_min, vector<vector<int>> &mintermvectorcomplement)
{
    
    for(int i=0; i < number_of_functions; i++)
    {
        for(int j=0; j < max_min; j++)
        {
            int count=0;
            
            //Checking if minterm is already in the normal function
            for(unsigned int e=0; e < mintermvector[i].size(); e++)
            {
                if(j == mintermvector[i][e])
                    count++;
            }
            
            //If not, then it is added to the complement function
            if(count==0)
                mintermvectorcomplement[i].push_back(j);
        }
    }
    
    
}

//Function to obtain the minimized minterm function
void function_minimization(vector<vector<int>> &mintermvector, int func, int minterm, int numstring, int number_of_functions, 
                            int number_of_variables, vector<vector<string>> &minimizedbinary)
{
    
    vector<vector<int>> sorted (func, vector<int>(minterm,0));
    
    //Clearing the declared vector
    clearing_2dvector_int(sorted);
    
    //Sorting the minterms in the function by the number of 1's
    sorting_vector(mintermvector,sorted, func, minterm, number_of_functions, number_of_variables);
    
    vector<vector<string>> binary (func, vector<string>(minterm, ""));
    
    //Clearing the declared vector
    clearing_2dvector(binary);
    
    //Converting the minterms to their binary strings e.g. (1 -> 0001)
    for(int i=0; i<number_of_functions; i++)
    {
        for(unsigned int j=0; j<sorted[i].size(); j++)
        {
            binary[i].push_back(decimal_to_binary(sorted[i][j],number_of_variables));
        }
    }
    
    vector<vector<string>> stringminterm (func, vector<string>(minterm, ""));
    
    //Clearing the declared vector
    clearing_2dvector(stringminterm);
    
    //Converting the minterms to their corresponding number strings e.g. (1(int) -> 1(string))
    for(int i=0; i<number_of_functions; i++)
    {
        for(unsigned int j=0; j<sorted[i].size(); j++)
        {
            stringminterm[i].push_back(to_string<int>(sorted[i][j]));
        }
    }
    
    vector<vector<string>> tablerecord (func, vector<string>(minterm, ""));
    vector<vector<string>> historyrecord (func, vector<string>(minterm, ""));
    vector<vector<string>> tablerecordleftover (func, vector<string>(minterm, ""));
    vector<vector<string>> historyrecordleftover (func, vector<string>(minterm, ""));
    
    //Clearing the declared vectors
    clearing_2dvector(tablerecord);
    clearing_2dvector(historyrecord);
    clearing_2dvector(tablerecordleftover);
    clearing_2dvector(historyrecordleftover);
    
    //Inserting the first column of binary numbers
    for(int i = 0;i<number_of_functions; i++)
        for(unsigned int j = 0; j<sorted[i].size(); j++)
            tablerecord[i].push_back(binary[i][j]);
       
    //Inserting the first column of number strings
    for(int i = 0;i<number_of_functions; i++)
        for(unsigned int j = 0; j<sorted[i].size(); j++)
            historyrecord[i].push_back(stringminterm[i][j]);

    //Performing the minimization through Tabular Method or Quine-McCluskey Algorithm
    performing_tabulation(number_of_functions, number_of_variables, tablerecord, historyrecord, tablerecordleftover, historyrecordleftover, binary, stringminterm);
    
    //Removing any duplicate terms in the vectors
    removing_duplicates(tablerecordleftover,number_of_functions, func, minterm);
    removing_duplicates(historyrecordleftover,number_of_functions, func, minterm);
    
    //Removing the common terms so as to get only unticked terms in the tablerecord and historyrecord
    removing_common_terms(tablerecord,tablerecordleftover, number_of_functions);
    removing_common_terms(historyrecord,historyrecordleftover, number_of_functions);
 
    vector<vector<string>> essentialprimeimplicants (func, vector<string>(minterm, ""));
    vector<vector<string>> essentialprimehistory (func, vector<string>(minterm, ""));
    vector<vector<vector<string>>> stringhistory (func, vector<vector<string>>(minterm, vector<string>(numstring, "")));
    vector<vector<vector<string>>> stringhistoryessential (func, vector<vector<string>>(minterm, vector<string>(numstring, "")));
    vector<vector<vector<string>>> stringhistorynonessential (func, vector<vector<string>>(minterm, vector<string>(numstring, "")));

    //Clearing the declared vectors
    clearing_2dvector(essentialprimeimplicants);
    clearing_2dvector(essentialprimehistory);
    
    //Clearing the declared vectors
    clearing_3dvector(stringhistory);
    clearing_3dvector(stringhistoryessential);
    clearing_3dvector(stringhistorynonessential);
        
    //Obtaining the minterms of each combination as individual strings
    seperating_strings(stringhistory,historyrecord,number_of_functions);

    //Obtaining the essential terms
    obtaining_essential_terms(number_of_functions, sorted, stringhistory, stringminterm, essentialprimeimplicants, essentialprimehistory, tablerecord, historyrecord);
    
    //Obtaining the minterms of each combination of essential minterms as individual strings
    seperating_strings(stringhistoryessential,essentialprimehistory,number_of_functions);

    vector<vector<string>> nonessential (func, vector<string>(minterm, ""));
    vector<vector<string>> nonessentialhistory (func, vector<string>(minterm, ""));
    vector<vector<string>> notcovered (func, vector<string>(minterm, ""));

    //Clearing the declared vectors
    clearing_2dvector(nonessential);
    clearing_2dvector(nonessentialhistory);
    clearing_2dvector(notcovered);
    
    //Removing the common terms so as to get non-essential terms in the tablerecord and historyrecord
    removing_common_terms(tablerecord,essentialprimeimplicants,number_of_functions);
    removing_common_terms(historyrecord,essentialprimehistory,number_of_functions);
 
    for(int i=0; i<number_of_functions; i++)
    for(unsigned int j=0 ; j<tablerecord[i].size(); j++)
        nonessential[i].push_back(tablerecord[i][j]);

    for(int i=0; i<number_of_functions; i++)
    for(unsigned int j=0 ; j<historyrecord[i].size(); j++)
        nonessentialhistory[i].push_back(historyrecord[i][j]);

    //Obtaining the minterms of each combination of non-essential minterms as individual strings
    seperating_strings(stringhistorynonessential,nonessentialhistory,number_of_functions);

    //Obtaining the minterms that are not covered by essential terms
    obtaining_notcovered_minterms(number_of_functions, sorted, stringhistoryessential, stringminterm, notcovered);

    vector<vector<string>> finalnonessential (func, vector<string>(minterm, ""));
    vector<vector<int>> countmeter (func, vector<int>(minterm, 0));
    
    //Clearing the declared vectors
    clearing_2dvector(finalnonessential);
    clearing_2dvector_int(countmeter);
    
    //Selecting the best combination of the non-essential minterms to ensure minimized function
    minimizing_nonessential_minterms(number_of_functions, stringhistorynonessential, notcovered, countmeter, 
                                        nonessential, finalnonessential, func, minterm, numstring);
    
    for (int i=0; i<number_of_functions; i++)
        for(unsigned int j = 0 ; j<finalnonessential[i].size(); j++)
            essentialprimeimplicants[i].push_back(finalnonessential[i][j]);
      
    //Allocating the final terms (essential + non-essential) in one vector
    for(int i=0 ; i<number_of_functions; i++)
        for(unsigned j=0; j<essentialprimeimplicants[i].size(); j++)
            minimizedbinary[i].push_back(essentialprimeimplicants[i][j]);
            
            
}

//Function to remove the duplicate terms in a vector
void removing_duplicates(vector<vector<string>> &vectorfunction, int number_of_functions, int func, int minterm)
{
    //Declaring temporary array
    vector<vector<string>> temparray(func, vector<string>(minterm,""));
    clearing_2dvector(temparray);
    
    vector<int> arraysize (func,0);
    arraysize.clear();
    
    for(int i=0; i<number_of_functions; i++)
    {   
        arraysize[i] = 0;
    }
    
    //Copying elements to the temporary array
    for(int i=0; i<number_of_functions; i++)
    {
        for(unsigned int j=0; j<vectorfunction[i].size(); j++)
        {
            temparray[i].push_back(vectorfunction[i][j]);
            arraysize[i]++;
        }
    }
    
    //Removing the duplicate terms in the temporary array
    for(int l=0; l<number_of_functions ; l++)
    {
        for(int i=0;i<arraysize[l];++i)
        {
            for(int j=i+1;j<arraysize[l];)
            {
                if(temparray[l][i]==temparray[l][j])
                {
                    for(int k=j;k<arraysize[l]-1;++k)
                        temparray[l][k]=temparray[l][k+1];

                    --arraysize[l];
                }
                else
                    ++j;
            }
        }
    }   

    //Clearing the original vector
    for(int i = 0;i<number_of_functions; i++)
    {
        vectorfunction[i].clear();
    }
    
    //Copying the unique elements into the original vector
    for(int i=0; i<number_of_functions; i++)
    {
        for(int j=0; j<arraysize[i]; j++)
        {
            vectorfunction[i].push_back(temparray[i][j]);
        }
    }
    
    
}
  
//Function to remove those terms in a vector that it has in common with another vector
void removing_common_terms(vector<vector<string>> &bigvector,vector<vector<string>> &smallvector, int number_of_functions)
{
    //Creating a vector to store the indexes of the common terms
    vector<vector<int>> delete_index (10,vector<int> (64,0));
    for(int i = 0;i<number_of_functions; i++)
        delete_index[i].clear();
    
    //Storing the indexes of the common terms
    for(int i=0; i<number_of_functions; i++)
    {
        for(unsigned int j=0; j<bigvector[i].size(); j++)
        {
            for(unsigned int k=0; k<smallvector[i].size(); k++)
            {
                if(bigvector[i][j] == smallvector[i][k])
                {
                    delete_index[i].push_back(j);
                }
            }
        }
    }
    
    //Erasing the common terms from the vector
    for(int i=0; i<number_of_functions; i++)
    {
        for(unsigned int j=0; j<delete_index[i].size(); j++)
        {
             bigvector[i].erase(bigvector[i].begin() + (delete_index[i][j] - j));
        }
    }
    
}

//Function to seperate the minterm combinations of a term into individual strings
void seperating_strings(vector<vector<vector<string>>> &individualstring, vector<vector<string>> &stringfunction, int number_of_functions)
{
    //Seperating the string combination through the commas e.g. (8,16) -> 8 and 16 (as individual strings)
    for(int i=0; i<number_of_functions; i++)
    {
        for(unsigned int j=0; j<stringfunction[i].size(); j++)
        {
            int count=0, comma =0;
            for(unsigned int k=0; k<stringfunction[i][j].size(); k++)
            {
                //Keep increasing till comma appears, and then take the substring within the given indexes
                if(stringfunction[i][j][k] != ',')
                    count++;
                if(stringfunction[i][j][k] == ',')
                {  
                    comma++;
                    if(comma == 1)
                    {
                        individualstring[i][j].push_back(stringfunction[i][j].substr(0,count));
                    }
                    else
                    {
                        individualstring[i][j].push_back(stringfunction[i][j].substr(k-count,count));
                    }
                    count=0;
                        
                }
                if(k==stringfunction[i][j].size()-1)
                {
                    individualstring[i][j].push_back(stringfunction[i][j].substr(k-count+1,count));
                    count=0;
                }
            }
        }
    }
    
    
}

//Function to count the number of 1's in a binary number using bitwise operators
unsigned int count_ones(unsigned int number)
{ 

    unsigned int count = 0; 
    
    while (number) 
    { 
        count += number & 1; 
        number >>= 1; 
    } 
    
    return count; 

} 

//Function to convert an element into a string using a template
template <class T>
string to_string(T s)
{

    ostringstream oss;
    
    oss << s;
    
    return oss.str();

}

//Function to return the letters based on the number of variables
vector<string> variable_letters(int number_of_variables)
{
   
    vector<string> v;
    string letters[]={"A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"};
    
    //Storing the required number of letters in a vector and returning it
    for(int i=0;i<number_of_variables;i++)
        v.push_back(letters[i]);
    
    return v;

}

//Function to represent a binary term through the corresponding variables
string convert_to_letters(string a, int number_of_variables)
{

    string temp="";
    vector<string> vars=variable_letters(number_of_variables);
    
    //If 0, then A'. If 1, then A, If - then -;
    for(unsigned int i=0;i<a.length();i++)
    {
        if(a[i]!='-')
        {
            if(a[i]=='0')
                temp=temp+vars[i]+"\'";
            else
                temp=temp+vars[i];
        }
    }
    
    return temp;

}

//Function to convert a decimal number into its binary form
string decimal_to_binary(int number, int max_digits)
{
    
    string temp;
    
    for(int i=0;i<max_digits;i++)
    {
        
        if(number%2==1)
            temp = temp + "1";
        else 
            temp = temp + "0";
            
        number = number/2;
    }
    
    //Reversing the string e.g. 1 is converted to 1000 -> 0001
    reverse(temp.begin(), temp.end()); 
    
    return temp;
    
}

//Function to combine the terms that differ by one bit at index p
string combiningstrings(string a, string b, int p, int number_of_variables)
{
    
    string temp;
    
    //If index is not p, keep copying. When index is p, replace by a -
    for(int i=0; i<number_of_variables; i++)
    {
        if(i!=p)
            temp = temp + a[i];
        else
            temp = temp + "-";
    }
    
    return temp;
    
}

//Function to combine the minterm combination (history) of the terms that differ by one bit
string combiningrecordstrings(string a, string b, int count)
{
    
    string temp;
    
    temp = a + "," + b ;
    
    return temp;
    
}

//Function to clear a 2d vector of strings
void clearing_2dvector(vector<vector<string>> &vector2d)
{
    
    for(unsigned int i = 0;i<vector2d.size(); i++)
        vector2d[i].clear();

}

//Function to clear a 2d vector of integers
void clearing_2dvector_int(vector<vector<int>> &vector2d)
{
    
    for(unsigned int i = 0;i<vector2d.size(); i++)
        vector2d[i].clear();
        
}

//Function to clear a 3d vector of strings
void clearing_3dvector(vector<vector<vector<string>>> &vector3d)
{
    
    for(unsigned int i = 0;i<vector3d.size(); i++)
        for(unsigned int j = 0;j<vector3d[i].size(); j++)
            vector3d[i][j].clear();
            
}

//Function to sort the minterms based on the number of 1's in the binary form
void sorting_vector(vector<vector<int>> &tosort, vector<vector<int>> &sorted, int func, int minterm, int number_of_functions, int number_of_variables)
{
    
    int dataentry[func][minterm];
    
    //Storing the number of 1's in each term of a function into an array
    for(int i=0; i<number_of_functions; i++)
    {
        for(unsigned int j=0 ; j<tosort[i].size(); j++)
        {
            dataentry[i][j] = count_ones(tosort[i][j]);
        }
    }
    
    //Sorting the array based on the number of 1's in each term using the sort function
    for(int i=0; i<number_of_functions; i++)
    {
        sort(dataentry[i],dataentry[i] + tosort[i].size());
    }
    
    //Storing the sorted elements into a new vector
    for(int i=0; i<number_of_functions; i++)
    {
        for(int k=0 ; k<=number_of_variables; k++)
        {
            for(unsigned int j=0 ; j<tosort[i].size(); j++)
            {
                if(static_cast<signed>(count_ones(tosort[i][j])) == k)
                {
                    sorted[i].push_back(tosort[i][j]);
                }
            }
        }
    }
    
    
}

//Function to implement the Tabular Method and obtain the prime implicants
void performing_tabulation(int number_of_functions, int number_of_variables, 
                            vector<vector<string>> &tablerecord, vector<vector<string>> &historyrecord, 
                            vector<vector<string>> &tablerecordleftover, vector<vector<string>> &historyrecordleftover, 
                            vector<vector<string>> &binary, vector<vector<string>> &stringminterm)
{

    for(int i=0; i<number_of_functions; i++)
    {
        for(int count =0;count<number_of_variables; count++)
        {  
            int size = tablerecord[i].size();
            for(int j=0; j<size-1; j++)
            {   
                for(int l=j+1; l<size; l++)
                {
                    int flag = 0,p;
                    
                    //For the first column of the table
                    if (count == 0 )
                    {
                        //Checking if the terms differ by one bit and storing the index
                        for(int k=0; k<number_of_variables; k++)
                        {
                            if(binary[i][j][k] != binary[i][l][k])
                            {    
                                flag++;
                                p = k;
                            }
                        }
                        if(flag == 1)   
                        {
                            //Storing the ticked terms and their corresponding minterm combinations (history)
                            tablerecordleftover.at(i).push_back(binary[i][j]);
                            tablerecordleftover.at(i).push_back(binary[i][l]);
                            historyrecordleftover.at(i).push_back(stringminterm[i][j]);
                            historyrecordleftover.at(i).push_back(stringminterm[i][l]);
                        
                            string combined,combinedrecord;
                            
                            //Combining and storing the terms and their corresponding minterm combinations (history)
                            combined = combiningstrings(binary[i][j],binary[i][l],p,number_of_variables);
                            tablerecord.at(i).push_back(combined); 
                            combinedrecord = combiningrecordstrings(stringminterm[i][j], stringminterm[i][l], count);        
                            historyrecord.at(i).push_back(combinedrecord);
                        
                        }
                        
                    }
                    
                    //For the subsequent columns of the table
                    if(count>0)
                    {    
                        flag=0;
                        int ctr= 0;
                        
                        //Checking if the dashes coincide
                        for(int k=0; k<number_of_variables; k++)
                        {
                            if(tablerecord[i][j][k] == '-' && tablerecord[i][l][k] == '-')
                            {
                                ctr++; 
                            }
                        }
                        
                        if(ctr==count)
                        {   
                            //Checking if the terms differ by one bit and storing the index
                            for(int k=0; k<number_of_variables; k++)
                            {
                                if(tablerecord[i][j][k] =='0' && tablerecord[i][l][k] == '1')
                                {    
                                    flag++;
                                    p = k;
                                }
                                else if(tablerecord[i][j][k] == '1' && tablerecord[i][l][k] == '0')
                                {    
                                    flag++;
                                    p = k;
                                }   
                            }
                        }
                    
                        if( flag == 1)   
                        {
                            string combined,combinedrecord;
                            int c2=0;
                            
                            combined = combiningstrings(tablerecord[i][j],tablerecord[i][l],p,number_of_variables);
                            
                            //Storing the ticked terms and their corresponding minterm combinations (history)
                            tablerecordleftover.at(i).push_back(tablerecord[i][j]);     
                            tablerecordleftover.at(i).push_back(tablerecord[i][l]);
                            historyrecordleftover.at(i).push_back(historyrecord[i][j]);
                            historyrecordleftover.at(i).push_back(historyrecord[i][l]);
                            
                            //Checking if the combination has already been entered into the table
                            for(unsigned int e=0; e<tablerecord[i].size(); e++)
                            {
                                if(tablerecord[i][e] == combined)
                                    c2++;
                            }
                            
                            if(c2==0)
                            {
                                //Combining and storing the terms and their corresponding minterm combinations (history)
                                tablerecord.at(i).push_back(combined);
                                combinedrecord = combiningrecordstrings(historyrecord[i][j], historyrecord[i][l], count);
                                historyrecord.at(i).push_back(combinedrecord);

                            }      
                        }
                    }
                }
            }
        }
    }
    
                        
}

//Function to obtain the essential prime implicants and their history
void obtaining_essential_terms(int number_of_functions, vector<vector<int>> &sorted, vector<vector<vector<string>>> &stringhistory, 
                            vector<vector<string>> &stringminterm, vector<vector<string>> &essentialprimeimplicants,
                            vector<vector<string>> &essentialprimehistory, vector<vector<string>> &tablerecord, vector<vector<string>> &historyrecord)
{
    
    for(int i=0; i<number_of_functions; i++)
    {
        for(unsigned int l=0; l<sorted[i].size(); l++)
        {
            int count=0;
            int p;
            //Checking if a minterm is present in only one combined term's history
            for(unsigned int j=0; j<stringhistory[i].size(); j++)
            {
                for(unsigned int k=0; k<stringhistory[i][j].size(); k++)
                {
                    //If yes, then store the index of the combined term in p
                    if(stringhistory[i][j][k] == stringminterm[i][l])
                        {
                            count++;
                            p = j;
                        }
                }
            }
            
            //Pushing the implicants and their history to their respective vectors
            if(count == 1)
            {
                int c2=0;
                for(unsigned int e=0; e<essentialprimeimplicants[i].size(); e++)
                {
                    //Checking if the essential prime implicant has already been entered
                    if(essentialprimeimplicants[i][e] == tablerecord[i][p])
                        c2++;
                }   
                            
                if(c2==0)
                {
                    essentialprimeimplicants[i].push_back(tablerecord[i][p]);
                    essentialprimehistory[i].push_back(historyrecord[i][p]);
                }
            }
        }
    }
    
                                
}

//Function to obtain the minterms that are not covered by the essential prime implicants
void obtaining_notcovered_minterms(int number_of_functions, vector<vector<int>> &sorted, vector<vector<vector<string>>> &stringhistoryessential, 
                            vector<vector<string>> &stringminterm, vector<vector<string>> &notcovered)
{
       
    for(int i=0; i<number_of_functions; i++)
    {
        for(unsigned int l=0; l<sorted[i].size(); l++)
        {
            int c2=0;
            
            //Checking if the minterm is in the history of the essential prime implicants
            for(unsigned int u=0 ; u<stringhistoryessential[i].size(); u++)
            {
                for(unsigned int e=0; e<stringhistoryessential[i][u].size(); e++)
                {
                    if(stringminterm[i][l] == stringhistoryessential[i][u][e])
                        c2++;
                }
            }
            
            //If not, then it is added to the vector
            if(c2==0)
            {
                notcovered[i].push_back(stringminterm[i][l]);
            }
            
        }
    }
 
                               
}

//Function to obtain the least number of non-essential prime implicants that cover the remaining minterms
void minimizing_nonessential_minterms(int number_of_functions, vector<vector<vector<string>>> &stringhistorynonessential, 
                            vector<vector<string>> &notcovered, vector<vector<int>> &countmeter, vector<vector<string>> &nonessential, 
                            vector<vector<string>> &finalnonessential, int func, int minterm, int numstring)
{
    //Creating a vector to check if the non-essential prime implicants are covering the same terms
    vector<vector<vector<string>>> checking(func, vector<vector<string>>(minterm, vector<string>(numstring, "")));
    
    //Clearing the declared vector
    clearing_3dvector(checking);
    
    //Keeping a record of how many leftover minterms each non-essential prime implicant covers
    for(int i=0; i<number_of_functions; i++)
    {
        for(unsigned int j=0 ; j <stringhistorynonessential[i].size(); j++)
        {
            int count=0;
            for(unsigned int k=0; k<stringhistorynonessential[i][j].size(); k++)
            {
                for(unsigned int l=0; l<notcovered[i].size(); l++)
                {
                    if(stringhistorynonessential[i][j][k] == notcovered[i][l])
                    {
                        count++;
                        checking[i][j].push_back(notcovered[i][l]);
                        
                    }
                    
                    if(j>0)
                    {
                        int ctr=0;
                        
                        //Checking if the minterms covered by the current implicant are the same as that of the previous implicant's
                        for(unsigned int q=0; q<countmeter[i].size() ; q++)
                        {   
                            if(count == countmeter[i][q])
                            {
                                for(unsigned int u = 0; u<checking[i][j].size(); u++)
                                {
                                    for(unsigned int y =0; y<checking[i][q].size() ; y++)
                                    {    
                                        if(checking[i][j][u] == checking[i][q][y])
                                            ctr++;
                                    }
                                }
                                
                                //If same minterms are covered
                                if(ctr==count)
                                {
                                    //Giving least priority and jumping to the next iteration
                                    countmeter[i].push_back(0);
                                    goto out;
                                }
                            }
                        }
                    }
                    
                }
            }
            
            //Storing the number of minterms a non-essential prime implicant covers
            countmeter[i].push_back(count);
            out:
            ;
        }
    }

    //Obtaining the least number of non-essential prime implicants
    for(int i=0; i<number_of_functions; i++)
    {
        unsigned int count=0; 
        
        //Giving highest priority to the minterm with maximum value of count
        for(int c=0; c<*max_element(countmeter[i].begin(), countmeter[i].end()); c++)
        {      
            for(unsigned int j=0 ; j <nonessential[i].size(); j++)
            {
                if(countmeter[i][j] == static_cast<signed>(*max_element(countmeter[i].begin(), countmeter[i].end())) - c)
                {
                    finalnonessential[i].push_back(nonessential[i][j]);
                    
                    //Checking if all the leftover minterms have been covered
                    for(unsigned int k=0; k<stringhistorynonessential[i][j].size(); k++)
                    {
                        for(unsigned int l=0; l<notcovered[i].size(); l++)
                        {
                            if(stringhistorynonessential[i][j][k] == notcovered[i][l])
                            {
                                count++;
                            }
                        }
                    }
                }
                
                //If yes, then no need to check for remaining implicants
                if(count == notcovered[i].size())
                {
                     goto outside;
                }
            }
        }
    }
    
    outside:
    ;
                                
}

//Function to return the function names based on the number of functions
vector<string> function_letters(int number_of_functions)
{
    
    vector<string> v;
    string function[]={"F1","F2","F3","F4","F5","F6","F7","F8","F9","F10"};
    string complement[] = {"F1'","F2'","F3'","F4'","F5'","F6'","F7'","F8'","F9'","F10'"};
    
    //Storing the required number of function names in a vector and returning it
    for(int i=0;i<number_of_functions;i++)
        v.push_back(function[i]);
    
    for(int i=0;i<number_of_functions;i++)
        v.push_back(complement[i]);
        
    return v;

}

//Function to call the combination function 
void making_combinations(vector<string> &funcstring, int n, int num, vector<vector<int>> &combinations, vector<vector<string>> &funcname)  
{  
    //Declaring a vector to store the names of the functions involved in a combination
    vector<string> data (num,"");
    
    //Clearing the declared vector
    data.clear();
    
    //Declaring an array to store the indexes of the functions involved in a combination
    int dataindex[10];
    int o=0;
    
    //Calling the combination function
    combination_function(funcstring, data, 0, n-1, 0, num, dataindex, combinations, o, funcname);  
}  

//Function to make all possible combinations of the minterm functions and their complements
void combination_function(vector<string> &funcstring, vector<string> &data, int start, int end, int index, int num, int dataindex[], 
                            vector<vector<int>> &combinations, int &o, vector<vector<string>> &funcname)  
{ 
    //If all the functions are present in the combination, store the required combination
    if (index == num)  
    {  
        //Storing indexes of the functions
        for(int j = 0; j < num; j++)  
            combinations[o].push_back(dataindex[j]);  
            
        //Storing the function names
        for(int j = 0; j < num; j++)
            funcname[o].push_back(data[j]);
        
        //Incrementing for the next combination
        o++;
        
        return;  
    }  
    
    //Recursively making the combinations
    for (int i = start; i <= end &&  end - i + 1 >= num - index; i++)  
    {  
        dataindex[index]= i;
        
        //Condition to prevent a functions complement from being combined with it e.g. F2 and F2' cannot be in the same combination
        for(int j=0; j<num; j++)
        {
            //If encountered, then skip the iteration
            if(dataindex[j] + num == i )
                goto out;
        }
        
        data[index] = funcstring[i];  
        
        //Calling the combination function again 
        combination_function(funcstring, data, i+1,  
                        end, index+1, num,dataindex, combinations, o, funcname);
        out:
        ;
    } 
    

}  

//Main Function to find the best combination of all functions
void best_function_combination(vector<vector<string>> &minimizedbinary, vector<vector<string>> &minimizedbinarycomplement, 
                                int number_of_functions, int func, int minterm, vector<string> &finalterms, vector<string> &finalfunctions, 
                                vector<vector<string>> &functionrecord, int number_of_variables)
{
                                    
    vector<vector<string>> allfunctions;
    
    //Clearing the declared the vector
    clearing_2dvector(allfunctions);
    
    //Storing all the functions in a single vector
    for(int i=0; i<number_of_functions; i++)
    {
        allfunctions.push_back(minimizedbinary[i]);
    }
    for(int i=0; i<number_of_functions; i++)
    {
        allfunctions.push_back(minimizedbinarycomplement[i]);
    }

    int x = pow(2,func);
    int y = minterm*3;
         
    //Storing function names in a vector e.g. F1 and F2
    vector<string> funcstring = function_letters(number_of_functions);
    
    vector<vector<int>> combinations(x,vector<int>(func,0));
    vector<vector<string>> funcname(x,vector<string>(func,""));
    
    //Clearing the declared vectors
    clearing_2dvector_int(combinations);
    clearing_2dvector(funcname);
    
    int num = number_of_functions;
    int n = funcstring.size(); 
    
    //Calling the function to make all combinations
    making_combinations(funcstring, n, num, combinations, funcname);  

    vector<vector<string>> allcombinations(x, vector<string>(y,""));
    vector<vector<string>> allfunctionsletters(2*func, vector<string>(minterm,""));
    
    //Clearing the declared vectors
    clearing_2dvector(allcombinations);
    clearing_2dvector(allfunctionsletters);

    //Combining the implicants based on the results of the function combination
    for(unsigned int i=0; i<combinations.size(); i++)
    {
        for(unsigned int j=0; j<combinations[i].size(); j++)
        {
            for(unsigned int k=0; k<allfunctions[combinations[i][j]].size(); k++)
                allcombinations[i].push_back(allfunctions[combinations[i][j]][k]);
        }
    }

    //Removing any duplicate implicants
    removing_duplicates(allcombinations,pow(2,number_of_functions),x, y);
    
    //Finding the combination with the least number of implicants
    smallest_combination(allcombinations, funcname, number_of_functions, finalterms, finalfunctions, functionrecord, allfunctions, funcstring);     

    //Converting the binary forms to letters
    for(unsigned int i=0; i<allfunctions.size(); i++)
    {
        for(unsigned int j=0; j<allfunctions[i].size(); j++)
        {
            allfunctionsletters[i].push_back(convert_to_letters(allfunctions[i][j],number_of_variables));
        }
    }
    
    cout<<endl;
    cout<<endl;
    cout<<endl;
    
    //Printing all the functions and their names
    cout<<"____________________________________________________________"<<endl;
    cout<<"The minimized functions are:"<<endl;
     for(unsigned int i=0; i<allfunctions.size(); i++)
    {
        cout<<funcstring[i]<<" = ";
        for(unsigned int j=0; j<allfunctions[i].size(); j++)
        {
            if(j+1 == allfunctions[i].size())
                cout<<allfunctionsletters[i][j];
            else
                cout<<allfunctionsletters[i][j]<<" + ";
        }
        cout<<endl;
    }
    cout<<"____________________________________________________________"<<endl;
    
    

}

//Function to find the combination with the least number of implicants
void smallest_combination(vector<vector<string>> &allcombinations, vector<vector<string>> &funcname, int number_of_functions, 
                            vector<string> &finalterms, vector<string> &finalfunctions, vector<vector<string>> &functionrecord, 
                            vector<vector<string>> &allfunctions, vector<string> &funcstring)
{

    vector<string> smallest = allcombinations[0];
    int position=0;
    
    //Finding smallest combination
    for(int i=1; i<pow(2,number_of_functions); i++)
    {
        if(smallest.size() > allcombinations[i].size())
        {
            smallest = allcombinations[i];
            position = i;
        }   
    }
    
    //Storing the final implicants
    for(unsigned int i=0; i<allcombinations[position].size(); i++)
    {
        finalterms.push_back(allcombinations[position][i]);
    }
    
    //Storing the final function names
    for(unsigned int i=0; i<funcname[position].size(); i++)
    {
        finalfunctions.push_back(funcname[position][i]);
    }
    
    //Keeping a record of which implicants are in which functions
    for(unsigned int i=0; i<finalterms.size(); i++)
    {
        for(unsigned int j=0; j<allfunctions.size(); j++)
        {
            for(unsigned int k=0; k<allfunctions[j].size(); k++)
            {
                if(allfunctions[j][k] == finalterms[i])
                    functionrecord[i].push_back(funcstring[j]);
            }
        }
    }
    
                            
}

//Function to print the Progammable Logic Array or PLA
void printing_PLA(int number_of_functions, int number_of_variables, vector<string> &finalfunctions, vector<string> &finalterms, 
                    vector<vector<string>> &functionrecord, vector<string> &finaltermsletters)
{
    //Creating a vector to store the required letters
    vector<string> varletters=variable_letters(number_of_variables);
    
    int lineamount = number_of_variables+14+8+(number_of_variables+1)*2+16+(number_of_variables+1)*2+16+10;
    
    //Displaying the main heading
    cout<<setw(lineamount/3+25)<<"PLA PROGRAMMING TABLE\n";
    
    for(int i=0; i<lineamount+7; i++)
        cout<<"_";
    cout<<endl;
    
    //Displaying the first row of sub-headings
    cout<<setw(number_of_variables+14)<<"Product Terms"<<setw(8)<<"S.No.";
    cout<<setw((number_of_variables+1)*2+16)<<"Input Variables"<<setw((number_of_variables+1)*2+16)<<"Output Functions";
    cout<<endl;
    
    for(int i=0; i<lineamount+7; i++)
        cout<<"_";
    cout<<endl;
    
    //Displaying the second row of sub-headings
    for(int i=0; i<number_of_variables+14+8; i++)
        cout<<" ";
    cout<<setw(7)<<" ";
    
    for(int j=0; j<number_of_variables; j++)
        cout<<setw(5)<<varletters[j];
    cout<<setw(8)<<" ";
    
    for(unsigned int i=0; i<finalfunctions.size(); i++)
        cout<<setw(5)<<finalfunctions[i];
    cout<<setw(8)<<" ";
    cout<<endl;
    
    for(int i=0; i<lineamount+7; i++)
        cout<<"_";
    cout<<endl;
    
    //Displaying the calculated data
    for(unsigned int i=0; i<finaltermsletters.size(); i++)
    {
        string count="-";
        cout<<setw(number_of_variables+14)<<finaltermsletters[i]<<setw(8)<<i+1;
        cout<<setw(7)<<" ";
        for(unsigned int j=0; j<finalterms[i].size(); j++)
            {
                cout<<setw(5)<<finalterms[i][j];
            }
        cout<<setw(8)<<" ";
        for(unsigned int k=0; k<finalfunctions.size(); k++)
        {
            for(unsigned int o=0; o<functionrecord[i].size(); o++)
            {
                if(finalfunctions[k] == functionrecord[i][o])
                    count = "1";
            }
            cout<<setw(5)<<count;
            count = "-";
        }
        cout<<setw(8)<<" ";
        cout<<endl;
    }
    
    for(int i=0; i<lineamount+7; i++)
        cout<<"_";
    cout<<endl;
    
                            
}

//Function to take the inputs and call the main functions
void inputs_and_mainfunctioncalls(int func, int minterm, int numstring)
{
    int number_of_functions, number_of_minterms, number_of_variables;
    long long max_min, max_func;
    
    //Setting limits and taking the inputs
    again1:
    cout<<"Enter the number of variables (1-26): ";
    cin>>number_of_variables;
    if(number_of_variables==0 || number_of_variables>26)
        {
            cout<<"PLA Table cannot be formed. Please enter a value between 1 and 26."<<endl<<endl;
            goto again1;
        }

    max_func = pow(2,pow(2,number_of_variables));
    max_min = pow(2,number_of_variables);

    again2:
    cout<<"Enter the number of functions (1 to "<<max_func<<"): ";
    cin>>number_of_functions;
    if(number_of_functions==0 || number_of_functions>max_func)
        {
            cout<<"PLA Table cannot be formed. Please enter a value between 1 and "<<max_func<<"."<<endl<<endl;
            goto again2;
        } 
        
    vector<vector<int>> mintermvector (func, vector<int>(minterm,0));
    
    //Clearing the declared vector
    clearing_2dvector_int(mintermvector);
    
    for(int i=0; i < number_of_functions; i++)
    {
        again3:
        cout<<"Enter the number of minterms for function "<<i+1<<" (1 - "<<max_min<<"): ";
        cin>>number_of_minterms;
        
        if(number_of_minterms==0 || number_of_minterms>max_min)
        {
            cout<<"PLA Table cannot be formed. Please enter a value between 1 and "<<max_min<<"."<<endl<<endl;
            goto again3;
        }
        cout<<"Enter the minterms for function "<<i+1<<",seperated by spaces (0 to "<<max_min-1<<"): ";
        
        for(int j=0; j < number_of_minterms; j++)
        {   
            int input;
            again4:
            cin>>input;
            if(input<0 || input>max_min)
            {
                cout<<"Minterms must lie in the range of (0 to "<<max_min-1<<")."<<endl<<endl;
                goto again4;
            }
            mintermvector[i].push_back(input);
        }

    }
    
    cout<<endl;
    cout<<"Thank you for the inputs.";
    
    vector<vector<int>> mintermvectorcomplement (func, vector<int>(minterm,0));
    
    //Clearing the declared vector
    clearing_2dvector_int(mintermvectorcomplement);
    
    //Calling the function to obtain the complement functions
    obtaining_complements(mintermvector, number_of_functions, max_min, mintermvectorcomplement);
  
    vector<vector<string>> minimizedbinary (func, vector<string>(minterm, ""));
    vector<vector<string>> minimizedbinarycomplement (func, vector<string>(minterm, ""));
    
    //Clearing the declared vectors
    clearing_2dvector(minimizedbinary);
    clearing_2dvector(minimizedbinarycomplement);
   
    //Minimizing the functions and their complements
    function_minimization(mintermvector, func, minterm, numstring, number_of_functions, number_of_variables, minimizedbinary);
    function_minimization(mintermvectorcomplement, func, minterm, numstring, number_of_functions, number_of_variables, minimizedbinarycomplement);
    
    vector<string> finalterms (func*minterm, "");
    vector<string> finalfunctions (func, "");
    vector<vector<string>> functionrecord(func*minterm, vector<string>(func, ""));
    
    //Clearing the declared vectors
    finalterms.clear();
    finalfunctions.clear();
    clearing_2dvector(functionrecord);
    
    //Finding the function combination with the least number of implicants
    best_function_combination(minimizedbinary, minimizedbinarycomplement, number_of_functions, func, minterm, finalterms, 
                                finalfunctions, functionrecord, number_of_variables);
    
    vector<string> finaltermsletters (func*minterm, "");
    
    //Clearing the declared vector
    finaltermsletters.clear();
    
    for(unsigned int i=0; i<finalterms.size(); i++)
        finaltermsletters.push_back(convert_to_letters(finalterms[i], number_of_variables));
    
    //Displaying the final function combination and implicants
    cout<<"The best combination of functions is: ";
    for(unsigned int i=0; i<finalfunctions.size(); i++)
        cout<<finalfunctions[i]<<"  ";
    cout<<endl;
        
    cout<<"The final minterms are: ";
    for(unsigned int i=0; i<finaltermsletters.size(); i++)
    {
        if( i+1 == finaltermsletters.size() )
            cout<<finaltermsletters[i];
        else
            cout<<finaltermsletters[i]<<" + ";
            
    }
    cout<<endl;
    cout<<"____________________________________________________________"<<endl;
    
    cout<<endl;
    cout<<endl;
    cout<<endl;

    //Displaying the PLA
    printing_PLA(number_of_functions, number_of_variables, finalfunctions, finalterms, functionrecord, finaltermsletters);
    
    cout<<endl;
    cout<<endl;

    
}