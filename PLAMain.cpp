/* This is a C++ program to print a Programmable Logic Array or PLA Table.
 * The functions are minimized through the Quine-McCluskey Algorithm in SOP or Sum Of Products form.
 * The program does not account for 'don't care' conditions.
*/

/* MAIN PROGRAM FILE */

//Including the required libraries and header files
#include <iostream>
#include "PLA.h"

//Including the required entities/identifiers
using std::cout;
using std::cin;
using std::endl;

int main()
{
    //Declaring program constants
    const int func = 10;
    const int minterm = 200;
    const int numstring = 500;
    
    //Displaying menu
    cout<<"-------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Welcome!"<<endl<<"This is a C++ program to print a Programmable Logic Array or PLA Table."<<endl;
    cout<<"The functions are minimized through the Quine-McCluskey Algorithm in SOP or Sum Of Products form."<<endl;
    cout<<"The program does not account for 'don't care' conditions."<<endl;
    cout<<"-------------------------------------------------------------------------------------------------"<<endl;
    
    char option;
    
    //Calling the main overall function
    do
    {
        inputs_and_mainfunctioncalls(func, minterm, numstring);
        
        cout<<"Do you want to quit(Q/q) or print another PLA(P/p)? ";
        cin>>option;
        
    } while(option == 'P' || option == 'p');
    
    cout<<"Thank you."<<endl;
    cout<<endl;
    return 0;
    
}
