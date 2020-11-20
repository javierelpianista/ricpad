#include <ricpad.hpp>
#include <fstream>
#include <boost/algorithm/string.hpp>

mpfr_float E0;

void read_input(string filename) {
    map<string, string> parsed_options;
    map<string, string>::iterator it;

    ifstream read_file(filename);
    string line;
    vector<string> tokens;

    while ( getline( read_file, line ) ) {
        boost::split(tokens, line, boost::is_any_of(" "));

        auto p = tokens.begin();
        string what = *p++;
        string value;

        if ( what == "D" ) {
            Dmin = stoi(*p++);
            Dmax = stoi(*p);
        } else if ( what == "d" ) {
            d = stoi(*p);
        } else if ( what == "E0" ) {
            E0 = mpfr_float(*p);
        } else {
            while ( p != tokens.end() ) {
                value += *p++;
            }

            parsed_options[what] = value;
        }
    }

    it = parsed_options.find("var");
    if ( it == parsed_options.end() ) {
        cout << "The variable should be defined with var <VAR>" << endl;
        throw 10;
    } 
    
    x = gi::symbol(it->second);
    gi::symtab table;
    table[it->second] = x;
    gi::parser reader(table);

    it = parsed_options.find("pot");
    if ( it == parsed_options.end() ) {
        cout << "Potential should be defined with pot <POT>" << endl;
        throw 10;
    }

    potential = reader(it->second);
}
