#ifndef UTILITIES.H
#define UTILITIES.H

#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <math.h>

using std::string;
using std::cin;
using std::cout;
using std::endl;
using std::vector;


class Utilities {
    public:
        template <typename T> void print_vector(const T& v) const {
            // Print a vector of any type :)
            for (const auto& elem: v) {
                cout << elem << " ";
            }
            cout << endl;
        }

        template <typename T> void print_matrix(const T& m) const {
            // Print a matrix of any type :)
            for (const auto& vec : m ) {
                print_vector(vec);
            }
        }

        template <typename T> T requestInput(char type, T min, T max, string m) const {
            T input;
            
            bool condition_satisfied = false;
            bool type_satisfied = false;
            
            double temp;
            
            while (!condition_satisfied || !type_satisfied) {
                cout << m << endl;
                
                if(type == 'i') {
                    cin >> temp;
                    if (!cin.fail()) {
                        
                        if ((temp - floor(temp) != 0)) type_satisfied = false;
                        else { 
                            type_satisfied = true;
                            input = static_cast<int>(temp); 
                        }
                        
                    }
                    else type_satisfied = false;
                } else {
                    cin >> input;  
                    
                    type_satisfied = !cin.fail(); 
                }

                if (input < min || input > max) condition_satisfied = false;
                else condition_satisfied = true;    

                cin.clear();
                cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            
            return input;
        }

};

#endif