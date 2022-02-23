#ifndef __UTILITIES_H_
#define __UTILITIES_H_

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

// Utility functions to be used throughout the code
class Utilities {
    public:
        template <typename T> void print_vector(const T& v) const {
            // @ Print a vector of any type :)
            for (const auto& elem: v) {
                cout << elem << " ";
            }
            cout << endl;
        }

        template <typename T> void print_matrix(const T& m) const {
            // @ Print a matrix of any type :)
            for (const auto& vec : m ) {
                print_vector(vec);
            }
        }

        template <typename T> T requestInput(char type, T min, T max, string m) const {
            // @ request input from the user with some restrictions. will keep asking for the same input
            //   unless the entered value is correct (i.e. satifies the restrictions)

            T input;
            
            // conditions to keep asking for the input 
            bool condition_satisfied = false;
            bool type_satisfied = false;
            
            double temp;
            
            while (!condition_satisfied || !type_satisfied) {
                // output the message of the required input
                cout << m;
                
                /* == Integer case == */
                if(type == 'i') {
                    cin >> temp;
                    if (!cin.fail()) {
                        
                        if ((temp - floor(temp)) != 0) type_satisfied = false;
                        else { 
                            type_satisfied = true;
                            input = static_cast<int>(temp); 
                        }
                        
                    }
                    else type_satisfied = false;
                } 
                /* == Other case (string or double) == */
                else {
                    cin >> input; 
                    type_satisfied = !cin.fail();
                }

                if (input < min || input > max) condition_satisfied = false;
                else condition_satisfied = true;    

                // clear current input
                cin.clear();
                //clear buffer before taking a new line
                cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            
            return input;
        }

};

#endif