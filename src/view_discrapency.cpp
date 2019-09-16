#include <iostream>
#include "edge.h"
int main(){
    std::cout << "discrapency" <<std::endl;
    bcclean::pair_map<int> a; 
    a[std::pair<int,int>(1,2)]=1;
    
    return 0;
}