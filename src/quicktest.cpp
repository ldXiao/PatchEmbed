//
// Created by Lind Xiao on 7/30/19.
//
#include <map>
#include <vector>
#include <iostream>
int main(){
    std::cout << "ahha"<< std::endl;
    std::map<int, std::vector<int> > dict;
    dict[0].push_back(1);
    dict[0].push_back(2);
    dict[1].push_back(3);
}
