#include <algorithm>
#include <iostream>
#include <vector>
#include <utility>
#include <regex>
#include <string>
#include <filesystem>
#include <unordered_map>

class DisjointSet{
    std::unordered_map<int, int> root;
    //constructor
public:
    void initialize(const std::vector<int> & universe){
        for(int i: universe){
            root[i]=i;
        }
    }

    int find_root(int k){
        while(root[k]!=k){
            k = root[k];
        }
        return k;
    }

    void make_union(int i, int j){
        int x = find_root(i);
        int y = find_root(j);

        root[x]=y;
    }
};

std::vector<int> get_vertices(const std::vector<std::pair<int,std::pair<int, int> > > frame_graph)
{
    std::vector<int> universe;
    std::unordered_map<int, int> set;
    for(auto edg: frame_graph)
    {
        if(set.find(edg.second.first)== set.end())
        {
            universe.push_back(edg.second.first);
            set[edg.second.first] = edg.first;
        }
        if(set.find(edg.second.second)== set.end())
        {
            universe.push_back(edg.second.second);
            set[edg.second.second] = edg.first;
        }
    }    
    return universe;
}

std::vector<std::pair<int, std::pair<int, int> > > Kruskal_MST(const std::vector<std::pair<int, std::pair<int, int> > > & frame_graph)
{
    std::vector<std::pair<int, std::pair<int, int> > > frame_tree;
    std::vector<std::pair<int, std::pair<int, int> > > frame_graph_copy = frame_graph;
    std::sort(frame_graph_copy.begin(), frame_graph_copy.end(),
        [](const std::pair<int, std::pair<int, int> > & a, const std::pair<int, std::pair<int, int> > & b)
        {
            return a.first < b.first;
        });
    std::vector<int> universe = get_vertices(frame_graph_copy);
    DisjointSet DS;
    DS.initialize(universe);
    for(auto edg: frame_graph_copy)
    {
        int x = edg.second.first;
        int y = edg.second.second;
        if(DS.find_root(x)!= DS.find_root(y))
        {
            DS.make_union(x,y);
            frame_tree.push_back(edg);
        }
    }
    return frame_tree;
}

int main()
{
    std::vector<std::pair<int, std::pair<int,int> > > frame_graph;
    frame_graph.push_back(std::make_pair(1, std::make_pair(1, 2)));
    frame_graph.push_back(std::make_pair(2, std::make_pair(2, 3)));
    frame_graph.push_back(std::make_pair(3, std::make_pair(4, 3)));
    frame_graph.push_back(std::make_pair(4, std::make_pair(1, 4)));
    frame_graph.push_back(std::make_pair(5, std::make_pair(1, 5)));
    frame_graph.push_back(std::make_pair(6, std::make_pair(6, 2)));
    frame_graph.push_back(std::make_pair(7, std::make_pair(5, 7)));
    frame_graph.push_back(std::make_pair(8, std::make_pair(6, 8)));
    frame_graph.push_back(std::make_pair(9, std::make_pair(3, 8)));
    frame_graph.push_back(std::make_pair(10, std::make_pair(7, 8)));
    frame_graph.push_back(std::make_pair(11, std::make_pair(7, 4)));
    frame_graph.push_back(std::make_pair(12, std::make_pair(5, 6)));
    std::vector<std::pair<int, std::pair<int, int> > > ff = Kruskal_MST(frame_graph);
    for(auto edg: ff)
    {
        std::cout<< edg.first << ":" << edg.second.first <<"," << edg.second.second << std::endl;
    }
    std::regex r(".*\\.obj");
    std::regex exclude(".*(test|good\\.mesh__sf)\\.obj");
    std::string data_root = "/Users/vector_cat/gits/bcclean_jupyters/data/2";
    std::string good_mesh_file, bad_mesh_file;
    for (const auto & entry : std::filesystem::directory_iterator(data_root))
    {
        if(std::regex_match(entry.path().c_str(), r))
        {
            if(!std::regex_match(entry.path().c_str(), exclude))
            {
                bad_mesh_file = entry.path().c_str();
                std::printf("%s\n",bad_mesh_file.c_str());
            }
            else
            {
                std::printf("wrong%s\n",entry.path().c_str());
            }
            
        }
    }
    
    return 0;
}