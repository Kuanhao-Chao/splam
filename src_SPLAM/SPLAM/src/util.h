#include <string>
#include <iostream>

std::string get_full_path(std::string fname){
    const char *cur_path = fname.c_str();
    char *actualpath;


    actualpath = realpath(cur_path, NULL);
    if (actualpath != NULL){
        return std::string(actualpath);
        free(actualpath);
    }
    else{
        std::cerr<<"could not resolve path: "<<fname<<std::endl;
        exit(-1);
    }
}