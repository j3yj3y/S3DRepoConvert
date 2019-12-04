#include <iostream>
#include <string>
#include <algorithm>
#include <filesystem>

#include "range_scan_reader.h"

template < typename T >
T convert( const std::string &str )
{
    std::istringstream ss(str);
    T num;
    ss >> num;
    return num;
}

bool optionSet( char **start, char **end, const std::string &optionName )
{
    return std::find(start, end, optionName) != end;
}

template < typename T >
T getOption( char **start, char **end, const std::string &optionName, const std::string &abbreviation, T defaultValue )
{
    char **itr = std::find(start, end, optionName);

    if (itr != end && ++itr != end)
        return convert<T>(*itr);

    if (!abbreviation.empty())
    {
        itr = std::find(start, end, abbreviation);
        if (itr != end && ++itr != end)
            return convert<T>(*itr);
    }
    return defaultValue;
}

void print_usage()
{
    std::cout << "Usage: " << std::endl;
    std::cout << "\t" << "stanford3drepoconvertert <input-path> <output-path>" << std::endl << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "\t" << "-i <input-path> \tInput 'conf' file containing range-scan" << std::endl;
    std::cout << "\t" << "-o <output-path> \tResulting output 'ply' file containing transformed and aligned scans" << std::endl;
    std::cout << "\t" << "--help prints out this usage page" << std::endl;
}



int main(int argc, char** argv)
{
    //if (argc < 3 || optionSet(argv, argv + argc, "--help")) {
    //    print_usage();
    //    return EXIT_SUCCESS;
    //}
    //
    //const std::string conf_file   = getOption(argv, argv + argc, "-i", "", std::string(""));
    //const std::string output_file = getOption(argv, argv + argc, "-o", "", std::string(""));

    const std::filesystem::path conf_file{"/data/datasets/Stanford3dScanRepository/happy_buddha/happy_scans/happy_back/happyBackRight.conf"};
    stanford_repo::range_scan::transform(conf_file, "./test.ply");

    return EXIT_SUCCESS;
}



