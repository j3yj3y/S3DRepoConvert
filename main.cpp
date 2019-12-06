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
    std::cout << "\t" << "-i <input-path> \tInput directory containing stanford conf files (subdirectories are allowed too)" << std::endl;
    std::cout << "\t" << "-o <output-path> \tOutput directory containing transformed and aligned scans" << std::endl;
    std::cout << "\t" << "--help prints out this usage page" << std::endl;
}



int main(int argc, char** argv)
{
    if (argc < 3 || optionSet(argv, argv + argc, "--help")) {
        print_usage();
        return EXIT_SUCCESS;
    }

    const std::filesystem::path input_dir  = getOption(argv, argv + argc, "-i", "", std::string(""));
    const std::filesystem::path output_dir = getOption(argv, argv + argc, "-o", "", std::string(""));

    if (!std::filesystem::exists(output_dir))
        std::filesystem::create_directory(output_dir);

    std::vector<std::filesystem::path> conf_files;
    for (const auto& p : std::filesystem::recursive_directory_iterator(input_dir))
        if(p.path().extension() == ".conf")
            conf_files.push_back(p);

    for (auto& conf : conf_files)
        stanford_repo::range_scan::transform(conf, output_dir / conf.filename().replace_extension(".ply"));

    return EXIT_SUCCESS;
}



