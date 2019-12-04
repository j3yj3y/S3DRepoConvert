#ifndef RANGE_SCAN_READER_H
#define RANGE_SCAN_READER_H

#include <filesystem>


namespace stanford_repo { namespace range_scan {

bool transform( const std::filesystem::path& conf_file,
                const std::filesystem::path& output );


} // end namespace range_scan
} // end namespace stanford_repo


#endif // RANGE_SCAN_READER_H

