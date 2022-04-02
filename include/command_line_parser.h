#ifndef COMMAND_LINE_PARSER_H
#define COMMAND_LINE_PARSER_H
#include "options.h"
#include <utility>
namespace multipole_conv {
std::pair<std::size_t, MPOptions> cmd_parser(int argc, char *argv[]);
}
#endif
