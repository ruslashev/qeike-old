#include "load_world.hh"
#include "../engine/utils.hh"

void load_world(const std::string &filename) {
  std::vector<char> file;
  qke::utils::read_file_to_vector(filename, file);
}

