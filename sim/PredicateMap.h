#ifndef HEADER_CONDITIONAL_STRUCTURE
#define HEADER_CONDITIONAL_STRUCTURE


#include <stdio.h>
#include <string>
#include <map>
#include <list>
#include <unordered_map>
#include <unordered_set>
//#include "bt9_reader.h"

namespace bt9 {
  class BT9Reader;
}

struct cond_struct_t
{
  // cond_struct_t(uint64_t last_pc_included, uint32_t instruction_number, uint32_t taken_path, uint32_t not_taken_path, 
  // uint32_t dest_node, )
  uint64_t last_pc_included = 0; // included in predicated block
  uint32_t instruction_number = 0; // total instructions in the control flow path
  uint32_t taken_path = 0;
  uint32_t not_taken_path = 0;
  uint32_t dest_node = 0;
  bool is_taken_predicated = false;
  bool is_not_taken_predicated = false;
  bool is_default = false;
};

class PredicateMap
{
public:
  static cond_struct_t getPredicate(uint32_t key);
  static void setFlowGraph(bt9::BT9Reader* new_flow_graph);

private:
  static cond_struct_t getDefault();

  static bt9::BT9Reader* flowGraph;
  static std::unordered_map<uint32_t, cond_struct_t> predicated_blocks;
  
};

#endif
// cond_struct_t PredicateMap::default_struct{.is_default = true};