
#include <stdio.h>
#include <string>
#include <map>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include "PredicateMap.h"
#include "bt9_reader.h"


bt9::BT9Reader* PredicateMap::flowGraph;
std::unordered_map<uint32_t, cond_struct_t> PredicateMap::predicated_blocks;
  

void PredicateMap::setFlowGraph(bt9::BT9Reader* new_flow_graph) {
  // assert (flowGraph == nullptr);
  flowGraph = new_flow_graph;
}

cond_struct_t PredicateMap::getDefault(){
  cond_struct_t result;
  result.is_default = true;
  printf("getDefault\n");
  return result;
}

cond_struct_t PredicateMap::getPredicate(uint32_t key) {
  if (predicated_blocks.find(key) != predicated_blocks.end())
    return predicated_blocks.find(key)->second;
  
  if (flowGraph == nullptr){
    printf("no flow graph -- ");
    return getDefault();
  }

  // from the given branch instruction
  bt9::BT9Reader::NodeTableIterator initial_node(flowGraph, key);
  bt9::BT9ReaderEdgeRecord *taken_path = initial_node.getNextEdge(true);
  bt9::BT9ReaderEdgeRecord *not_taken_path = initial_node.getNextEdge(false);

  if (!initial_node->brClassConditionalityIs("CND") 
      || taken_path == nullptr 
      || not_taken_path == nullptr
      || taken_path->brVirtualTarget() < initial_node->brVirtualAddr()) // backward branch
  {
    // printf("unvalid intital -- ");
    // return getDefault();
    cond_struct_t result;
    result.is_default = true;
    // printf("%d %d %d %d %08lx\n", !initial_node->brClassConditionalityIs("CND"), taken_path == nullptr ,
    //   not_taken_path == nullptr, initial_node.getBranchTarget() < initial_node->brVirtualAddr(), initial_node->brVirtualAddr());
    return result;
  }
  bt9::BT9Reader::NodeTableIterator nt_end_node(flowGraph, not_taken_path->destNodeIndex());
  bt9::BT9Reader::NodeTableIterator  t_end_node(flowGraph, taken_path->destNodeIndex());


  // then-block browsing
  while (nt_end_node->brVirtualAddr() < t_end_node->brVirtualAddr()) {
    if (nt_end_node->brClassConditionalityIs("UCD")) {
      if (nt_end_node.getNextEdge(true)->brVirtualTarget() < taken_path->brVirtualTarget()) {
        // can't build a correct representation
        return getDefault();
      }
      
      if (nt_end_node->brVirtualAddr() +4 < taken_path->brVirtualTarget()) {
        // can't build a correct representation
        return getDefault();
      }

      // collect data
      printf("move on (nt) (UCD): %u -> %u\n", nt_end_node->brNodeIndex(), nt_end_node.getNextEdge(true)->destNodeIndex());
      // move forward the node
      nt_end_node += nt_end_node.getNextEdge(true)->destNodeIndex() - nt_end_node->brNodeIndex();
    }

    else if (nt_end_node->brClassConditionalityIs("UCD")) {
      // search for predicated construction
      cond_struct_t path = getPredicate(nt_end_node->brNodeIndex()); // also checks if forward 

      if (path.is_default) { // no result
        return getDefault();
      }

      // collect data
      
      printf("move on (nt) (CND): %u -> %u\n", nt_end_node->brNodeIndex(), nt_end_node.getNextEdge(true)->destNodeIndex());
      // move forward the node
      nt_end_node += path.dest_node - nt_end_node->brNodeIndex();
    }

    else {
      // conditionnality error
      return getDefault();
    }
  }


  printf("got nt to %u (%0lx)\n", nt_end_node->brNodeIndex(), nt_end_node->brVirtualAddr());
  
  
  // now that not taken flow is further than the taken one
  // browse the taken path until find a complex structure
  while (t_end_node->brVirtualAddr() < nt_end_node->brVirtualAddr()) {
    if (t_end_node->brClassConditionalityIs("UCD")) {
      // no favorable issue
      return getDefault();
    }

    else if (t_end_node->brClassConditionalityIs("CND")) {
      cond_struct_t path = getPredicate(t_end_node->brNodeIndex()); // also checks if forward 

      if (path.is_default) { // no result
        return getDefault();
      }

      // collect data
      printf("move on  (t) (CND): %u -> %u\n", t_end_node->brNodeIndex(), t_end_node.getNextEdge(true)->destNodeIndex());
      
      // move forward the node
      t_end_node += path.dest_node - t_end_node->brNodeIndex();
    }
    
    else {
      // conditionnality error
      return getDefault();
    }
  }

  if (t_end_node->brVirtualAddr() == nt_end_node->brVirtualAddr()) {
    // found a valid structure

    cond_struct_t result;

    result.last_pc_included = taken_path->brVirtualTarget();
    result.instruction_number = not_taken_path->nonBrInstCnt() + 1;
    result.taken_path = taken_path->edgeIndex();
    result.not_taken_path = not_taken_path->edgeIndex();
    result.dest_node = t_end_node->brNodeIndex();

    predicated_blocks[key] = result;
    return result;

  }
 



  return getDefault();
}
















  // while (t_end_node->brNodeIndex() != nt_end_node->brNodeIndex()) {
    
  //   bt9::BT9Reader::NodeTableIterator* earlier_node = &nt_end_node;
  //   if (nt_end_node->brVirtualAddr() > t_end_node->brVirtualAddr())
  //     earlier_node = &t_end_node;
    
  //   if (nt_end_node->brClassConditionalityIs("UCD")) {
  //     if (nt_end_node.getNextEdge(true)->brVirtualTarget() < t_end_node->brVirtualAddr()) {
  //       // can't build a correct representation
  //       return getDefault();
  //     }
    
  //     (*earlier_node) = bt9::BT9Reader::NodeTableIterator(flowGraph, earlier_node->getNextEdge(true)->destNodeIndex());

  //   }

  //   if ((*earlier_node)->brClassConditionalityIs("CND")) {
  //     cond_struct_t path = getPredicate(earlier_node->getNextEdge(true)->destNodeIndex());

      

  //   }


  // }






  // // try to find a if_then
  // if (taken_path->destNodeIndex() == not_taken_path->destNodeIndex()) 
  // {
  //   // return an ift block
  //   cond_struct_t result;

  //   result.last_pc_included = taken_path->brVirtualTarget();
  //   result.instruction_number = not_taken_path->nonBrInstCnt() + 1;
  //   result.taken_path = taken_path->edgeIndex();
  //   result.not_taken_path = not_taken_path->edgeIndex();
  //   result.dest_node = taken_path->destNodeIndex();

  //   predicated_blocks[key] = result;
  //   return result;
  // }


  // // or if_then_else structure
  // bt9::BT9Reader::NodeTableIterator nt_dest(flowGraph, not_taken_path->destNodeIndex());
  // bt9::BT9ReaderEdgeRecord *over_taken_path = nt_dest.getNextEdge(true);
  // if (nt_dest->brClassConditionalityIs("UCD") 
  //     && over_taken_path->destNodeIndex() == taken_path->destNodeIndex()) 
  // {
  //   // return an ifte block
  //   cond_struct_t result;

  //   result.last_pc_included = over_taken_path->brVirtualTarget();
  //   result.instruction_number = not_taken_path->nonBrInstCnt() + taken_path->nonBrInstCnt() + 1;
  //   result.taken_path = taken_path->edgeIndex();
  //   result.not_taken_path = not_taken_path->edgeIndex();
  //   result.dest_node = taken_path->destNodeIndex();

  //   predicated_blocks[key] = result;
  //   return result;
  // }



  // // 