///////////////////////////////////////////////////////////////////////
//  Copyright 2015 Samsung Austin Semiconductor, LLC.                //
///////////////////////////////////////////////////////////////////////

// Description : Main file for CBP2016

#include <assert.h>
#include <stdlib.h>
#include <string.h>
// #include <map>
#include <list>
#include <unordered_map>
#include <map>
using namespace std;

#include "utils.h"
//#include "bt9.h"
#include "bt9_reader.h"
//#include "predictor.cc"
#include "predictor.h"

#define COUNTER unsigned long long
//#define NB_PRE_PRED 6

// HOTSPOT STUDY ////////////////////////////////////////////////
#define HOTSPOT_KEY_SIZE 16
#define HOTSPOT_TAG_SIZE 10
typedef uint16_t hotspot_t;
#define HOTSPOT_OFFSET 4
#define HOTSPOT_ASSOCIATIVITY 4
#define HOTSPOT_AVERAGE_MAX 256 //?
#define GHR_LENGTH 6

// #define PRINT_HOTSPOT
// #define PRINT_TRACES
#define PREDICT_TRACE
#define TRACE_LENGTH 6
#define TRACE_PRED_TRIG_THRES 128
#define TRACE_INST_LENG_THRES 32
#define TRACE_INST_LENG_MAX 1024
#define TRACE_CONF_THRES 0.6
#define TRACE_CONF_EVICT 0.1

#define TRACE_CACHE_SIZE 3 // because we only have three address between two instruction addresses

typedef struct
{
  // uint32_t cond_br[TRACE_LENGTH] = {0};
  // uint64_t nb_early_exit[TRACE_LENGTH] = {0};
  map<int, uint32_t> block_length;
  uint32_t nb_instr_tot = 0;
  // uint64_t trace_id = 0;
  // double precision = 0.0;
  float confidence = 0.0;
  uint16_t length = 0;
  uint32_t pred;
  uint32_t count_use = 0;

} trace_t;

bool getpred(trace_t trace, int i)
{
  if (i >= trace.length)
    return false;
  return (trace.pred >> i) & 1;
}

bool operator==(const trace_t &lhs, const trace_t &rhs)
{
  // if (lhs.trace_id != rhs.trace_id)
  //   return false;
  // if (lhs.confidence != rhs.confidence)
  //   return false;
  if (lhs.length != rhs.length)
    return false;
  if (lhs.pred != rhs.pred)
    return false;
  // if (lhs.count_use != rhs.count_use)
  //   return false;
  return true;
}

/////////////////////////////////////////////////////////////////

// int misprepred[NB_PRE_PRED] = {0};
// int prepreds[NB_PRE_PRED][NB_PRE_PRED+1] = {{0}};
// float trace_mispred_ponder = 0.f;
// uint32_t id_circ_ppred = 0;

void CheckHeartBeat(UINT64 numIter, UINT64 numMispred)
{
  UINT64 dotInterval = 1000000;
  UINT64 lineInterval = 30 * dotInterval;

  UINT64 d1K = 1000;
  UINT64 d10K = 10000;
  UINT64 d100K = 100000;
  UINT64 d1M = 1000000;
  UINT64 d10M = 10000000;
  UINT64 d30M = 30000000;
  UINT64 d60M = 60000000;
  UINT64 d100M = 100000000;
  UINT64 d300M = 300000000;
  UINT64 d600M = 600000000;
  UINT64 d1B = 1000000000;
  UINT64 d10B = 10000000000;

  bool heartbeat = false;
  //  if(numIter % lineInterval == 0){ //prints line every 30 million branches
  //    printf("\n");
  //    fflush(stdout);
  //  }
  if (numIter == d1K)
  { // prints MPKI after 100K branches
    printf("  MPKBr_1K         \t : %10.4f", 1000.0 * (double)(numMispred) / (double)(numIter));
    fflush(stdout);
    heartbeat = true;
  }

  if (numIter == d10K)
  { // prints MPKI after 100K branches
    printf("  MPKBr_10K         \t : %10.4f", 1000.0 * (double)(numMispred) / (double)(numIter));
    fflush(stdout);
    heartbeat = true;
  }

  if (numIter == d100K)
  { // prints MPKI after 100K branches
    printf("  MPKBr_100K         \t : %10.4f", 1000.0 * (double)(numMispred) / (double)(numIter));
    fflush(stdout);
    heartbeat = true;
  }
  if (numIter == d1M)
  {
    printf("  MPKBr_1M         \t : %10.4f", 1000.0 * (double)(numMispred) / (double)(numIter));
    fflush(stdout);
    heartbeat = true;
  }

  if (numIter == d10M)
  { // prints MPKI after 100K branches
    printf("  MPKBr_10M         \t : %10.4f", 1000.0 * (double)(numMispred) / (double)(numIter));
    fflush(stdout);
    heartbeat = true;
  }

  if (numIter == d30M)
  { // prints MPKI after 100K branches
    printf("  MPKBr_30M         \t : %10.4f", 1000.0 * (double)(numMispred) / (double)(numIter));
    fflush(stdout);
    heartbeat = true;
  }

  if (numIter == d60M)
  { // prints MPKI after 100K branches
    printf("  MPKBr_60M         \t : %10.4f", 1000.0 * (double)(numMispred) / (double)(numIter));
    fflush(stdout);
    heartbeat = true;
  }

  if (numIter == d100M)
  { // prints MPKI after 100K branches
    printf("  MPKBr_100M         \t : %10.4f", 1000.0 * (double)(numMispred) / (double)(numIter));
    fflush(stdout);
    heartbeat = true;
  }

  if (numIter == d300M)
  { // prints MPKI after 100K branches
    printf("  MPKBr_300M         \t : %10.4f", 1000.0 * (double)(numMispred) / (double)(numIter));
    fflush(stdout);
    heartbeat = true;
  }

  if (numIter == d600M)
  { // prints MPKI after 100K branches
    printf("  MPKBr_600M         \t : %10.4f", 1000.0 * (double)(numMispred) / (double)(numIter));
    fflush(stdout);
    heartbeat = true;
  }

  if (numIter == d1B)
  { // prints MPKI after 100K branches
    printf("  MPKBr_1B         \t : %10.4f", 1000.0 * (double)(numMispred) / (double)(numIter));
    fflush(stdout);
    heartbeat = true;
  }

  if (numIter == d10B)
  { // prints MPKI after 100K branches
    printf("  MPKBr_10B         \t : %10.4f", 1000.0 * (double)(numMispred) / (double)(numIter));
    fflush(stdout);
    heartbeat = true;
  }

  if (heartbeat)
  {
    time_t timer;
    char buffer[28];
    struct tm *tm_info;

    timer = time(NULL);
    tm_info = localtime(&timer);

    strftime(buffer, 28, " %Y-%m-%d %H:%M:%S", tm_info);
    puts(buffer);
    fflush(stdout);
  }
  // if (heartbeat) {
  //   //print prepreds
  //   printf("\n");
  //   for (int i = 0; i < NB_PRE_PRED; i ++) {
  //     for (int j = 0; j <= NB_PRE_PRED; j ++ )
  //       printf("%d\t", prepreds[i][j]);
  //     printf("\n");
  //   }
  //   printf("\n");
  //   for(int i = 0; i < NB_PRE_PRED; i ++) {
  //     printf("%d\t", misprepred[i]);
  //   }
  //   printf("\n %d \n", id_circ_ppred);

  // }

} // void CheckHeartBeat

uint64_t mask64(unsigned int length)
{
  return ((uint64_t)1 << length) - 1;
}
uint32_t mask32(unsigned int length)
{
  return ((uint32_t)1 << length) - 1;
}

uint64_t traceID(uint64_t PC, uint32_t GHR)
{
  //   |            |  |  | 16 bits
  //   |------------|--|    PC
  //   |------------|  |--| GHR
  return PC;
  // uint64_t mask_GHR = ((GHR & mask64(GHR_LENGTH))>>2)<<4 | (GHR&3);
  // return (PC<<2) ^ mask_GHR;

  // return PC<<GHR_LENGTH | GHR;
}

uint32_t updateGHR(uint32_t GHR, bool taken)
{
  return ((GHR << 1) | taken) & mask32(GHR_LENGTH + TRACE_LENGTH);
}

int main(int argc, char *argv[])
{

  if (argc != 2)
  {
    printf("usage: %s <trace>\n", argv[0]);
    exit(-1);
  }

  ///////////////////////////////////////////////
  // Init variables
  ///////////////////////////////////////////////

  PREDICTOR brpred; // this instantiates the predictor code
#ifdef NB_PRE_PRED
  long unsigned int misprepred[NB_PRE_PRED] = {0};
  long unsigned int not_reached[NB_PRE_PRED] = {0};
  uint32_t prepred_trace_mis_id[NB_PRE_PRED] = {0};
  long int prepreds[NB_PRE_PRED][NB_PRE_PRED] = {{0}};
  long trace_length[NB_PRE_PRED][NB_PRE_PRED] = {{0}};
  long unsigned int entirely_well_pred = 0;
  double trace_mispred_ponder = 0.0;
  uint32_t id_circ_ppred = 0;
#endif
  hotspot_t hotspotness[(uint64_t)1 << HOTSPOT_KEY_SIZE][HOTSPOT_ASSOCIATIVITY];
  hotspot_t hotspot_tag[(uint64_t)1 << HOTSPOT_KEY_SIZE][HOTSPOT_ASSOCIATIVITY];

  unordered_map<uint64_t, trace_t[TRACE_CACHE_SIZE]> trace_pred(0);
  int count_fail_conf = 0;
  int count_fail_size = 0;
  int count_fail_already_here = 0;
  int count_fail = 0;
  uint32_t GHR = 0;

  for (long unsigned int i = 0; i < (uint64_t)1 << HOTSPOT_KEY_SIZE; i++)
  {
    for (int j = 0; j < HOTSPOT_ASSOCIATIVITY; j++)
      hotspotness[(uint32_t)i][j] = 0;
  }
  ///////////////////////////////////////////////
  // read each trace recrod, simulate until done
  ///////////////////////////////////////////////

  std::string trace_path;
  trace_path = argv[1];
  bt9::BT9Reader bt9_reader(trace_path);

  std::string key = "total_instruction_count:";
  std::string value;
  bt9_reader.header.getFieldValueStr(key, value);
  UINT64 total_instruction_counter = std::stoull(value, nullptr, 0);
  UINT64 current_instruction_counter = 0;
  key = "branch_instruction_count:";
  bt9_reader.header.getFieldValueStr(key, value);
  UINT64 branch_instruction_counter = std::stoull(value, nullptr, 0);
  UINT64 numMispred = 0;
  UINT64 numMispredTrace = 0;
  // ver2    UINT64     numMispred_btbMISS =0;
  // ver2    UINT64     numMispred_btbANSF =0;
  // ver2    UINT64     numMispred_btbATSF =0;
  // ver2    UINT64     numMispred_btbDYN =0;

  UINT64 total_instruction_counter_bis = 0;
  UINT64 cond_branch_instruction_counter = 0;
  // ver2     UINT64 btb_ansf_cond_branch_instruction_counter=0;
  // ver2     UINT64 btb_atsf_cond_branch_instruction_counter=0;
  // ver2     UINT64 btb_dyn_cond_branch_instruction_counter=0;
  // ver2     UINT64 btb_miss_cond_branch_instruction_counter=0;
  UINT64 uncond_branch_instruction_counter = 0;

  // ver2    ///////////////////////////////////////////////
  // ver2    // model simple branch marking structure
  // ver2    ///////////////////////////////////////////////
  // ver2    std::map<UINT64, UINT32> myBtb;
  // ver2    map<UINT64, UINT32>::iterator myBtbIterator;
  // ver2
  // ver2    myBtb.clear();

  ///////////////////////////////////////////////
  // read each trace record, simulate until done
  ///////////////////////////////////////////////
#define TAGE_PRED_MAX 7
  PREDICTOR snd_pred;
  // int pred_count = 0;

  OpType opType;
  UINT64 PC;
  bool branchTaken;
  UINT64 branchTarget;
  UINT64 numIter = 0;

  UINT64 global_trace_instr_nb = 0;
  UINT64 global_trace_use = 0;
  uint64_t global_nb_trace = 0;
  uint64_t global_nb_unused_trace = 0;
  double global_trace_precision = 0;
  uint64_t global_trace_evicted = 0;
  uint64_t hotspot_eviction = 0;

  long long unsigned int trace_instruction_counter = 0;
  trace_t *current_trace = nullptr;
  short unsigned int current_trace_it = 0;

  for (auto it = bt9_reader.begin(); it != bt9_reader.end(); ++it)
  {
    CheckHeartBeat(++numIter, numMispred); // Here numIter will be equal to number of branches read
    // char jauge[JAUGE_LENGTH+1];
    // jauge[JAUGE_LENGTH] = '\0';
    // for (unsigned int i = 0; i < JAUGE_LENGTH; i ++)
    //   if (i <= numIter * JAUGE_LENGTH / branch_instruction_counter)
    //     jauge[i] = '#';
    //   else jauge[i] = ' ';
    // printf("\tproceeds \t : (%s) %10lld\t%6ld\r", jauge, numIter, trace_pred.size());

    try
    {
      bt9::BrClass br_class = it->getSrcNode()->brClass();

      //          bool dirDynamic = (it->getSrcNode()->brObservedTakenCnt() > 0) && (it->getSrcNode()->brObservedNotTakenCnt() > 0); //JD2_2_2016
      //          bool dirNeverTkn = (it->getSrcNode()->brObservedTakenCnt() == 0) && (it->getSrcNode()->brObservedNotTakenCnt() > 0); //JD2_2_2016

      // JD2_2_2016 break down branch instructions into all possible types
      opType = OPTYPE_ERROR;

      if ((br_class.type == bt9::BrClass::Type::UNKNOWN) && (it->getSrcNode()->brNodeIndex()))
      {                        // only fault if it isn't the first node in the graph (fake branch)
        opType = OPTYPE_ERROR; // sanity check
      }
      // NOTE unconditional could be part of an IT block that is resolved not-taken
      //           else if (dirNeverTkn && (br_class.conditionality == bt9::BrClass::Conditionality::UNCONDITIONAL)) {
      //             opType = OPTYPE_ERROR; //sanity check
      //           }
      // JD_2_22 There is a bug in the instruction decoder used to generate the traces
      //           else if (dirDynamic && (br_class.conditionality == bt9::BrClass::Conditionality::UNCONDITIONAL)) {
      //             opType = OPTYPE_ERROR; //sanity check
      //           }
      else if (br_class.type == bt9::BrClass::Type::RET)
      {
        if (br_class.conditionality == bt9::BrClass::Conditionality::CONDITIONAL)
          opType = OPTYPE_RET_COND;
        else if (br_class.conditionality == bt9::BrClass::Conditionality::UNCONDITIONAL)
          opType = OPTYPE_RET_UNCOND;
        else
        {
          opType = OPTYPE_ERROR;
        }
      }
      else if (br_class.directness == bt9::BrClass::Directness::INDIRECT)
      {
        if (br_class.type == bt9::BrClass::Type::CALL)
        {
          if (br_class.conditionality == bt9::BrClass::Conditionality::CONDITIONAL)
            opType = OPTYPE_CALL_INDIRECT_COND;
          else if (br_class.conditionality == bt9::BrClass::Conditionality::UNCONDITIONAL)
            opType = OPTYPE_CALL_INDIRECT_UNCOND;
          else
          {
            opType = OPTYPE_ERROR;
          }
        }
        else if (br_class.type == bt9::BrClass::Type::JMP)
        {
          if (br_class.conditionality == bt9::BrClass::Conditionality::CONDITIONAL)
            opType = OPTYPE_JMP_INDIRECT_COND;
          else if (br_class.conditionality == bt9::BrClass::Conditionality::UNCONDITIONAL)
            opType = OPTYPE_JMP_INDIRECT_UNCOND;
          else
          {
            opType = OPTYPE_ERROR;
          }
        }
        else
        {
          opType = OPTYPE_ERROR;
        }
      }
      else if (br_class.directness == bt9::BrClass::Directness::DIRECT)
      {
        if (br_class.type == bt9::BrClass::Type::CALL)
        {
          if (br_class.conditionality == bt9::BrClass::Conditionality::CONDITIONAL)
          {
            opType = OPTYPE_CALL_DIRECT_COND;
          }
          else if (br_class.conditionality == bt9::BrClass::Conditionality::UNCONDITIONAL)
          {
            opType = OPTYPE_CALL_DIRECT_UNCOND;
          }
          else
          {
            opType = OPTYPE_ERROR;
          }
        }
        else if (br_class.type == bt9::BrClass::Type::JMP)
        {
          if (br_class.conditionality == bt9::BrClass::Conditionality::CONDITIONAL)
          {
            opType = OPTYPE_JMP_DIRECT_COND;
          }
          else if (br_class.conditionality == bt9::BrClass::Conditionality::UNCONDITIONAL)
          {
            opType = OPTYPE_JMP_DIRECT_UNCOND;
          }
          else
          {
            opType = OPTYPE_ERROR;
          }
        }
        else
        {
          opType = OPTYPE_ERROR;
        }
      }
      else
      {
        opType = OPTYPE_ERROR;
      }

      PC = it->getSrcNode()->brVirtualAddr();

      branchTaken = it->getEdge()->isTakenPath();
      branchTarget = it->getEdge()->brVirtualTarget();
      total_instruction_counter_bis += it->getEdge()->nonBrInstCnt();
      // printf("PC: %llx type: %x T %d N %d outcome: %d", PC, (UINT32)opType, it->getSrcNode()->brObservedTakenCnt(), it->getSrcNode()->brObservedNotTakenCnt(), branchTaken);

      /************************************************************************************************************/

      if (opType == OPTYPE_ERROR)
      {
        if (it->getSrcNode()->brNodeIndex())
        { // only fault if it isn't the first node in the graph (fake branch)
          fprintf(stderr, "%s: OPTYPE_ERROR\n", trace_path.c_str());
          printf("%s: OPTYPE_ERROR\n", trace_path.c_str());
          exit(-1); // this should never happen, if it does please email CBP org chair.
        }
      }
      else if (br_class.conditionality == bt9::BrClass::Conditionality::CONDITIONAL)
      {

#ifdef TRACE_LENGTH
        // HOTSPOT TRACE PREDICTION STAT ///////////////////////////////////////////////////
        // printf(" %10lld\t%6ld\t%1u\r", numIter, global_nb_trace - global_trace_evicted, GHR);

        if (current_trace != nullptr)
        {
          // if mispred -> quit trace

          bool trace_br = getpred(*current_trace, current_trace_it);
          bool is_taken = it->getEdge()->isTakenPath();
          if (trace_br != is_taken)
          {
            // TODO all sorts of metrics
            // current_trace->nb_early_exit[current_trace_it]++;

            // update trace
            current_trace->confidence *= 1.0 - 1.0 / (double)(1 << (current_trace_it));

            // if (current_trace->confidence < TRACE_CONF_EVICT) {  //evict : not relevant ?
            //   uint64_t key = (current_trace->trace_id >> HOTSPOT_OFFSET) & mask64(HOTSPOT_KEY_SIZE);
            //   trace_pred[key].remove(*current_trace);
            //   global_trace_evicted ++;
            // }

            // leave trace
            current_trace = nullptr;
            current_trace_it = 0;
          }
          else
          {
            // trace_instruction_counter += current_trace->block_length[current_trace_it];
            current_trace->confidence += max(1.0f, current_trace->confidence * (1 + 1 / (1 << current_trace_it)));
            global_trace_precision += ((double)current_trace->block_length[current_trace_it])/(double)current_trace->nb_instr_tot;
            current_trace_it++;
            if (current_trace_it >= current_trace->length)
            {
              current_trace = nullptr;
              current_trace_it = 0;
            }
          }
        }
#endif

        // CONSTRUCT THE GOOD TRACE PREDICTION /////////////////////////////////////////////
        int max_length = 0;
        for (int i = 0; i < TRACE_CACHE_SIZE; i++)
        {
          int length = trace_pred[PC][i].length;
          if (length > max_length)
            max_length = length;
        }
        // read the future //////
        bt9::BT9Reader::BranchInstanceIterator it_br = it;
        list<bool> future_cond;
        int future_length = max_length;

        future_cond.push_back(it_br->getEdge()->isTakenPath());
        it_br++;
        for (int i = 0; i < future_length && it_br != bt9_reader.end(); i++)
        {
          while (it_br != bt9_reader.end() && it_br->getSrcNode()->brClass().conditionality == bt9::BrClass::Conditionality::UNCONDITIONAL)
          {
            if (it_br->getSrcNode()->brClass().directness == bt9::BrClass::Directness::INDIRECT)
            {
              future_length = i + 1;
              break;
            }

            it_br++;
          }
          future_cond.push_back(it_br->getEdge()->isTakenPath());
          it_br++;
        }

        // choose best pred /////
        char best_pred = 0;
        char best_length = 0;

        for (int i = 0; i < TRACE_CACHE_SIZE; i++)
        {

          trace_t trace = trace_pred[PC][i];
          if (trace.confidence < TRACE_CONF_EVICT)
            continue;

          int length = 0;
          auto it_future = future_cond.begin();
          while (length < future_length && getpred(trace, length) == *it_future)
            length++;

          if (length > best_length)
          {
            best_pred = i + 1;
            best_length = length;
          }
        }

        // PREDICTION //////////////////////////////////////////////////////////////////////
        bool predDir = false;

        predDir = brpred.GetPrediction(PC);
        brpred.UpdatePredictor(PC, branchTaken, predDir, branchTarget);

        int branch_to = 0;
        float confidence = 0;

        for (int i = 1; i <= TRACE_CACHE_SIZE; i++)
        {
          if (trace_pred[PC][i - 1].confidence < TRACE_CONF_EVICT)
            continue;

          bool predDir_bis = brpred.GetPrediction(PC + i);
          float conf = brpred.getConfidence();
          brpred.UpdatePredictor(PC, best_pred == i, predDir_bis, branchTarget);
          numMispredTrace += predDir_bis && (best_pred != i);
          if (predDir_bis && confidence < conf)
          {
            branch_to = i; // hope for only one, other wise the last will be taken
            confidence = conf;
          }
        }

        if (current_trace == nullptr && branch_to > 0)
        { // not in a trace and predict to jump in trace (correct to the future)
          current_trace = &(trace_pred[PC][branch_to - 1]);
          current_trace_it = 1;
          // trace_instruction_counter += current_trace->block_length[0];
          global_trace_use++;
        }

        numMispred += predDir && (best_pred != 0);

        // if (predDir != branchTaken)
        // {
        //   numMispred++; // update mispred stats
        // }

        cond_branch_instruction_counter++;

#ifdef TRACE_LENGTH
        // HOTSPOT /////////////////////////////////////////////////////////////////////////

        // | branch target                                                 |
        // |         null         |      tag      |      key      |  null  |

        hotspot_t tag = (branchTarget >> (HOTSPOT_KEY_SIZE + HOTSPOT_OFFSET)) & mask64(HOTSPOT_TAG_SIZE);
        uint32_t key = (branchTarget >> HOTSPOT_OFFSET) & mask64(HOTSPOT_KEY_SIZE);

        // TRACE ID /////
        uint64_t targetID = PC; // traceID(branchTarget, GHR);
        uint32_t trace_key = (targetID >> HOTSPOT_OFFSET) & mask64(HOTSPOT_KEY_SIZE);

        GHR = updateGHR(GHR, branchTaken);

        // find corresponding way in hotspot cache
        int way = -1;
        hotspot_t counter_sum = 0;
        for (int it_way = 0; it_way < HOTSPOT_ASSOCIATIVITY; it_way++)
        {
          if (hotspot_tag[key][it_way] == tag)
          {
            way = it_way;
            // if (hotspotness[key][it_way] + 1 != 0) // not max yet, redondant with lfu
            hotspotness[key][it_way]++;
          }
          counter_sum += hotspotness[key][it_way];
        }
        // LFU
        if (counter_sum / HOTSPOT_ASSOCIATIVITY > HOTSPOT_AVERAGE_MAX)
          for (int it_way = 0; it_way < HOTSPOT_ASSOCIATIVITY; it_way++)
            hotspotness[key][it_way] = (hotspotness[key][it_way] + 1) / 2;
        // ceiled half

        // HOTSPOT CACHE ALLOCATION (and eviction)
        if (way == -1)
        { // is not in cache
          bool free_way = false;
          for (int i = 0; i < HOTSPOT_ASSOCIATIVITY; i++)
            if (hotspotness[key][i] == 0)
            { // place it in this entry
              hotspot_tag[key][i] = tag;
              hotspotness[key][i] = 1;
              free_way = true;
            }

          if (!free_way)
          { // eviction
            uint16_t min_counter = hotspotness[key][0];
            uint16_t evict_way = 0;

            for (int i = 0; i < HOTSPOT_ASSOCIATIVITY; i++)
              if (hotspotness[key][i] < min_counter)
              {
                evict_way = i;
                min_counter = hotspotness[key][i];
              }

            hotspot_eviction++;
            // evict selected way
            hotspotness[key][evict_way] = 1;
            hotspot_tag[key][evict_way] = tag;
          }
        } // end of allocation (and eviction)

        // search for a trace from here
        bool room_for_a_trace = false;
        int trace_evict_id = 0;
        double min_conf = 1.0;
        for (int i = 0; i < TRACE_CACHE_SIZE; i++)
        {
          if (trace_pred[targetID][i].confidence < TRACE_CONF_EVICT)
          {
            room_for_a_trace = true;
            if (min_conf > trace_pred[targetID][i].confidence)
            {
              trace_evict_id = i;
              min_conf = trace_pred[targetID][i].confidence;
            }
          }
        }

        if (hotspotness[key][way] >= TRACE_PRED_TRIG_THRES && current_trace == nullptr && room_for_a_trace)
        // trigger trace construction
        {

          // TRACE CONSTRUCTION ///////////////////////////////////////////////////////////////
          trace_t new_trace;
          new_trace.confidence = 1.0;

          bool useful = true;
          uint64_t trace_length_instr = 0;

          bool prepred_dir = branchTaken;
          bt9::BT9Reader::NodeTableIterator node_it = bt9_reader.node_table.begin();
          node_it += it->getSrcNode()->brNodeIndex();
          // node_it.nextConditionalNode(prepred_dir); // target the next conditional br

          uint64_t pc_pred = node_it->brVirtualAddr();

          snd_pred = brpred; // copy_ctor
          // snd_pred.UpdatePredictor(PC, predDir, branchTaken, 0);
          // take true information as if the prediction is made after branching
          
          int i = 0;
          while (new_trace.confidence > TRACE_CONF_THRES)
          {
            prepred_dir = snd_pred.GetPrediction(pc_pred);
            int next_edge = node_it.nextConditionalNode(prepred_dir);

            if (next_edge < 0 || trace_length_instr + node_it.getPathInstrucCount() > TRACE_INST_LENG_MAX)
            { // program's end or indirect path
              // therefore unpredictible

              // or too many instruction in trace
              break; // stops construction
            }

            new_trace.confidence *= snd_pred.getConfidence();
            // new_trace.cond_br[i] = (unsigned)next_edge;
            new_trace.pred |= prepred_dir << i;
            new_trace.block_length[i] = node_it.getPathInstrucCount();
            trace_length_instr += node_it.getPathInstrucCount();
            new_trace.length = i + 1;

            snd_pred.UpdatePredictor(pc_pred, prepred_dir, prepred_dir, 0);
            pc_pred = node_it->brVirtualAddr();
            i++;
          }

          // if (trace_length_instr >= TRACE_INST_LENG_THRES && new_trace.confidence >= TRACE_CONF_THRES)
          if (new_trace.confidence >= TRACE_CONF_THRES && trace_length_instr > TRACE_INST_LENG_THRES)
          {

            bool already_here = false;
            for (int i = 0; i < TRACE_CACHE_SIZE; i++)
              if (new_trace == trace_pred[targetID][trace_evict_id])
                already_here = true;
            if (already_here)
            {
              count_fail_already_here++;
              count_fail++;
            }
            else
            {
              // insertion
              new_trace.nb_instr_tot = trace_length_instr;
              trace_pred[targetID][trace_evict_id] = new_trace;
              global_trace_instr_nb += trace_length_instr;
              global_nb_trace++;
            }
          }
          else
            count_fail++;

          if (trace_length_instr < TRACE_INST_LENG_THRES)
            count_fail_size++;
          if (new_trace.confidence < TRACE_CONF_THRES)
            count_fail_conf++;
          // END TRACE CONSTRUCTION ///////////////////////////////////////////////////////////
        }
#endif

        // ver2            if (btbDYN)
        // ver2              btb_dyn_cond_branch_instruction_counter++; //number of branches that have been N at least once after being T at least once
        // ver2            else if (btbATSF)
        // ver2              btb_atsf_cond_branch_instruction_counter++; //number of branches that have been T at least once, but have not yet seen a N after the first T
        // ver2            else if (btbANSF)
        // ver2              btb_ansf_cond_branch_instruction_counter++; //number of cond branches that have not yet been observed T
        // ver2            else
        // ver2              btb_miss_cond_branch_instruction_counter++; //number of cond branches that have not yet been observed T
      }
      else if (br_class.conditionality == bt9::BrClass::Conditionality::UNCONDITIONAL)
      { // for predictors that want to track unconditional branches
        uncond_branch_instruction_counter++;
        brpred.TrackOtherInst(PC, opType, branchTarget);
      }
      else
      {
        fprintf(stderr, "CONDITIONALITY ERROR\n");
        printf("CONDITIONALITY ERROR\n");
        exit(-1); // this should never happen, if it does please email CBP org chair.
      }

      if (current_trace != nullptr){
        trace_instruction_counter += it->getEdge()->nonBrInstCnt() + 1;
      }
      
      /************************************************************************************************************/
    }
    catch (const std::out_of_range &ex)
    {
      std::cout << ex.what() << '\n';
      break;
    }
    catch (const std::exception &ex)
    {
      std::cout << "\ndunno exception\n"
                << ex.what() << '\n';
      break;
    }

  } // for (auto it = bt9_reader.begin(); it != bt9_reader.end(); ++it)

  total_instruction_counter_bis += cond_branch_instruction_counter + uncond_branch_instruction_counter;

  ///////////////////////////////////////////
  // print_stats
  ///////////////////////////////////////////
  // NOTE: competitors are judged solely on MISPRED_PER_1K_INST. The additional stats are just for tuning your predictors.
#ifdef PRINT_HOTSPOT
  for (uint64_t i = 0; i < ((unsigned long int)1 << HOTSPOT_KEY_SIZE); i++)
  {
    for (int j = 0; j < HOTSPOT_ASSOCIATIVITY; j++)
      if (hotspotness[(uint32_t)i][j] > 0)
      {
        printf("0x%08x\t", i << HOTSPOT_OFFSET + (hotspot_tag[i][j]) << (HOTSPOT_OFFSET + HOTSPOT_KEY_SIZE));
        printf("%u\n", hotspotness[(uint32_t)i][j]);
      }
  }
#endif
  printf("  TRACE \t : %s\n", trace_path.c_str());
  printf("  NUM_INSTRUCTIONS            \t : %10llu (%10llu)\n", total_instruction_counter, total_instruction_counter_bis);
  printf("  NUM_BR                      \t : %10llu (%10llu)\n", branch_instruction_counter - 1, cond_branch_instruction_counter + uncond_branch_instruction_counter); // JD2_2_2016 NOTE there is a dummy branch at the beginning of the trace...
  printf("  NUM_UNCOND_BR               \t : %10llu\n", uncond_branch_instruction_counter);
  printf("  NUM_CONDITIONAL_BR          \t : %10llu\n", cond_branch_instruction_counter);
  // ver2      printf("  NUM_CONDITIONAL_BR_BTB_MISS \t : %10llu",   btb_miss_cond_branch_instruction_counter);
  // ver2      printf("  NUM_CONDITIONAL_BR_BTB_ANSF \t : %10llu",   btb_ansf_cond_branch_instruction_counter);
  // ver2      printf("  NUM_CONDITIONAL_BR_BTB_ATSF \t : %10llu",   btb_atsf_cond_branch_instruction_counter);
  // ver2      printf("  NUM_CONDITIONAL_BR_BTB_DYN  \t : %10llu",   btb_dyn_cond_branch_instruction_counter);
  printf("  NUM_MISPREDICTIONS          \t : %10llu (%10llu)\n", numMispred, numMispredTrace);

  // ver2      printf("  NUM_MISPREDICTIONS_BTB_MISS \t : %10llu",   numMispred_btbMISS);
  // ver2      printf("  NUM_MISPREDICTIONS_BTB_ANSF \t : %10llu",   numMispred_btbANSF);
  // ver2      printf("  NUM_MISPREDICTIONS_BTB_ATSF \t : %10llu",   numMispred_btbATSF);
  // ver2      printf("  NUM_MISPREDICTIONS_BTB_DYN  \t : %10llu",   numMispred_btbDYN);
  printf("  MISPRED_PER_1K_INST         \t : %10.6f (%10.6f)\n", 1000.0 * (double)(numMispred) / (double)(total_instruction_counter), 1000.0 * (double)(numMispredTrace) / (double)(total_instruction_counter));
#ifdef NB_PRE_PRED
  for (int i = 0; i < NB_PRE_PRED; i++)
  {
    printf("    MISPREPRED_PER_1K_INST %2d\t : %10.6f\t %3.4f% \t %3.4f% \n", i + 1,
           1000.0 * (double)(misprepred[i]) / (double)(total_instruction_counter),
           100.0 - 100.0 * (double)(not_reached[i]) / (double)(cond_branch_instruction_counter),
           100.0 * (double)(misprepred[i]) / (double)(cond_branch_instruction_counter - not_reached[i]));
  }
  printf("  WELL_PRED_TRACE             \t : %10.4f %\n", 100.0 - 100.0 * (double)(trace_mispred_ponder) / (double)(cond_branch_instruction_counter));
  printf("  ENTIRELY_WELL_PRED_TRACE    \t : %10lu\n", entirely_well_pred);
#endif
  printf("\n");

#ifdef TRACE_LENGTH

  unordered_map<uint64_t, trace_t[3]>::iterator traces = trace_pred.begin();
  // long long unsigned int trace_instr_exec = 0;
  uint32_t nb_nonnull_trace = 0;
  uint32_t nb_cache_entry = 0;
  while (traces != trace_pred.end())
  {
    nb_cache_entry++;
    for (int it_trace = 0; it_trace < TRACE_CACHE_SIZE; it_trace++)
    {
      if (traces->second[it_trace].count_use == 0 && traces->second[it_trace].confidence > 0)
        global_nb_unused_trace++;
      nb_nonnull_trace += traces->second[it_trace].confidence > 0;
#ifdef PRINT_TRACES
      printf("    TRACE Ox%016lx \t", it_trace->trace_id);
      // for (int i = 0; i < it_trace->length; i++)
      // {
      //   printf("%5d ", it_trace->pred[i]);
      // }
      // for (int i = 0; i < TRACE_LENGTH - it_trace->length; i++)
      // {
      //   printf("     ");
      // }

      printf(": (CONF:%1.4f) (USE:%6d) ", it_trace->confidence, it_trace->count_use);

      printf("\n");
#endif
    }
    traces++;
  }

  double prec_display = 0.0;
  double instr_display = 0.0;
  double aver_nb_trace_cache = 0.0;

  if (global_nb_trace > 0)
  { // just in case ...
    prec_display = (double)(global_trace_precision) / (double)(global_trace_use);
    instr_display = (double)(global_trace_instr_nb) / (double)(global_nb_trace);
    aver_nb_trace_cache = (double)(nb_nonnull_trace) / (double)(nb_cache_entry);
  }

  printf("  TRACE_NUMBER                \t : %10ld\n", (global_nb_trace));
  printf("  CONSTRUCTION_FAILS          \t : %10d\n", count_fail);
  printf("  CONSTRUCTION_FAILS_CONF     \t : %10d\n", count_fail_conf);
  printf("  CONSTRUCTION_FAILS_SIZE     \t : %10d\n", count_fail_size);
  printf("  CONSTRUCTION_FAILS_HERE     \t : %10d\n", count_fail_already_here);
  printf("  HOTSPOT_EVICTION            \t : %10ld\n", hotspot_eviction);
  if (!trace_pred.empty())
  {
    printf("  MEAN_TRACE_INSTRUCTION      \t : %10.2f\n", instr_display);
    printf("  MEAN_TRACE_PRECISION        \t : %10.6f\n", prec_display);
    printf("  MEAN_TRACE_PER_CACHE        \t : %10.6f\n", aver_nb_trace_cache);
    printf("  MEAN_TRACE_USAGE            \t : %10.6f\n", (double)trace_instruction_counter / (double)total_instruction_counter_bis);
    printf("  UNUSED_TRACE_PER            \t : %10.6f\n", 100 * (double)global_nb_unused_trace / (double)(global_nb_trace));
  }
#endif

  // ver2      printf("  MISPRED_PER_1K_INST_BTB_MISS\t : %10.4f",   1000.0*(double)(numMispred_btbMISS)/(double)(total_instruction_counter));
  // ver2      printf("  MISPRED_PER_1K_INST_BTB_ANSF\t : %10.4f",   1000.0*(double)(numMispred_btbANSF)/(double)(total_instruction_counter));
  // ver2      printf("  MISPRED_PER_1K_INST_BTB_ATSF\t : %10.4f",   1000.0*(double)(numMispred_btbATSF)/(double)(total_instruction_counter));
  // ver2      printf("  MISPRED_PER_1K_INST_BTB_DYN \t : %10.4f",   1000.0*(double)(numMispred_btbDYN)/(double)(total_instruction_counter));
  printf("\n");
}
