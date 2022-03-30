///////////////////////////////////////////////////////////////////////
//  Copyright 2015 Samsung Austin Semiconductor, LLC.                //
///////////////////////////////////////////////////////////////////////

// Description : Main file for CBP2016

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <list>
using namespace std;

#include "utils.h"
//#include "bt9.h"
#include "bt9_reader.h"
//#include "predictor.cc"
#include "predictor.h"

#define COUNTER unsigned long long
//#define NB_PRE_PRED 6

// HOTSPOT STUDY ////////////////////////////////////////////////
#define TRACE_LENGTH 6
#define HOTSPOT_KEY_SIZE 10
#define HOTSPOT_TAG_SIZE 10
typedef uint16_t hotspot_t;
#define HOTSPOT_OFFSET 4
#define HOTSPOT_ASSOCIATIVITY 2
#define HOTSPOT_AVERAGE_MAX 20000 //?

// #define PRINT_HOTSPOT
#define PRINT_TRACES
#define PREDICT_TRACE
#define TRACE_PRED_TRIG_THRES 5000
#define TRACE_INST_LENG_THRES 256
#define TRACE_CONF_THRES 0.7

typedef struct trace_descr
{
  uint32_t cond_br[TRACE_LENGTH] = {0};
  uint32_t block_length[TRACE_LENGTH] = {0};
  uint64_t nb_early_exit[TRACE_LENGTH] = {0}; // [0] -> nb prediction (not sure)
  uint64_t trace_id = 0;
  double precision = 0.0;
  float confidence = 1.0;
  uint32_t nb_instr_tot = 0;
  uint16_t length = 0;
  bool pred[TRACE_LENGTH] = {0};
  uint32_t count_use = 0;
} trace_t;

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
// usage: predictor <trace>

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

  list<trace_t> trace_pred(0);
  int count_fail = 0;

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

  // ver2    UINT64     numMispred_btbMISS =0;
  // ver2    UINT64     numMispred_btbANSF =0;
  // ver2    UINT64     numMispred_btbATSF =0;
  // ver2    UINT64     numMispred_btbDYN =0;

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
  int pred_count = 0;

  OpType opType;
  UINT64 PC;
  bool branchTaken;
  UINT64 branchTarget;
  UINT64 numIter = 0;
  trace_t *current_trace = nullptr;
  short unsigned int current_trace_it = 0;

  for (auto it = bt9_reader.begin(); it != bt9_reader.end(); ++it)
  {
    CheckHeartBeat(++numIter, numMispred); // Here numIter will be equal to number of branches read

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

      // printf("PC: %llx type: %x T %d N %d outcome: %d", PC, (UINT32)opType, it->getSrcNode()->brObservedTakenCnt(), it->getSrcNode()->brObservedNotTakenCnt(), branchTaken);

      /************************************************************************************************************/

      if (opType == OPTYPE_ERROR)
      {
        if (it->getSrcNode()->brNodeIndex())
        { // only fault if it isn't the first node in the graph (fake branch)
          fprintf(stderr, "OPTYPE_ERROR\n");
          printf("OPTYPE_ERROR\n");
          exit(-1); // this should never happen, if it does please email CBP org chair.
        }
      }
      else if (br_class.conditionality == bt9::BrClass::Conditionality::CONDITIONAL)
      { // JD2_17_2016 call UpdatePredictor() for all branches that decode as conditional
        // printf("COND ");

        // NOTE: contestants are NOT allowed to use the btb* information from ver2 of the infrastructure below:
        // ver2             myBtbIterator = myBtb.find(PC); //check BTB for a hit
        // ver2            bool btbATSF = false;
        // ver2            bool btbANSF = false;
        // ver2            bool btbDYN = false;
        // ver2
        // ver2            if (myBtbIterator == myBtb.end()) { //miss -> we have no history for the branch in the marking structure
        // ver2              //printf("BTB miss ");
        // ver2              myBtb.insert(pair<UINT64, UINT32>(PC, (UINT32)branchTaken)); //on a miss insert with outcome (N->btbANSF, T->btbATSF)
        // ver2              predDir = brpred->GetPrediction(PC, btbANSF, btbATSF, btbDYN);
        // ver2              brpred->UpdatePredictor(PC, opType, branchTaken, predDir, branchTarget, btbANSF, btbATSF, btbDYN);
        // ver2            }
        // ver2            else {
        // ver2              btbANSF = (myBtbIterator->second == 0);
        // ver2              btbATSF = (myBtbIterator->second == 1);
        // ver2              btbDYN = (myBtbIterator->second == 2);
        // ver2              //printf("BTB hit ANSF: %d ATSF: %d DYN: %d ", btbANSF, btbATSF, btbDYN);
        // ver2
        // ver2              predDir = brpred->GetPrediction(PC, btbANSF, btbATSF, btbDYN);
        // ver2              brpred->UpdatePredictor(PC, opType, branchTaken, predDir, branchTarget, btbANSF, btbATSF, btbDYN);
        // ver2
        // ver2              if (  (btbANSF && branchTaken)   // only exhibited N until now and we just got a T -> upgrade to dynamic conditional
        // ver2                 || (btbATSF && !branchTaken)  // only exhibited T until now and we just got a N -> upgrade to dynamic conditional
        // ver2                 ) {
        // ver2                myBtbIterator->second = 2; //2-> dynamic conditional (has exhibited both taken and not-taken in the past)
        // ver2              }
        // ver2            }
        // ver2            //puts("");

        // if (cond_branch_instruction_counter > 0) {
        //   if (id_srcNode != it->getSrcNode()->brNodeIndex()) {
        //     wrongsrcnode++;
        //   }
        // }
#ifdef NB_PRE_PRED
        // STATS COLLECT ///////////////////////////////////////////////////////////////////
        for (int i = 0; i < NB_PRE_PRED && i < cond_branch_instruction_counter; i++)
        {

          uint32_t id_j = (id_circ_ppred - i - 1 + NB_PRE_PRED) % NB_PRE_PRED;

          if (prepred_trace_mis_id[i] == NB_PRE_PRED)
          { // no mispred yet
            long int pred = prepreds[i][id_j];

            if (pred != it->getEdge()->edgeIndex())
            {
              misprepred[id_j]++;
              prepred_trace_mis_id[i] = id_j;
            }
          }
          else
          {
            misprepred[id_j]++;
          }
        }

        uint32_t num = 0, den = 0;
        if (prepred_trace_mis_id[id_circ_ppred] < NB_PRE_PRED)
        { // there is a misprediction
          for (int i = 0; i < NB_PRE_PRED; i++)
          {

            uint32_t edgesize = trace_length[id_circ_ppred][i];

            den += edgesize;
            if (i >= prepred_trace_mis_id[id_circ_ppred])
              num += edgesize;

            trace_length[id_circ_ppred][i] = 0;
            prepreds[id_circ_ppred][i] = (long)0; // reset
          }
          if (den > 0)
            trace_mispred_ponder += ((double)num / (double)den);
        }
        else
          entirely_well_pred++;

        prepred_trace_mis_id[id_circ_ppred] = NB_PRE_PRED; // reset
#endif
#ifdef TRACE_LENGTH
        // HOTSPOT TRACE PREDICTION STAT ///////////////////////////////////////////////////

        if (current_trace != nullptr)
        {
          // if mispred -> quit trace
          
          uint32_t trace_br = current_trace->cond_br[current_trace_it];
          uint32_t it_br = it->getEdge()->edgeIndex();
          if (trace_br != it_br)
          {
            // TODO all sorts of metrics
            current_trace->nb_early_exit[current_trace_it]++;
            // leave trace
            current_trace = nullptr;
            current_trace_it = 0;
          }
          else
          {
            current_trace_it++;
            if (current_trace_it >= current_trace->length)
            {
              current_trace->precision += 1;
              current_trace = nullptr;
              current_trace_it = 0;
            }
          }
          //
          // count
        }
#endif

        // PREDICTION //////////////////////////////////////////////////////////////////////
        bool predDir = false;

        predDir = brpred.GetPrediction(PC);

#ifdef NB_PRE_PRED
        // prepredictions
        //  init
        bool prepred_dir = predDir;
        bt9::BT9Reader::NodeTableIterator node_it = bt9_reader.node_table.begin();
        node_it += it->getSrcNode()->brNodeIndex();
        node_it.nextConditionalNode(prepred_dir); // target the next conditional br

        uint64_t pc_pred = node_it->brVirtualAddr();

        snd_pred = brpred;
        snd_pred.UpdatePredictor(PC, predDir, predDir, 0); // assume correct prediction

        for (int i = 0; i < NB_PRE_PRED; i++)
        {
          prepred_dir = snd_pred.GetPrediction(pc_pred);
          prepreds[id_circ_ppred][i] = (long)node_it.nextConditionalNode(prepred_dir);
          trace_length[id_circ_ppred][i] = node_it.getPathInstrucCount();

          if (prepreds[id_circ_ppred][i] == -2)
          { // program's end or wrong path
            not_reached[i]++;
            while (++i < NB_PRE_PRED)
            {
              prepreds[id_circ_ppred][i] = -2;
              not_reached[i]++;
            }
            break; // stops prepredictions
          }

          if (i < NB_PRE_PRED - 1)
          {
            snd_pred.UpdatePredictor(pc_pred, prepred_dir, prepred_dir, 0);
            pc_pred = node_it->brVirtualAddr();
          }
        }

        id_circ_ppred = (id_circ_ppred + 1) % NB_PRE_PRED;
        // prepredictions end
#endif

        brpred.UpdatePredictor(PC, branchTaken, predDir, branchTarget);

        if (predDir != branchTaken)
        {

          numMispred++; // update mispred stats
                        // ver2              if(btbATSF)
          // ver2                numMispred_btbATSF++; // update mispred stats
          // ver2              else if(btbANSF)
          // ver2                numMispred_btbANSF++; // update mispred stats
          // ver2              else if(btbDYN)
          // ver2                numMispred_btbDYN++; // update mispred stats
          // ver2              else
          // ver2                numMispred_btbMISS++; // update mispred stats
        }
        cond_branch_instruction_counter++;

#ifdef TRACE_LENGTH
        // HOTSPOT /////////////////////////////////////////////////////////////////////////

        // | branch target                                                 |
        // |         null         |      tag      |      key      |  null  |

        hotspot_t tag = (branchTarget >> (HOTSPOT_KEY_SIZE + HOTSPOT_OFFSET)) & mask64(HOTSPOT_TAG_SIZE);
        uint32_t key = (branchTarget >> HOTSPOT_OFFSET) & mask64(HOTSPOT_KEY_SIZE);

        // find corresponding way in hotspot cache
        int way = -1;
        for (int it_way = 0; it_way < HOTSPOT_ASSOCIATIVITY; it_way++)
        {
          hotspot_t counter_sum = 0;
          if (hotspot_tag[key][it_way] == tag)
          {
            way = it_way;
            if (hotspotness[key][it_way] + 1 != 0) // not max yet
              hotspotness[key][it_way]++;

            counter_sum += hotspotness[key][it_way];
            break;
          }
          if (counter_sum / HOTSPOT_ASSOCIATIVITY > HOTSPOT_AVERAGE_MAX)
            for (int it_way = 0; it_way < HOTSPOT_ASSOCIATIVITY; it_way++)
              hotspotness[key][it_way] = (hotspotness[key][it_way] + 1) / 2;
          // ceiled half
        }

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

            // TODO evict trace if any ?
            // trace evictions screw the metrics
            // evict useless traces
            hotspot_t evict_tag = hotspot_tag[key][evict_way];
            list<trace_t>::iterator evict_it = trace_pred.begin();
            while (evict_it != trace_pred.end())
            {
              uint64_t evict_id = (evict_it->trace_id >> HOTSPOT_OFFSET) & mask64(HOTSPOT_KEY_SIZE + HOTSPOT_TAG_SIZE);

              if (evict_id == (key | evict_tag << HOTSPOT_KEY_SIZE) &&
                  evict_it->count_use == 0) // 0 implies it is not the current trace
              {
                evict_it = trace_pred.erase(evict_it);
              }
              else
              {
                evict_it++;
              }
            }

            // evict selected way
            hotspotness[key][evict_way] = 1;
            hotspot_tag[key][evict_way] = tag;
          }
        }

        for (int way = 0; way < HOTSPOT_ASSOCIATIVITY; way++)
        {
          // trigger trace prediction
          if (hotspotness[key][way] >= TRACE_PRED_TRIG_THRES)
          {

            // search for a prediction from here
            bool no_pred_from_here = true;
            list<trace_t>::iterator it_test = trace_pred.begin();
            while (it_test != trace_pred.end())
            {
              if (it_test->trace_id == branchTarget + branchTaken)
              {

                if (current_trace == nullptr)
                {
                  // branch into the trace
                  current_trace = &*it_test;
                  current_trace_it = 0;
                  if (current_trace->count_use +1 != 0) //saturation
                    current_trace->count_use++;
                }

                no_pred_from_here = false;
                break;
              }
              it_test++;
            }

            if (no_pred_from_here)
            { // no prediction here yet

              trace_t new_trace;
              new_trace.trace_id = branchTarget + branchTaken;
              // unique if the br instrs measure more than a byte

              bool useful = true;

              // trace construction
              bool prepred_dir = branchTaken;
              bt9::BT9Reader::NodeTableIterator node_it = bt9_reader.node_table.begin();
              node_it += it->getSrcNode()->brNodeIndex();
              node_it.nextConditionalNode(prepred_dir); // target the next conditional br

              uint64_t pc_pred = node_it->brVirtualAddr();

              snd_pred = brpred;
              snd_pred.UpdatePredictor(PC, predDir, branchTaken, 0);
              // take true information as if the parediction was made after branching

              for (int i = 0; i < TRACE_LENGTH; i++)
              {
                prepred_dir = snd_pred.GetPrediction(pc_pred);
                int next_edge = node_it.nextConditionalNode(prepred_dir);

                if (next_edge < 0)
                { // program's end or indirect path
                  // therefore unpredictible
                  if (new_trace.nb_instr_tot < TRACE_INST_LENG_THRES)
                  {
                    useful = false;
                    count_fail++;
                    hotspotness[key][way] /= 2;
                  } // if trace is long enough, it is kept

                  // cannot go any further
                  break; // stops construction
                }

                new_trace.confidence *= snd_pred.getConfidence();
                new_trace.pred[i] = prepred_dir;
                new_trace.cond_br[i] = (unsigned)next_edge;
                new_trace.block_length[i] = node_it.getPathInstrucCount();
                new_trace.nb_instr_tot += new_trace.block_length[i];
                new_trace.length = i + 1;

                if (i < TRACE_LENGTH - 1)
                {
                  snd_pred.UpdatePredictor(pc_pred, prepred_dir, prepred_dir, 0);
                  pc_pred = node_it->brVirtualAddr();
                }
              }

              if (useful && new_trace.confidence >= TRACE_CONF_THRES)
              {

                trace_pred.push_back(new_trace);

                // branch immediatly in it (if able) or consider a schedule time
              }
            }
          }
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
  printf("%d\n", count_fail);
  printf("  TRACE \t : %s", trace_path.c_str());
  printf("  NUM_INSTRUCTIONS            \t : %10llu\n", total_instruction_counter);
  printf("  NUM_BR                      \t : %10llu\n", branch_instruction_counter - 1); // JD2_2_2016 NOTE there is a dummy branch at the beginning of the trace...
  printf("  NUM_UNCOND_BR               \t : %10llu\n", uncond_branch_instruction_counter);
  printf("  NUM_CONDITIONAL_BR          \t : %10llu\n", cond_branch_instruction_counter);
  // ver2      printf("  NUM_CONDITIONAL_BR_BTB_MISS \t : %10llu",   btb_miss_cond_branch_instruction_counter);
  // ver2      printf("  NUM_CONDITIONAL_BR_BTB_ANSF \t : %10llu",   btb_ansf_cond_branch_instruction_counter);
  // ver2      printf("  NUM_CONDITIONAL_BR_BTB_ATSF \t : %10llu",   btb_atsf_cond_branch_instruction_counter);
  // ver2      printf("  NUM_CONDITIONAL_BR_BTB_DYN  \t : %10llu",   btb_dyn_cond_branch_instruction_counter);
  printf("  NUM_MISPREDICTIONS          \t : %10llu\n", numMispred);

#ifdef NB_PRE_PRED
  for (int i = 0; i < NB_PRE_PRED; i++)
  {
    printf("    NUM_MISPREPREDICTIONS %2d \t : %10lu\t%10lu\n", i + 1, misprepred[i], not_reached[i]);
  }
#endif

  // ver2      printf("  NUM_MISPREDICTIONS_BTB_MISS \t : %10llu",   numMispred_btbMISS);
  // ver2      printf("  NUM_MISPREDICTIONS_BTB_ANSF \t : %10llu",   numMispred_btbANSF);
  // ver2      printf("  NUM_MISPREDICTIONS_BTB_ATSF \t : %10llu",   numMispred_btbATSF);
  // ver2      printf("  NUM_MISPREDICTIONS_BTB_DYN  \t : %10llu",   numMispred_btbDYN);
  printf("  MISPRED_PER_1K_INST         \t : %10.6f\n", 1000.0 * (double)(numMispred) / (double)(total_instruction_counter));
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
  uint32_t instr_sum = 0;
  double prec_sum = 0;
  uint64_t nb_early_exit[TRACE_LENGTH] = {0}; // [0] -> nb prediction

  list<trace_t>::iterator it_trace = trace_pred.begin();
  long long unsigned int trace_instr_exec = 0;
  unsigned int unused_trace = 0;
  while (it_trace != trace_pred.end())
  {
    instr_sum += it_trace->nb_instr_tot;

    double prec = 0;
    if (it_trace->count_use > 0)
      prec = (it_trace->precision) / (double)(it_trace->count_use);
    prec_sum += prec;

    // printf("    TRACE %16lx ", it_trace->trace_id);
    int factor = it_trace->count_use;
    unused_trace += (it_trace->count_use == 0);

    for (int i = 0; i < it_trace->length; i++)
    {
      factor -= it_trace->nb_early_exit[i];
      trace_instr_exec += factor * it_trace->block_length[i];
      nb_early_exit[i] += it_trace->nb_early_exit[i];
    }

#ifdef PRINT_TRACES
    printf("    TRACE ");
    for (int i = 0; i < it_trace->length; i++)
    {
      printf("%5d ", it_trace->cond_br[i]);
    }
    for (int i = 0; i < TRACE_LENGTH - it_trace->length; i++)
    {
      printf("     ");
    }

    printf(": (PREC:%1.4f) (CONF:%1.4f) (USE:%6d) ", prec, it_trace->confidence, it_trace->count_use);

    for (int i = 0; i < it_trace->length; i++)
    {

      printf("%8ld ", it_trace->nb_early_exit[i]);
    }
    printf("\n");
#endif

    it_trace++;
  }
  printf("NB_EARLY_EXIT :");
  for (int i = 0; i < TRACE_LENGTH; i++)
  {
    printf(" %8ld", nb_early_exit[i]);
  }
  printf("\n");

  if (!trace_pred.empty())
  {
    printf("  MEAN_TRACE_INSTRUCTION         \t : %10.6f\n", (double)(instr_sum) / (double)(trace_pred.size()));
    printf("  MEAN_TRACE_PRECISION           \t : %10.6f\n", (double)(prec_sum) / (double)(trace_pred.size()));
    printf("  MEAN_TRACE_USAGE               \t : %10.6f\n", (double)trace_instr_exec / (double)total_instruction_counter);
    printf("  UNUSED_TRACE_PER               \t : %10.6f\n", 100 * (double)unused_trace / (double)(trace_pred.size()));
  }
#endif

  // ver2      printf("  MISPRED_PER_1K_INST_BTB_MISS\t : %10.4f",   1000.0*(double)(numMispred_btbMISS)/(double)(total_instruction_counter));
  // ver2      printf("  MISPRED_PER_1K_INST_BTB_ANSF\t : %10.4f",   1000.0*(double)(numMispred_btbANSF)/(double)(total_instruction_counter));
  // ver2      printf("  MISPRED_PER_1K_INST_BTB_ATSF\t : %10.4f",   1000.0*(double)(numMispred_btbATSF)/(double)(total_instruction_counter));
  // ver2      printf("  MISPRED_PER_1K_INST_BTB_DYN \t : %10.4f",   1000.0*(double)(numMispred_btbDYN)/(double)(total_instruction_counter));
  printf("\n");
}
