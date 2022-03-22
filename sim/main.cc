///////////////////////////////////////////////////////////////////////
//  Copyright 2015 Samsung Austin Semiconductor, LLC.                //
///////////////////////////////////////////////////////////////////////

// Description : Main file for CBP2016

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <deque>
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
#define HOTSPOT_TAG_SIZE 20
typedef uint16_t hotspot_t ;
#define HOTSPOT_OFFSET 8

// #define PRINT_HOTSPOT
#define PREDICT_TRACE
#define TRACE_PRED_THRES 4000

typedef struct trace_descr{
  uint32_t cond_br[TRACE_LENGTH] = {0};
  uint32_t block_length[TRACE_LENGTH] = {0};
  uint64_t nb_early_exit[TRACE_LENGTH] = {0}; // [0] -> nb prediction
  uint64_t trace_id = 0;
  double precision = 0.0;
  uint32_t nb_instr_tot = 0;
  bool pred[TRACE_LENGTH] = {0};
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

uint64_t mask64(unsigned int length) {
  return ((uint64_t)1 << length)-1;
}
uint32_t mask32(unsigned int length) {
  return ((uint32_t)1 << length)-1;
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
  hotspot_t hotspotness[(uint64_t)1 << HOTSPOT_TAG_SIZE];
  uint64_t trace_trigger_address = 0;
  vector<trace_t> trace_pred(0);
  std::deque<uint32_t> cond_branch_histo(TRACE_LENGTH);
  int count_fail = 0;
  
  for (long unsigned int i = 0; i < (uint64_t)1 << HOTSPOT_TAG_SIZE; i ++) {
    hotspotness[(uint32_t)i] = 0;
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

  PREDICTOR snd_pred;
  int pred_count = 0;

  OpType opType;
  UINT64 PC;
  bool branchTaken;
  UINT64 branchTarget;
  UINT64 numIter = 0;

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

        // HOTSPOT TRACE PREDICTION STAT ///////////////////////////////////////////////////


        if (cond_branch_histo.size() >= TRACE_LENGTH) {
          cond_branch_histo.pop_front();
        }
        cond_branch_histo.push_back(it->getEdge()->edgeIndex());

        vector<trace_t>::iterator it_trace = trace_pred.begin();
        while (it_trace != trace_pred.end()) {
          
          if (it_trace->cond_br[0] == cond_branch_histo[0]) { // only compare on the whole histo, no partial overlap
            uint32_t num = it_trace->block_length[0];

            it_trace->nb_early_exit[0]++;
            
            // for (int id_br = 0; id_br < cond_branch_histo.size(); id_br ++) {
            for (int id_br = 1; id_br < TRACE_LENGTH; id_br ++) {
                            
              if (it_trace->cond_br[id_br] != cond_branch_histo[id_br]) {
                it_trace->nb_early_exit[id_br] ++;
                break;
              } else {
                num += it_trace->block_length[id_br];
              }
            }
            if (it_trace->nb_instr_tot > 0)
              it_trace->precision+= (double)num/(double)it_trace->nb_instr_tot;

          } // end of matching trace


          it_trace ++;
        }
        
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

        // HOTSPOT /////////////////////////////////////////////////////////////////////////

        // | branch target                                                 |
        // |            null               |         tag          |  null  |

        // hotspot_t key = (branchTarget>>(HOTSPOT_TAG_SIZE + HOTSPOT_OFFSET)) & mask64(HOTSPOT_KEY_SIZE);
        uint32_t  tag = (branchTarget>>HOTSPOT_OFFSET) & mask64(HOTSPOT_TAG_SIZE);//(((uint64_t)1 << HOTSPOT_TAG_SIZE) -1);
        
        // if (hotspot_key[tag] == key) {
        if (hotspotness[tag] + 1 != 0)
          hotspotness[tag] ++;

        // trigger trace prediction
        if (hotspotness[tag] >= TRACE_PRED_THRES) {

          // search for a prediction from here
          bool no_pred_from_here = true;
          std::vector<trace_t>::iterator it_test = trace_pred.begin();
          while (it_test != trace_pred.end()) {
            if (it_test->trace_id == PC + branchTaken) { 
              no_pred_from_here = false;
              break;
            }
            it_test ++;
          }

          if (no_pred_from_here) { // no prediction here yet

            trace_t new_trace;
            new_trace.trace_id = PC + branchTaken;
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

            // new_trace.cond_br[0] = it->getEdge()->edgeIndex();
            // new_trace.block_length[0] = node_it.getPathInstrucCount();
            // new_trace.nb_instr_tot += new_trace.block_length[0];
            // for (int i = 1; i < TRACE_LENGTH; i++)

            for (int i = 0; i < TRACE_LENGTH; i++)
            {
              prepred_dir = snd_pred.GetPrediction(pc_pred);
              int next_edge = node_it.nextConditionalNode(prepred_dir);

              if (next_edge < 0) { // program's end or indirect path
                // therefore unpredictible
                useful = false;
                count_fail++;
                hotspotness[tag] /= 2;
                break; // stops construction
              }

              new_trace.cond_br[i] = (unsigned)next_edge;
              new_trace.block_length[i] = node_it.getPathInstrucCount();
              new_trace.nb_instr_tot += new_trace.block_length[i];

              if (i < TRACE_LENGTH - 1)
              {
                snd_pred.UpdatePredictor(pc_pred, prepred_dir, prepred_dir, 0);
                pc_pred = node_it->brVirtualAddr();
              }
            }

            if (useful) 
              trace_pred.push_back(new_trace);
          }
        }
        // } else {
          
        //   uint64_t i;
        //   for (i = 0; i < ((unsigned long int)1<<HOTSPOT_TAG_SIZE); i++) {
        //     hotspotness[i]/=2;
        //   }
        //   // attribute entry
        //   if (hotspotness[tag] == 0) {
        //     hotspotness[tag] ++;
        //     hotspot_key[tag] = key;
        //   }
        // }
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
  for (uint64_t i = 0; i < ((unsigned long int)1<<HOTSPOT_TAG_SIZE); i ++) {
    if (hotspotness[(uint32_t)i] > 0){
      // printf("0x%016lx\t", ((unsigned long)hotspot_key[(uint32_t)i]<< (HOTSPOT_TAG_SIZE+HOTSPOT_OFFSET)) | i<<HOTSPOT_OFFSET);
      printf("0x%08x\t",  i<<HOTSPOT_OFFSET);
      printf("%u\n", hotspotness[(uint32_t)i]);
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

  uint32_t instr_sum = 0;
  double   prec_sum  = 0;
  uint64_t nb_early_exit[TRACE_LENGTH] = {0}; // [0] -> nb prediction

  std::vector<trace_t>::iterator it_trace = trace_pred.begin();
  while (it_trace != trace_pred.end()) {
    instr_sum += it_trace->nb_instr_tot;
    
    double prec = 0;
    if (it_trace->nb_early_exit[0] > 0)
      prec = (it_trace->precision) / (double)(it_trace->nb_early_exit[0]);
    prec_sum += prec;

    printf("    TRACE %16lx ", it_trace->trace_id);
    for (int i = 0; i < TRACE_LENGTH; i ++) 
      printf("%4d ", it_trace->cond_br[i]);
    printf("\t: (PRECISION: %1.6f) ",  prec);
    
    for (int i = 0; i < TRACE_LENGTH; i ++) {
      nb_early_exit[i]+= it_trace->nb_early_exit[i];
      printf("%8ld ", it_trace->nb_early_exit[i]);
    }
    printf("\n");
    it_trace++;
  }

  if (!trace_pred.empty()) {
    printf("  MEAN_TRACE_INSTRUCTION         \t : %10.6f\n",  (double)(instr_sum) / (double)(trace_pred.size()));
    printf("  MEAN_TRACE_PRECISION           \t : %10.6f\n",  (double)(prec_sum) / (double)(trace_pred.size()));
  }
  // ver2      printf("  MISPRED_PER_1K_INST_BTB_MISS\t : %10.4f",   1000.0*(double)(numMispred_btbMISS)/(double)(total_instruction_counter));
  // ver2      printf("  MISPRED_PER_1K_INST_BTB_ANSF\t : %10.4f",   1000.0*(double)(numMispred_btbANSF)/(double)(total_instruction_counter));
  // ver2      printf("  MISPRED_PER_1K_INST_BTB_ATSF\t : %10.4f",   1000.0*(double)(numMispred_btbATSF)/(double)(total_instruction_counter));
  // ver2      printf("  MISPRED_PER_1K_INST_BTB_DYN \t : %10.4f",   1000.0*(double)(numMispred_btbDYN)/(double)(total_instruction_counter));
  printf("\n");
}
