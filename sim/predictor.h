// Jayson Boubin, Dec 2017
#ifndef _PREDICTOR_H_
#define _PREDICTOR_H_

#include "utils.h"
// #include "tracer.h"
#include <bitset>
#include <vector>

#define NUM_TAGE_TABLES 12

#define UINT16 unsigned short int

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
// Circular Shift Register for folding purposes
typedef struct csr
{
	UINT16 val;
	UINT16 origLen;
	UINT16 newLen;
} csr_t;

typedef struct bimodVal
{
	UINT32 pred; // prediction (2 bits)
							 // bool m; //metapredictor (1 bit) (eliminated)
} bimodVal_t;

typedef struct tagVal
{
	UINT16 pred;
	UINT16 tag;
	UINT16 u;
} tagVal_t;

typedef struct prediction
{
	bool pred;
	bool altPred;
	int table;
	int altTable;
	UINT32 index;
	UINT32 altIndex;
} prediction_t;

typedef struct loopVal
{
	UINT32 loopCount;		// loop count?
	UINT32 currentIter; // current iteration of the loop
	UINT32 tag;					// n bit tag
	UINT32 conf;				// 2 bit confidence counter
	UINT32 age;
	bool pred;
	bool used;
} loopVal_t;

class PREDICTOR
{

	// The state is defined for Gshare, change for your design

private:
	bitset<1001> GHR; // global history register
	UINT32 PHR;				// path history

	// tables
	vector<bimodVal_t> bimodal; // bimodal table
	// UINT32 numBimodalEntries; // number of entries in bimodal table
	vector<vector<tagVal_t>> tagTables; // TAGE table
	// UINT32 tageTableSize;	              //number of entries in TAGE table
	vector<loopVal_t> loopTable; // loop table
	// UINT32 loopTableSize; // number of loop table entries

	vector<UINT32> tageTableSize;
	vector<UINT32> tageTagSize;

	vector<UINT32> tageHistory;		// number ofGHR bits examined by CSR to index a given table
	vector<csr_t> csrIndex;				// circular shift register for indices
	vector<vector<csr_t>> csrTag; // 2 circular shift registers for tags

	prediction_t pred; // global prediction

	vector<UINT32> tageIndex; // index calculated for a given table
	vector<UINT32> tageTag;		// tag calculated for a given table
	UINT32 clock;							// global clock
	bool clockState;					// clocl flip it
	INT32 altBetterCount;			// number of times altpred is better than prd
	time_t rng_seed;

public:
	// The interface to the four functions below CAN NOT be changed

	PREDICTOR(void);
	PREDICTOR(const PREDICTOR &src) = default;
	~PREDICTOR() = default;
	bool GetPrediction(UINT64 PC);

	void UpdatePredictor(UINT64 PC, bool resolveDir, bool predDir, UINT64 branchTarget);
	void TrackOtherInst(UINT64 PC, OpType opType, UINT64 branchTarget);

	// void steal(UINT32 PC, UINT32 table, UINT32 index, UINT32 bimodalIndex, bool predDir);

	UINT32 getTag(UINT32 PC, int table, UINT64 tagSize);
	UINT32 getIndex(UINT32 PC, int table, UINT64 tagSize, UINT64 phrOffset);
	void initFold(csr_t *shift, UINT64 origLen, UINT64 newLen);
	void fold(csr_t *shift);

	// Contestants can define their own functions below
};

/***********************************************************/
#endif
