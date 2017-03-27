#ifndef VITERBI_2D_OCL_H
#define VITERBI_2D_OCL_H

#include <vector>
#include <string>
#include <stdint.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <data_view.h>

#include "proxyCL.h"

namespace ublas = boost::numeric::ublas;

static const int8_t MOVE_DIAG  = 0;
static const int8_t MOVE_RIGHT = 1;
static const int8_t MOVE_UP    = 2;
static const int8_t MOVE_UNDEF = 3;
static const int32_t ZERO_PROB_SCORE = -1000000000;
static const double MIN_EMISSION_SCORE = -20.0;


inline double square(double x) {
  return x * x;
}

inline int32_t prob2score(double x) {
  if (x < 0.0000000001) return -2400;
  return int32_t(100.0 * log(x));
}

/// Struct representing all nodes in the HMM
struct HmmNodesData {
    /* The NodeSlice struct represents a diagonal slice through the nodes which can be
     * scheduled simultaneously as they only have data dependencies on the previous two
     * slices (the left and down nodes are in the previous slice, the down-left diagonal
     * node is in the slice before that) and none on nodes within the current slice.
     *
     * The diagram below shows the nodes within a narrow band around an alignment of
     * two base sequences. The first 10 diagonals of nodes represented by the digits 0-9,
     * further diagonals are indicated by backslashes.
     * For example, slice 6 is represented as follows:
     * - It contains 5 nodes (.numNodes = 5)
     * - The first (top-left) node is at position (2, 4) (.index1 = 2, .index2 = 4)
     * - The top-left node is node #18 (.firstNode = 18) as the previous
     *   slices have 1, 2, 3, 4, 4, and 4 nodes, respectively.
     * - The first node has no valid left node but a valid SW diagonal node
     *   (.firstLeftValid = false, .firstDiagonalValid = true)
     * - The last node has no valid diagonal nor down nodes
     *   (.lastDownValid = false, .lastDiagonalValid = false)
     *
     *    ^
     *    |     \\\\\\\\\
     * s  |    9\\\\\\\\\
     * e 4|  6789\\\\\\
     * q  |3456789\\\\\
     * 2  |23456789\\
     *    |12345678
     *   0|0123456
     *    *------------------->
     *     0 2   
     *         sequence1
     */
    struct NodeSlice {
        int32_t numNodes;                   // Number of nodes in this slice
        int32_t firstNode;                  // Index of first node in slice
        int32_t index1;                     // Position of first node in slice along sequence 1
        int32_t index2;                     // Position of first node in slice along sequence 2
        bool firstLeftValid;                // Does the first node in this slice have a left node?
        bool lastDownValid;                 // Does the last node in this slice have a down node?
        bool firstDiagonalValid;            // Does the first node in this slice have a diagonal node?
        bool lastDiagonalValid;             // Does the last node in this slice have a diagonal node?
    };
    int32_t maxSliceSize;                   // Size of largest slice
    std::vector<NodeSlice> slices;          // List of slices to be scheduled separately
    ublas::matrix<int16_t> statePointers;   // Viterbi backtrace pointers.
    ublas::matrix<int8_t> dirPointers;      // NW alignment backtrace pointers.
};


/** Helper class for emission scores.
 *
 *  This class provides normal level emissions and gamma distributed noise emissions.
 *  Note that other emission objects can be substituted by changing the Emission typedef
 *  immediately following this class definition.
 */
class DefaultEmission {
private:
  std::vector<double> levels;
  std::vector<double> noises;
  std::vector<double> logNoises;
  std::vector<double> stayWeights;
  std::vector<double> emWeights;
  std::vector<double> modelLevels;
  std::vector<double> modelNoises;
  std::vector<double> offsets;
  std::vector<double> levelScales;
  std::vector<double> noiseScales;
  std::vector<double> noiseShapes;
  int numEvents;
  int numStates;
  bool useNoise;

public:
  /** Constructor.
   *  @param[in] mdlLevels Model current levels.
   *  @param[in] mdlLevelSpreads Spreads of model current levels.
   *  @params[in] mdlNoises Model noise levels.
   *  @param[in] mdlNoiseSpreads Spreads of model noise levels.
   *  @param[in] useSd Flag to specify whether to use noise levels in the basecall.
   */
    DefaultEmission(const std::vector<double>& mdlLevels, const std::vector<double>& mdlLevelSpreads,
                  const std::vector<double>& mdlNoises, const std::vector<double>& mdlNoiseSpreads,
                  bool useSd);

  /** Assign events to the object with vectors.
   *  @param[in] means Event current levels.
   *  @param[in] stdvs Event noise levels.
   *  @param[in] stayWts Event weights for modifying stay probabilities.
   *  @param[in] emWts Event weights for modifying emission probabilities.
   */
  void SetEvents(const std::vector<double>& means, const std::vector<double>& stdvs,
      const std::vector<double>& stayWts, const std::vector<double>& emWts);

  /// Set the number of events (for when SetEvents() will not be called.
  void SetNEvents(int n) {numEvents = n;}

  /// Returns the number of events.
  int NumEvents() const {return numEvents;}

  /// Returns the number of model states.
  int NumStates() const {return numStates;}

  /// Returns the model levels.
  const std::vector<double> GetModelLevels() const { return modelLevels; }

  /// Returns the stay weights.
  const std::vector<double> GetStayWeights() const { return stayWeights; }

  /// Returns the score for event i and state j.
  int32_t Score(int i, int j) const {
    double score = offsets[j] + levelScales[j] * square(levels[i] - modelLevels[j]);
    if (useNoise) score += (noiseShapes[j] - 1.0) * logNoises[i] - noiseScales[j] * noises[i];
    return int32_t(emWeights[i] * std::max(MIN_EMISSION_SCORE, score));
  }
};


typedef DefaultEmission Emission;
typedef std::vector<std::pair<int32_t, int32_t> > Alignment;


/// Worker class for performing 2D Viterbi basecall.
class Viterbi2Docl {
private:
  proxyCL &proxy_cl_;
  HmmNodesData nodes;                              // All HMM nodes, in the order they should be processed.
  std::vector<double> transProbs;                  // Nine transition probabilities (stay * dir, step * dir, skip * dir).
  ublas::matrix<int32_t> emScore1;                 // Pre-computed emissions for sequence 1.
  ublas::matrix<int32_t> emScore2;                 // Pre-computed emissions for sequence 2.
  std::vector<int32_t> viterbiScore;               // Viterbi scores for last node.
  int numStates;                                   // Number of states in the HMM.
  int numNodes;                                    // Total number of nodes to be processed.
  int numEvents1;                                  // Number of events in sequence 1.
  int numEvents2;                                  // Number of events in sequence 2.
  bool enable_fp64_;                               // Whether to use double floating point
  cl::Kernel kernelProcessNodes;                   // OpenCL kernel objects
  cl::Kernel kernelPickBest;

  void initNodes(const std::vector<int32_t>& bandStarts, const std::vector<int32_t>& bandEnds);
  void processNodes(const std::vector<double>& wts1, const std::vector<double>& wts2,
    const std::vector<int32_t>& priors);
  void backTrace(Alignment& alignment, std::vector<int16_t>& states);

public:
  /** Constructor.
   */
  Viterbi2Docl(proxyCL& proxy_cl);

  /** Perform the basecall with emission objects.
   *  @param[in] data1 Emission object for sequence 1.
   *  @param[in] data2 Emission object for sequence 2.
   *  @param[in] bandStarts For each event in sequence 2, the first candidate position in sequence 1.
   *  @param[in] bandEnds For each event in sequence 2, the last candidate position in sequence 1.
   *  @param[in] priors The prior scores for the "before alignment" node. All zeros means no prior.
   *  @param[out] alignment The final alignment of events.
   *  @param[out] states The final basecalled states.
   */
  void Call(const Emission& data1, const Emission& data2, const std::vector<int32_t>& bandStarts,
            const std::vector<int32_t>& bandEnds, const std::vector<int32_t>& priors,
            Alignment& alignment, std::vector<int16_t>& states);

  /** Perform the basecall with precomputed emissions.
   *  @param[in] data1 Precomputed emissions for sequence 1.
   *  @param[in] data2 Precomputed emissions for sequence 2.
   *  @param[in] stayWt1 Stay weights for sequence 1.
   *  @param[in] stayWt2 Stay weights for sequence 2.
   *  @param[in] bandStarts For each event in sequence 2, the first candidate position in sequence 1.
   *  @param[in] bandEnds For each event in sequence 2, the last candidate position in sequence 1.
   *  @param[in] priors The prior scores for the "before alignment" node. All zeros means no prior.
   *  @param[out] alignment The final alignment of events.
   *  @param[out] states The final basecalled states.
   */
  void Call(const MatView<float>& data1, const MatView<float>& data2,
            const VecView<double>& stayWt1, const VecView<double>& stayWt2,
            const std::vector<int32_t>& bandStarts, const std::vector<int32_t>& bandEnds,
            const std::vector<int32_t>& priors, Alignment& alignment, std::vector<int16_t>& states);

  /* Set default transition values and allocate memory
   *  @param[in] len The maximum number of events to support for either sequence.
   *  @param[in] states The number of states in the HMM.
   *  @param[in] trans The six transition probabilities (stay1, step1, skip1, stay2, step2, skip2).
   */
  void InitData(int len, int states, const std::vector<double>& trans);

  bool InitCL(const std::string& srcKernelDir, const std::string& binKernelDir,
    std::string &error, bool enable_fp64, size_t num_states, size_t work_group_size = 0);
};

#endif /* VITERBI_2D_OCL_H */
