#include <algorithm> // For std::ranges::fill (though assign is used for clarity here)
#include <gtest/gtest.h> // For Google Test framework
#include <unordered_map> // For std::unordered_map
#include <vector>        // For std::vector

// Define the Crack class
class Crack
{
public:
  // Public member vectors for simplicity in testing.
  // In a real application, these might be private with getter methods.
  std::vector<std::vector<int>> d_cfSegNodes;
  std::vector<std::vector<int>> d_cfSegPreIdx;
  std::vector<std::vector<int>> d_cfSegMinIdx;
  std::vector<std::vector<int>> d_cfSegMaxIdx;

  // Default constructor
  Crack() = default;

  /**
   * @brief Ensures the outer vectors (d_cfSegNodes, d_cfSegPreIdx, etc.) are
   * sized appropriately to contain the given segment index `m_idx`.
   * This prevents out-of-bounds access when accessing d_cfSegNodes[m_idx], etc.
   * @param m_idx The segment index to ensure capacity for.
   */
  void
  ensure_segment_vectors_sized(int m_idx)
  {
    if (m_idx < 0) {
      // Negative index is invalid; handle as an error or assert
      return;
    }
    // Resize all relevant member vectors if m_idx is beyond their current
    // capacity
    if (static_cast<size_t>(m_idx) >= d_cfSegNodes.size()) {
      d_cfSegNodes.resize(m_idx + 1);
    }
    if (static_cast<size_t>(m_idx) >= d_cfSegPreIdx.size()) {
      d_cfSegPreIdx.resize(m_idx + 1);
    }
    if (static_cast<size_t>(m_idx) >= d_cfSegMinIdx.size()) {
      d_cfSegMinIdx.resize(m_idx + 1);
    }
    if (static_cast<size_t>(m_idx) >= d_cfSegMaxIdx.size()) {
      d_cfSegMaxIdx.resize(m_idx + 1);
    }
  }

  /**
   * @brief Finds crack-front node indexes for a specific segment `m_idx`.
   * This includes:
   * 1. `d_cfSegPreIdx`: The previous index of a crack-front node
   * with the same value (if any, otherwise -1).
   * 2. `d_cfSegMinIdx` and `d_cfSegMaxIdx`: The minimum and maximum
   * node indexes of the crack-front segment to which the node belongs.
   * @param m_idx The index of the crack-front segment to process.
   */
  void
  FindCrackFrontNodeIndexes(const int& m_idx)
  {
    // Ensure the outer vectors are sufficiently sized for m_idx
    ensure_segment_vectors_sized(m_idx);
    if (m_idx < 0) {
      return;
    }

    // Get the number of nodes in the current segment
    int num = static_cast<int>(d_cfSegNodes[m_idx].size());

    // --- Part 1: Calculate d_cfSegPreIdx (Previous Node Index) ---
    // For each node, find the index of its last previous occurrence in the
    // segment. Optimized from O(N^2) to average O(N) using a hash map.
    // Initialize d_cfSegPreIdx[m_idx] with -1 and resize to `num`.
    d_cfSegPreIdx[m_idx].assign(num, -1);
    // Map to store the most recent index where each node value was encountered.
    std::unordered_map<int, int> last_seen_index;

    for (int i = 0; i < num; ++i) {
      int current_node_value = d_cfSegNodes[m_idx][i];
      // Check if this node value has been seen before (C++20 `contains`)
      if (last_seen_index.contains(current_node_value)) {
        // If yes, set d_cfSegPreIdx[m_idx][i] to the last seen index.
        d_cfSegPreIdx[m_idx][i] = last_seen_index[current_node_value];
      }
      // Always update the last seen index for the current node value to its
      // current position `i`.
      last_seen_index[current_node_value] = i;
    }

    // --- Part 2: Calculate d_cfSegMinIdx and d_cfSegMaxIdx (Segment Min/Max
    // Indexes) --- These identify the boundaries of logical crack-front
    // segments within d_cfSegNodes[m_idx]. Initialize both vectors with -1 and
    // resize to `num`.
    d_cfSegMinIdx[m_idx].assign(num, -1);
    d_cfSegMaxIdx[m_idx].assign(num, -1);

    // If the segment is empty, there's nothing more to do.
    if (num == 0) {
      return;
    }

    // `current_segment_min_idx` tracks the starting index of the segment
    // currently being defined.
    int current_segment_min_idx = 0;

    // Iterate through each node in the segment.
    // We determine segments block by block. Once a segment's min/max are found,
    // we fill values for all nodes within that segment.
    for (int i = 0; i < num; ++i) {
      // Only search for a new segment's boundaries if `i` is the actual start
      // of a new segment. This prevents redundant calculations for nodes
      // already assigned to a segment.
      if (i == current_segment_min_idx) {
        // `current_segment_max_idx` will store the end index of the current
        // segment. Initialize it to the end of the array, assuming it's a
        // single segment for now.
        int current_segment_max_idx = num - 1;

        // Determine `j_start_for_boundary_check` as per original logic:
        // If `i` is odd, `j_start` is `i`. If `i` is even, `j_start` is `i +
        // 1`. This means `j_start` will always be an odd index, or `num` if `i`
        // is `num-1` (and even).
        int j_start_for_boundary_check = (i % 2 != 0) ? i : (i + 1);

        // Iterate from `j_start_for_boundary_check` with a step of 2 to find
        // segment boundaries. This specific iteration pattern suggests that
        // segment breaks are determined by comparing nodes at `odd_index` and
        // `odd_index + 1`.
        for (int j = j_start_for_boundary_check; j < num; j += 2) {
          // Check if `j` is the last element. If so, the segment ends here.
          if (j == num - 1) {
            current_segment_max_idx = j;
            break; // End of array reached, segment ends here.
          }
          // Check if the node at `j` is different from the node at `j + 1`.
          // If they are different, it signifies a segment boundary.
          if (d_cfSegNodes[m_idx][j] != d_cfSegNodes[m_idx][j + 1]) {
            current_segment_max_idx = j;
            break; // Segment break detected.
          }
        }

        // After the inner loop, `current_segment_max_idx` holds the true end
        // index of the segment that started at `current_segment_min_idx`.

        // Now, assign these determined `minIdx` and `maxIdx` values to all
        // nodes within this identified segment.
        for (int k = current_segment_min_idx; k <= current_segment_max_idx;
             ++k) {
          d_cfSegMinIdx[m_idx][k] = current_segment_min_idx;
          d_cfSegMaxIdx[m_idx][k] = current_segment_max_idx;
        }

        // Advance `current_segment_min_idx` to the start of the next potential
        // segment.
        current_segment_min_idx = current_segment_max_idx + 1;
      }
      // The outer loop `i` will eventually catch up to
      // `current_segment_min_idx` to start processing the next block.
    }
  }

  // Determine how the crack-font nodes are connected
  void
  FindCrackFrontNodeIndexesOld(const int& m)
  {
    // Ensure the outer vectors are sufficiently sized for m_idx
    ensure_segment_vectors_sized(m);
    if (m < 0) {
      return;
    }
    
    // The previous node index of a crack-front node (d_cfSegPreIdx)
    // for which node[i]=node[preIdx] (preIdx<i)
    d_cfSegPreIdx[m].clear();
    int num = (int)d_cfSegNodes[m].size();
    d_cfSegPreIdx[m].resize(num);

    for (int i = 0; i < num; i++) {
      int preIdx   = -1;
      int thisNode = d_cfSegNodes[m][i];
      for (int j = i - 1; j >= 0; j--) {
        int preNode = d_cfSegNodes[m][j];
        if (thisNode == preNode) {
          preIdx = j;
          break;
        }
      }
      d_cfSegPreIdx[m][i] = preIdx;
    }

    // The minimum and maximum node indexes of the crack-front
    // on which the node resides: d_cfSegMaxIdx and d_cfSegMinIdx
    //d_cfSegMaxIdx[m].clear();
    //d_cfSegMinIdx[m].clear();
    d_cfSegMaxIdx[m].resize(num);
    d_cfSegMinIdx[m].resize(num);

    int maxIdx = -1, minIdx = 0;
    for (int i = 0; i < num; i++) {
      if (!(i >= minIdx && i <= maxIdx)) {
        for (int j = ((i % 2) != 0 ? i : i + 1); j < num; j += 2) {
          if (j == num - 1 ||
              (j < num - 1 && d_cfSegNodes[m][j] != d_cfSegNodes[m][j + 1])) {
            maxIdx = j;
            break;
          }
        }
      }
      d_cfSegMinIdx[m][i] = minIdx;
      d_cfSegMaxIdx[m][i] = maxIdx;
      if (i == maxIdx) {
        minIdx = maxIdx + 1;
      }
    }
  }
};

// --- Google Tests for Crack::FindCrackFrontNodeIndexes ---

// Test fixture for Crack class
struct CrackTest : public ::testing::Test
{
  Crack crack;
  Crack crack_old;
  // Common setup if needed for each test
  void
  SetUp() override
  {
    // Clear vectors before each test to ensure clean state
    crack.d_cfSegNodes.clear();
    crack.d_cfSegPreIdx.clear();
    crack.d_cfSegMinIdx.clear();
    crack.d_cfSegMaxIdx.clear();

    crack_old.d_cfSegNodes.clear();
    crack_old.d_cfSegPreIdx.clear();
    crack_old.d_cfSegMinIdx.clear();
    crack_old.d_cfSegMaxIdx.clear();
  }
};

// Test case for empty segment
TEST_F(CrackTest, HandlesEmptySegment)
{
  crack.d_cfSegNodes.resize(1); // Resize for segment 0
  crack.d_cfSegNodes[0] = {};   // Empty segment

  crack_old.d_cfSegNodes.resize(1); // Resize for segment 0
  crack_old.d_cfSegNodes[0] = {};   // Empty segment

  crack.FindCrackFrontNodeIndexes(0);
  crack_old.FindCrackFrontNodeIndexesOld(0); 

  ASSERT_EQ(crack.d_cfSegPreIdx[0].size(), 0);
  ASSERT_EQ(crack.d_cfSegMinIdx[0].size(), 0);
  ASSERT_EQ(crack.d_cfSegMaxIdx[0].size(), 0);

  ASSERT_EQ(crack_old.d_cfSegPreIdx[0].size(), 0);
  ASSERT_EQ(crack_old.d_cfSegMinIdx[0].size(), 0);
  ASSERT_EQ(crack_old.d_cfSegMaxIdx[0].size(), 0);
}

// Test case for a single node segment
TEST_F(CrackTest, HandlesSingleNodeSegment)
{
  crack.d_cfSegNodes.resize(1);
  crack.d_cfSegNodes[0] = { 100 };

  crack_old.d_cfSegNodes.resize(1);
  crack_old.d_cfSegNodes[0] = { 100 };

  crack.FindCrackFrontNodeIndexes(0);
  crack_old.FindCrackFrontNodeIndexesOld(0);

  ASSERT_EQ(crack.d_cfSegPreIdx[0], std::vector<int>{ -1 });
  ASSERT_EQ(crack.d_cfSegMinIdx[0], std::vector<int>{ 0 });
  ASSERT_EQ(crack.d_cfSegMaxIdx[0], std::vector<int>{ 0 });

  ASSERT_EQ(crack_old.d_cfSegPreIdx[0], std::vector<int>{ -1 });
  ASSERT_EQ(crack_old.d_cfSegMinIdx[0], std::vector<int>{ 0 });
  ASSERT_EQ(crack_old.d_cfSegMaxIdx[0], std::vector<int>{ 0 });
}

// --- Tests for d_cfSegPreIdx (Part 1) ---

TEST_F(CrackTest, PreIdx_NoDuplicates)
{
  crack.d_cfSegNodes.resize(1);
  crack.d_cfSegNodes[0] = { 1, 2, 3, 4, 5 };

  crack_old.d_cfSegNodes.resize(1);
  crack_old.d_cfSegNodes[0] = { 1, 2, 3, 4, 5 };

  crack.FindCrackFrontNodeIndexes(0);
  crack_old.FindCrackFrontNodeIndexesOld(0);

  ASSERT_EQ(crack.d_cfSegPreIdx[0], (std::vector<int>{ -1, -1, -1, -1, -1 }));
  ASSERT_EQ(crack_old.d_cfSegPreIdx[0], (std::vector<int>{ -1, -1, -1, -1, -1 }));
}

TEST_F(CrackTest, PreIdx_SimpleDuplicate)
{
  crack.d_cfSegNodes.resize(1);
  crack.d_cfSegNodes[0] = {
    1, 2, 1, 3, 2
  }; // 1 at index 0, 1 at index 2; 2 at index 1, 2 at index 4

  crack_old.d_cfSegNodes.resize(1);
  crack_old.d_cfSegNodes[0] = {
    1, 2, 1, 3, 2
  }; // 1 at index 0, 1 at index 2; 2 at index 1, 2 at index 4

  crack.FindCrackFrontNodeIndexes(0);
  crack_old.FindCrackFrontNodeIndexesOld(0);

  ASSERT_EQ(crack.d_cfSegPreIdx[0], (std::vector<int>{ -1, -1, 0, -1, 1 }));
  ASSERT_EQ(crack_old.d_cfSegPreIdx[0], (std::vector<int>{ -1, -1, 0, -1, 1 }));
}

TEST_F(CrackTest, PreIdx_MultipleDuplicates)
{
  crack.d_cfSegNodes.resize(1);
  crack.d_cfSegNodes[0] = { 5, 10, 5, 15, 10, 5 };

  crack_old.d_cfSegNodes.resize(1);
  crack_old.d_cfSegNodes[0] = { 5, 10, 5, 15, 10, 5 };

  crack.FindCrackFrontNodeIndexes(0);
  crack_old.FindCrackFrontNodeIndexesOld(0);

  // 5: -1, 0, 2
  // 10: -1, 1
  // 15: -1
  ASSERT_EQ(crack.d_cfSegPreIdx[0], (std::vector<int>{ -1, -1, 0, -1, 1, 2 }));
  ASSERT_EQ(crack_old.d_cfSegPreIdx[0], (std::vector<int>{ -1, -1, 0, -1, 1, 2 }));
}

TEST_F(CrackTest, PreIdx_AllSameValue)
{
  crack.d_cfSegNodes.resize(1);
  crack.d_cfSegNodes[0] = { 7, 7, 7, 7 };

  crack_old.d_cfSegNodes.resize(1);
  crack_old.d_cfSegNodes[0] = { 7, 7, 7, 7 };

  crack.FindCrackFrontNodeIndexes(0);
  crack_old.FindCrackFrontNodeIndexesOld(0);

  ASSERT_EQ(crack.d_cfSegPreIdx[0], (std::vector<int>{ -1, 0, 1, 2 }));
  ASSERT_EQ(crack_old.d_cfSegPreIdx[0], (std::vector<int>{ -1, 0, 1, 2 }));
}

TEST_F(CrackTest, PreIdx_LongSegmentMixedValues)
{
  crack.d_cfSegNodes.resize(1);
  crack.d_cfSegNodes[0] = { 10, 20, 10, 30, 40, 20, 10, 50 };

  crack_old.d_cfSegNodes.resize(1);
  crack_old.d_cfSegNodes[0] = { 10, 20, 10, 30, 40, 20, 10, 50 };

  crack.FindCrackFrontNodeIndexes(0);
  crack_old.FindCrackFrontNodeIndexesOld(0);

  // Expected:
  // 10: -1, 0, 2
  // 20: -1, 1
  // 30: -1
  // 40: -1
  // 50: -1
  ASSERT_EQ(crack.d_cfSegPreIdx[0],
            (std::vector<int>{ -1, -1, 0, -1, -1, 1, 2, -1 }));
  ASSERT_EQ(crack_old.d_cfSegPreIdx[0],
            (std::vector<int>{ -1, -1, 0, -1, -1, 1, 2, -1 }));
}

// --- Tests for d_cfSegMinIdx and d_cfSegMaxIdx (Part 2) ---

TEST_F(CrackTest, SegMinMax_SingleSegmentNoBreaks)
{
  crack.d_cfSegNodes.resize(1);
  crack.d_cfSegNodes[0] = {
    1, 1, 1, 1
  }; // No breaks using the `j` logic (odd indices)

  crack_old.d_cfSegNodes.resize(1);
  crack_old.d_cfSegNodes[0] = {
    1, 1, 1, 1
  }; // No breaks using the `j` logic (odd indices)

  crack.FindCrackFrontNodeIndexes(0);
  crack_old.FindCrackFrontNodeIndexesOld(0);

  // j_start for i=0 is 1. d_cfSegNodes[0][1]==d_cfSegNodes[0][2].
  // j_start for j=1 is 1.
  // j=1 (val 1) vs j+1=2 (val 1) -> equal. Continue.
  // j=3 (val 1) vs j+1=4 (out of bounds). j==num-1 (3==3). break.
  // current_segment_max_idx = 3. So all elements 0,1,2,3 should have min=0,
  // max=3.
  ASSERT_EQ(crack.d_cfSegMinIdx[0], (std::vector<int>{ 0, 0, 0, 0 }));
  ASSERT_EQ(crack.d_cfSegMaxIdx[0], (std::vector<int>{ 3, 3, 3, 3 }));
  ASSERT_EQ(crack_old.d_cfSegMinIdx[0], (std::vector<int>{ 0, 0, 0, 0 }));
  ASSERT_EQ(crack_old.d_cfSegMaxIdx[0], (std::vector<int>{ 3, 3, 3, 3 }));
}

TEST_F(CrackTest, SegMinMax_TwoPairedSegments)
{
  crack.d_cfSegNodes.resize(1);
  // Values:   10, 10, 20, 20, 30, 30
  // Indices:  0,  1,  2,  3,  4,  5
  crack.d_cfSegNodes[0] = { 10, 10, 20, 20, 30, 30 };

  crack_old.d_cfSegNodes.resize(1);
  // Values:   10, 10, 20, 20, 30, 30
  // Indices:  0,  1,  2,  3,  4,  5
  crack_old.d_cfSegNodes[0] = { 10, 10, 20, 20, 30, 30 };

  crack.FindCrackFrontNodeIndexes(0);
  crack_old.FindCrackFrontNodeIndexesOld(0);

  // Expected logic:
  // i=0: current_segment_min_idx = 0.
  //   j_start = 1.
  //   j=1 (node 10) vs j+1=2 (node 20). Different. current_segment_max_idx = 1.
  //   Fill indices 0..1: min=0, max=1.
  //   d_cfSegMinIdx = {0,0,...}
  //   d_cfSegMaxIdx = {1,1,...}
  //   current_segment_min_idx = 2.
  // i=1: Already processed (part of 0..1 segment).
  // i=2: current_segment_min_idx = 2.
  //   j_start = 3.
  //   j=3 (node 20) vs j+1=4 (node 30). Different. current_segment_max_idx = 3.
  //   Fill indices 2..3: min=2, max=3.
  //   d_cfSegMinIdx = {0,0,2,2,...}
  //   d_cfSegMaxIdx = {1,1,3,3,...}
  //   current_segment_min_idx = 4.
  // i=3: Already processed.
  // i=4: current_segment_min_idx = 4.
  //   j_start = 5.
  //   j=5 (node 30) vs j+1=6 (out of bounds). j==num-1 (5==5).
  //   current_segment_max_idx = 5. Fill indices 4..5: min=4, max=5.
  //   d_cfSegMinIdx = {...,4,4}
  //   d_cfSegMaxIdx = {...,5,5}
  //   current_segment_min_idx = 6.
  // i=5: Already processed.

  ASSERT_EQ(crack.d_cfSegMinIdx[0], (std::vector<int>{ 0, 0, 2, 2, 4, 4 }));
  ASSERT_EQ(crack.d_cfSegMaxIdx[0], (std::vector<int>{ 1, 1, 3, 3, 5, 5 }));
  ASSERT_EQ(crack_old.d_cfSegMinIdx[0], (std::vector<int>{ 0, 0, 2, 2, 4, 4 }));
  ASSERT_EQ(crack_old.d_cfSegMaxIdx[0], (std::vector<int>{ 1, 1, 3, 3, 5, 5 }));
}

TEST_F(CrackTest, SegMinMax_MixedSegmentsAndSinglesAtEnd)
{
  crack.d_cfSegNodes.resize(1);
  // Values:   10, 10, 20, 30, 40, 40, 50
  // Indices:  0,  1,  2,  3,  4,  5,  6
  crack.d_cfSegNodes[0] = { 10, 10, 20, 30, 40, 40, 50 };

  crack_old.d_cfSegNodes.resize(1);
  // Values:   10, 10, 20, 30, 40, 40, 50
  // Indices:  0,  1,  2,  3,  4,  5,  6
  crack_old.d_cfSegNodes[0] = { 10, 10, 20, 30, 40, 40, 50 };

  crack.FindCrackFrontNodeIndexes(0);
  crack_old.FindCrackFrontNodeIndexesOld(0);

  // Expected logic:
  // i=0: min=0. j_start=1. j=1(10) vs j+1=2(20). Diff. max=1. Fill 0..1 with
  // min=0, max=1. Next min=2. i=1: processed. i=2: min=2. j_start=3. j=3(30) vs
  // j+1=4(40). Diff. max=3. Fill 2..3 with min=2, max=3. Next min=4. i=3:
  // processed. i=4: min=4. j_start=5. j=5(40) vs j+1=6(50). Diff. max=5.
  // Fill 4..5 with min=4, max=5. Next min=6. i=5: processed. i=6: min=6.
  // j_start=7. j=7 (num). Max=num-1=6. Fill 6..6 with min=6, max=6. Next min=7.

  ASSERT_EQ(crack.d_cfSegMinIdx[0], (std::vector<int>{ 0, 0, 2, 2, 4, 4, 6 }));
  ASSERT_EQ(crack.d_cfSegMaxIdx[0], (std::vector<int>{ 1, 1, 3, 3, 5, 5, 6 }));
  ASSERT_EQ(crack_old.d_cfSegMinIdx[0], (std::vector<int>{ 0, 0, 2, 2, 4, 4, 6 }));
  ASSERT_EQ(crack_old.d_cfSegMaxIdx[0], (std::vector<int>{ 1, 1, 3, 3, 5, 5, 6 }));
}

TEST_F(CrackTest, SegMinMax_SingleElementBreaksThrough)
{
  crack.d_cfSegNodes.resize(1);
  // Values:   1, 2, 3, 4, 5
  // Indices:  0, 1, 2, 3, 4
  crack.d_cfSegNodes[0] = { 1, 2, 3, 4, 5 };

  crack_old.d_cfSegNodes.resize(1);
  // Values:   1, 2, 3, 4, 5
  // Indices:  0, 1, 2, 3, 4
  crack_old.d_cfSegNodes[0] = { 1, 2, 3, 4, 5 };

  crack.FindCrackFrontNodeIndexes(0);
  crack_old.FindCrackFrontNodeIndexesOld(0);

  // Expected logic for j = (i%2 != 0 ? i : i+1); j+=2:
  // i=0: min=0. j_start=1. j=1(2) vs j+1=2(3). Diff. max=1. Fill 0..1 with
  // min=0,max=1. Next min=2. i=1: processed. i=2: min=2. j_start=3. j=3(4) vs
  // j+1=4(5). Diff. max=3. Fill 2..3 with min=2,max=3. Next min=4. i=3:
  // processed. i=4: min=4. j_start=5. j=5 (num). max=num-1=4. Fill 4..4 with
  // min=4,max=4. Next min=5.

  ASSERT_EQ(crack.d_cfSegMinIdx[0], (std::vector<int>{ 0, 0, 2, 2, 4 }));
  ASSERT_EQ(crack.d_cfSegMaxIdx[0], (std::vector<int>{ 1, 1, 3, 3, 4 }));
  ASSERT_EQ(crack_old.d_cfSegMinIdx[0], (std::vector<int>{ 0, 0, 2, 2, 4 }));
  ASSERT_EQ(crack_old.d_cfSegMaxIdx[0], (std::vector<int>{ 1, 1, 3, 3, 4 }));
}

// Test with multiple segments (e.g., m=0 and m=1)
TEST_F(CrackTest, HandlesMultipleSegments)
{
  crack.d_cfSegNodes.resize(2);               // Two segments: m=0 and m=1
  crack.d_cfSegNodes[0] = { 10, 20, 10, 30 }; // Segment 0
  crack.d_cfSegNodes[1] = { 100, 100, 200, 200, 100 }; // Segment 1

  crack_old.d_cfSegNodes.resize(2);               // Two segments: m=0 and m=1
  crack_old.d_cfSegNodes[0] = { 10, 20, 10, 30 }; // Segment 0
  crack_old.d_cfSegNodes[1] = { 100, 100, 200, 200, 100 }; // Segment 1

  crack.FindCrackFrontNodeIndexes(0);
  crack.FindCrackFrontNodeIndexes(1);
  crack_old.FindCrackFrontNodeIndexesOld(0);
  crack_old.FindCrackFrontNodeIndexesOld(1);

  // Verify Segment 0 (10, 20, 10, 30)
  ASSERT_EQ(crack.d_cfSegPreIdx[0], (std::vector<int>{ -1, -1, 0, -1 }));
  // i=0: min=0. j_start=1. j=1(20) vs j+1=2(10). Diff. max=1. Fill 0..1. Next
  // min=2. i=1: processed. i=2: min=2. j_start=3. j=3(30) vs j+1=4(out). max=3.
  // Fill 2..3. Next min=4.
  ASSERT_EQ(crack.d_cfSegMinIdx[0], (std::vector<int>{ 0, 0, 2, 2 }));
  ASSERT_EQ(crack.d_cfSegMaxIdx[0], (std::vector<int>{ 1, 1, 3, 3 }));

  // Verify Segment 1 (100, 100, 200, 200, 100)
  ASSERT_EQ(
    crack.d_cfSegPreIdx[1],
    (std::vector<int>{
      -1, 0, -1, 2, 1 })); // For 100,100,200,200,100 -> P_idx: -1,0,-1,2,1
  // i=0: min=0. j_start=1. j=1(100) vs j+1=2(200). Diff. max=1. Fill 0..1. Next
  // min=2. i=1: processed. i=2: min=2. j_start=3. j=3(200) vs j+1=4(100). Diff.
  // max=3. Fill 2..3. Next min=4. i=3: processed. i=4: min=4. j_start=5 (out of
  // bounds). max=num-1=4. Fill 4..4. Next min=5.
  ASSERT_EQ(crack.d_cfSegMinIdx[1], (std::vector<int>{ 0, 0, 2, 2, 4 }));
  ASSERT_EQ(crack.d_cfSegMaxIdx[1], (std::vector<int>{ 1, 1, 3, 3, 4 }));
  // Verify Segment 0 (10, 20, 10, 30)
  ASSERT_EQ(crack_old.d_cfSegPreIdx[0], (std::vector<int>{ -1, -1, 0, -1 }));
  // i=0: min=0. j_start=1. j=1(20) vs j+1=2(10). Diff. max=1. Fill 0..1. Next
  // min=2. i=1: processed. i=2: min=2. j_start=3. j=3(30) vs j+1=4(out). max=3.
  // Fill 2..3. Next min=4.
  ASSERT_EQ(crack_old.d_cfSegMinIdx[0], (std::vector<int>{ 0, 0, 2, 2 }));
  ASSERT_EQ(crack_old.d_cfSegMaxIdx[0], (std::vector<int>{ 1, 1, 3, 3 }));

  // Verify Segment 1 (100, 100, 200, 200, 100)
  ASSERT_EQ(
    crack_old.d_cfSegPreIdx[1],
    (std::vector<int>{
      -1, 0, -1, 2, 1 })); // For 100,100,200,200,100 -> P_idx: -1,0,-1,2,1
  // i=0: min=0. j_start=1. j=1(100) vs j+1=2(200). Diff. max=1. Fill 0..1. Next
  // min=2. i=1: processed. i=2: min=2. j_start=3. j=3(200) vs j+1=4(100). Diff.
  // max=3. Fill 2..3. Next min=4. i=3: processed. i=4: min=4. j_start=5 (out of
  // bounds). max=num-1=4. Fill 4..4. Next min=5.
  ASSERT_EQ(crack_old.d_cfSegMinIdx[1], (std::vector<int>{ 0, 0, 2, 2, 4 }));
  ASSERT_EQ(crack_old.d_cfSegMaxIdx[1], (std::vector<int>{ 1, 1, 3, 3, 4 }));
}

// Test for edge case where j_start_for_boundary_check is beyond num
TEST_F(CrackTest, SegMinMax_JStartBeyondNum)
{
  crack.d_cfSegNodes.resize(1);
  crack.d_cfSegNodes[0] = { 1 }; // Single element, i=0, num=1

  crack_old.d_cfSegNodes.resize(1);
  crack_old.d_cfSegNodes[0] = { 1 }; // Single element, i=0, num=1

  crack.FindCrackFrontNodeIndexes(0);
  crack_old.FindCrackFrontNodeIndexesOld(0);

  // For i=0, j_start_for_boundary_check becomes 1. Loop `for (int j = 1; j < 1;
  // ...)` does not run. current_segment_max_idx remains num-1 (which is 0).
  // Correct.
  ASSERT_EQ(crack.d_cfSegMinIdx[0], std::vector<int>{ 0 });
  ASSERT_EQ(crack.d_cfSegMaxIdx[0], std::vector<int>{ 0 });
  ASSERT_EQ(crack_old.d_cfSegMinIdx[0], std::vector<int>{ 0 });
  ASSERT_EQ(crack_old.d_cfSegMaxIdx[0], std::vector<int>{ 0 });
}

// Test with invalid segment index m
TEST_F(CrackTest, HandlesInvalidSegmentIndex)
{
  // Calling with negative index
  crack.FindCrackFrontNodeIndexes(
    -1); // Should gracefully return without crashing.
  // No assertions on vector contents, as they might be empty or unchanged,
  // but the key is no crash/exception from the function itself due to
  // `ensure_segment_vectors_sized`.

  // Call with an index larger than any pre-existing size.
  // `ensure_segment_vectors_sized` handles resizing the outer vectors.
  crack.d_cfSegNodes.resize(1);
  crack.d_cfSegNodes[0] = { 1, 2 };
  crack_old.d_cfSegNodes.resize(1);
  crack_old.d_cfSegNodes[0] = { 1, 2 };
  crack.FindCrackFrontNodeIndexes(5); // This will resize vectors up to index 5
  crack_old.FindCrackFrontNodeIndexesOld(5);
  ASSERT_GE(crack.d_cfSegNodes.size(), 6);
  ASSERT_GE(crack.d_cfSegPreIdx.size(), 6);
  ASSERT_GE(crack.d_cfSegMinIdx.size(), 6);
  ASSERT_GE(crack.d_cfSegMaxIdx.size(), 6);
  // The newly created inner vector for index 5 will be empty.
  ASSERT_TRUE(crack.d_cfSegNodes[5].empty());
  // And its corresponding results should also be empty.
  ASSERT_TRUE(crack.d_cfSegPreIdx[5].empty());
  ASSERT_TRUE(crack.d_cfSegMinIdx[5].empty());
  ASSERT_TRUE(crack.d_cfSegMaxIdx[5].empty());
  ASSERT_GE(crack_old.d_cfSegNodes.size(), 6);
  ASSERT_GE(crack_old.d_cfSegPreIdx.size(), 6);
  ASSERT_GE(crack_old.d_cfSegMinIdx.size(), 6);
  ASSERT_GE(crack_old.d_cfSegMaxIdx.size(), 6);
  // The newly created inner vector for index 5 will be empty.
  ASSERT_TRUE(crack_old.d_cfSegNodes[5].empty());
  // And its corresponding results should also be empty.
  ASSERT_TRUE(crack_old.d_cfSegPreIdx[5].empty());
  ASSERT_TRUE(crack_old.d_cfSegMinIdx[5].empty());
  ASSERT_TRUE(crack_old.d_cfSegMaxIdx[5].empty());
}
