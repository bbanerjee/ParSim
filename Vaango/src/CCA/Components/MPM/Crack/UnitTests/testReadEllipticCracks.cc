#include <cstdint>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <vector>

// Mock classes for testing (you'll need to replace these with your actual
// types)
struct Point
{
  double x, y, z;
  Point(double x = 0.0, double y = 0.0, double z = 0.0)
    : x(x)
    , y(y)
    , z(z)
  {
  }
  bool
  operator==(const Point& other) const
  {
    return x == other.x && y == other.y && z == other.z;
  }
};

/**
 * @brief Mock class to simulate the internal data structure of ProblemSpecP.
 * It holds configurations for "ellipse" blocks that the function will read.
 */
class MockProblemSpec
{
public:
  /**
   * @brief Represents a single "ellipse" block within the mock XML structure.
   * Contains data points, resolution, and an optional crack front segment ID.
   */
  struct EllipseBlock
  {
    std::vector<Point>
      points; // Expects 3 points: point1_axis1, point_axis2, point2_axis1
    std::optional<int> resolution_circumference; // Required field
    std::optional<int> crack_front_segment_ID;   // Optional field
  };

  std::vector<EllipseBlock> ellipseBlocks; // All configured ellipse blocks
  size_t currentBlockIdx =
    0; // Tracks the current "block" being processed during traversal

  /**
   * @brief Default constructor creating an empty (null-like) data object.
   */
  MockProblemSpec()
    : ellipseBlocks()
    , currentBlockIdx(0)
  {
  }

  /**
   * @brief Constructor to initialize the mock data with predefined ellipse
   * blocks.
   * @param blocks A vector of EllipseBlock objects.
   */
  explicit MockProblemSpec(std::vector<EllipseBlock> blocks)
    : ellipseBlocks(std::move(blocks))
    , currentBlockIdx(0)
  {
  }

  /**
   * @brief Simulates `geom_ps->findBlock("ellipse")`.
   * Returns a new shared_ptr pointing to a configured MockProblemSpec
   * representing the first "ellipse" block. Returns `nullptr` if no such block
   * is found.
   * @param name The name of the block to find (e.g., "ellipse").
   * @return A shared_ptr to a MockProblemSpecrepresenting the block, or
   * nullptr.
   */
  std::shared_ptr<MockProblemSpec>
  findBlock(const std::string& name) const
  {
    if (name == "ellipse" && !ellipseBlocks.empty()) {
      // Create a new shared_ptr that copies the current block data
      // but is configured to point to the first block.
      auto new_data = std::make_shared<MockProblemSpec>(ellipseBlocks);
      new_data->currentBlockIdx = 0; // Set to the first block
      return new_data;
    }
    return nullptr; // Simulates `0` (null pointer)
  }

  /**
   * @brief Simulates `ellipse_ps->findNextBlock("ellipse")`.
   * Returns a new shared_ptr pointing to a configured MockProblemSpec
   * representing the next "ellipse" block in the sequence. Returns `nullptr` if
   * no more blocks are found.
   * @param name The name of the block to find (e.g., "ellipse").
   * @return A shared_ptr to a MockProblemSpecrepresenting the next block, or
   * nullptr.
   */
  std::shared_ptr<MockProblemSpec>
  findNextBlock(const std::string& name) const
  {
    // Ensure we are looking for "ellipse" and there is a next block
    if (name == "ellipse" && currentBlockIdx + 1 < ellipseBlocks.size()) {
      // Create a new shared_ptr for the next block
      auto new_data = std::make_shared<MockProblemSpec>(ellipseBlocks);
      new_data->currentBlockIdx = currentBlockIdx + 1; // Move to the next block
      return new_data;
    }
    return nullptr; // Simulates `0` (no more blocks)
  }

  /**
   * @brief Simulates `ellipse_ps->require("key", p)` for Point type.
   * Retrieves a required Point value based on the key. Throws
   * `std::runtime_error` if the current block is invalid, or if the required
   * point data is missing/incomplete.
   * @param key The key of the point (e.g., "point1_axis1").
   * @param p Reference to a Point object to store the retrieved value.
   * @throws std::runtime_error if data is missing or key is unknown.
   */
  void
  require(const std::string& key, Point& p) const
  {
    if (currentBlockIdx >= ellipseBlocks.size()) {
      throw std::runtime_error(
        "Error: 'require' called on an invalid or out-of-bounds block index.");
    }
    const auto& current_block = ellipseBlocks[currentBlockIdx];
    if (current_block.points.size() < 3) {
      throw std::runtime_error("Error: Not enough points (expected 3) defined "
                               "in the current ellipse block for 'require'.");
    }

    if (key == "point1_axis1") {
      p = current_block.points[0];
    } else if (key == "point_axis2") {
      p = current_block.points[1];
    } else if (key == "point2_axis1") {
      p = current_block.points[2];
    } else {
      throw std::runtime_error("Error: Unknown required point key '" + key +
                               "'.");
    }
  }

  /**
   * @brief Simulates `ellipse_ps->require("key", n)` for int type.
   * Retrieves a required int value based on the key. Throws
   * `std::runtime_error` if the current block is invalid, or if the required
   * int data is missing.
   * @param key The key of the integer (e.g., "resolution_circumference").
   * @param n Reference to an int to store the retrieved value.
   * @throws std::runtime_error if data is missing or key is unknown.
   */
  void
  require(const std::string& key, int& n) const
  {
    if (currentBlockIdx >= ellipseBlocks.size()) {
      throw std::runtime_error(
        "Error: 'require' called on an invalid or out-of-bounds block index.");
    }
    const auto& current_block = ellipseBlocks[currentBlockIdx];
    if (key == "resolution_circumference") {
      if (current_block.resolution_circumference.has_value()) {
        n = current_block.resolution_circumference.value();
      } else {
        throw std::runtime_error("Error: Required 'resolution_circumference' "
                                 "not present in the current block.");
      }
    } else {
      throw std::runtime_error("Error: Unknown required int key '" + key +
                               "'.");
    }
  }

  /**
   * @brief Simulates `ellipse_ps->get("key", cfsID)` for optional int type.
   * Retrieves an optional int value based on the key. Returns `true` if the key
   * is found and present, `false` otherwise.
   * @param key The key of the integer (e.g., "crack_front_segment_ID").
   * @param cfsID Reference to an int to store the retrieved value.
   * @return `true` if the value was retrieved, `false` otherwise.
   */
  bool
  get(const std::string& key, int& cfsID) const
  {
    if (currentBlockIdx >= ellipseBlocks.size()) {
      return false; // Cannot get from an invalid block
    }
    const auto& current_block = ellipseBlocks[currentBlockIdx];
    if (key == "crack_front_segment_ID") {
      if (current_block.crack_front_segment_ID.has_value()) {
        cfsID = current_block.crack_front_segment_ID.value();
        return true;
      }
    }
    return false; // Key not found or not present
  }
};

// ProblemSpecP is defined as a shared pointer to our mock data.
// This mimics the smart pointer behavior implied by the original code's `!= 0`
// check.
using ProblemSpecP = std::shared_ptr<MockProblemSpec>;

// Mock Crack class for testing
class TestCrack
{
public:
  std::vector<std::vector<std::vector<Point>>> d_ellipses;
  std::vector<std::vector<int>> d_ellipseNCells;
  std::vector<std::vector<int>> d_ellipseCrkFrtSegID;

  void
  ReadEllipticCracks(const int& m, const ProblemSpecP& geom_ps);

  // Helper methods for testing
  void
  resizeForMaterial(int m)
  {
    if (m >= static_cast<int>(d_ellipses.size())) {
      d_ellipses.resize(m + 1);
      d_ellipseNCells.resize(m + 1);
      d_ellipseCrkFrtSegID.resize(m + 1);
    }
  }
};

void
TestCrack::ReadEllipticCracks(const int& m, const ProblemSpecP& geom_ps)
{
  // Ensure vectors are properly sized
  resizeForMaterial(m);

  for (ProblemSpecP ellipse_ps = geom_ps->findBlock("ellipse");
       ellipse_ps != nullptr;
       ellipse_ps = ellipse_ps->findNextBlock("ellipse")) {
    // Three points on the arc
    Point p;
    std::vector<Point> thisEllipsePts;
    ellipse_ps->require("point1_axis1", p);
    thisEllipsePts.push_back(p);
    ellipse_ps->require("point_axis2", p);
    thisEllipsePts.push_back(p);
    ellipse_ps->require("point2_axis1", p);
    thisEllipsePts.push_back(p);
    d_ellipses[m].push_back(thisEllipsePts);
    thisEllipsePts.clear();

    // Resolution on circumference
    int n = 1;
    ellipse_ps->require("resolution_circumference", n);
    d_ellipseNCells[m].push_back(n);

    // Crack front segment ID, -1 by default
    int cfsID;
    if (!ellipse_ps->get("crack_front_segment_ID", cfsID)) {
      cfsID = -1;
    }
    d_ellipseCrkFrtSegID[m].push_back(cfsID);
  }
}

TEST(CrackTest, HandlesNoEllipses)
{
  TestCrack crack_obj;
  // Create a ProblemSpecP that contains no ellipse blocks.
  ProblemSpecP geom_ps = std::make_shared<MockProblemSpec>();
  int m                = 0;

  // The function should run without throwing any exceptions and result in empty
  // data structures.
  crack_obj.ReadEllipticCracks(m, geom_ps);

  // Assert that the data structures for 'm' are empty.
  ASSERT_TRUE(crack_obj.d_ellipses[m].empty());
  ASSERT_TRUE(crack_obj.d_ellipseNCells[m].empty());
  ASSERT_TRUE(crack_obj.d_ellipseCrkFrtSegID[m].empty());
}

TEST(CrackTest, HandlesSingleEllipse)
{
  TestCrack crack_obj;
  int m = 10;

  std::vector<MockProblemSpec::EllipseBlock> blocks;
  blocks.push_back({
    { { 1.0, 2.0, 3.0 },
      { 4.0, 5.0, 6.0 },
      { 7.0, 8.0, 9.0 } },   // Three required Point objects
    std::optional<int>(100), // Required resolution_circumference
    std::optional<int>(5)    // Optional crack_front_segment_ID
  });
  ProblemSpecP geom_ps = std::make_shared<MockProblemSpec>(blocks);

  crack_obj.ReadEllipticCracks(m, geom_ps);

  // Verify that data was correctly parsed and stored.
  ASSERT_EQ(crack_obj.d_ellipses[m].size(), 1);
  ASSERT_EQ(crack_obj.d_ellipseNCells[m].size(), 1);
  ASSERT_EQ(crack_obj.d_ellipseCrkFrtSegID[m].size(), 1);

  ASSERT_EQ(crack_obj.d_ellipseNCells[m][0], 100);
  ASSERT_EQ(crack_obj.d_ellipseCrkFrtSegID[m][0], 5);
  ASSERT_EQ(crack_obj.d_ellipses[m][0][0].x, 1.0); // Check a point coordinate
}

TEST(CrackTest, HandlesSingleEllipseNoCfsId)
{
  TestCrack crack_obj;
  int m = 20;

  std::vector<MockProblemSpec::EllipseBlock> blocks;
  blocks.push_back({
    { { 10.0, 20.0, 30.0 },
      { 40.0, 50.0, 60.0 },
      { 70.0, 80.0, 90.0 } }, // Points
    std::optional<int>(200),  // Required resolution_circumference
    std::nullopt              // No optional crack_front_segment_ID
  });
  ProblemSpecP geom_ps = std::make_shared<MockProblemSpec>(blocks);

  crack_obj.ReadEllipticCracks(m, geom_ps);

  // Verify data and that optional field defaults to -1.
  ASSERT_EQ(crack_obj.d_ellipses[m].size(), 1);
  ASSERT_EQ(crack_obj.d_ellipseNCells[m].size(), 1);
  ASSERT_EQ(crack_obj.d_ellipseCrkFrtSegID[m].size(), 1);

  ASSERT_EQ(crack_obj.d_ellipseNCells[m][0], 200);
  ASSERT_EQ(crack_obj.d_ellipseCrkFrtSegID[m][0],
            -1); // Should default to -1 as per the function logic
}

TEST(CrackTest, ThrowsOnMissingRequiredResolution)
{
  TestCrack crack_obj;
  int m = 0;

  std::vector<MockProblemSpec::EllipseBlock> blocks;
  MockProblemSpec::EllipseBlock block_missing_res = {
    { { 1.0, 1.0, 1.0 }, { 2.0, 2.0, 2.0 }, { 3.0, 3.0, 3.0 } }, // Valid points
    std::nullopt, // Missing resolution_circumference
    std::nullopt  // Optional cfsID also missing
  };
  blocks.push_back(block_missing_res);
  ProblemSpecP geom_ps = std::make_shared<MockProblemSpec>(blocks);

  // Expect a std::runtime_error to be thrown because "resolution_circumference"
  // is required but marked as std::nullopt in our mock data.
  ASSERT_THROW(crack_obj.ReadEllipticCracks(m, geom_ps), std::runtime_error);
}