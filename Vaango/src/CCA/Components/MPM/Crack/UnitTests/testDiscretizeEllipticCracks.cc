#include <gtest/gtest.h> // For basic Google Test framework (optional for fuzzing itself)

#include <cmath> // For std::sqrt, std::cos, std::sin, M_PI
#include <format>
#include <iostream>
#include <map>
#include <memory>   // For std::shared_ptr
#include <optional> // For std::optional in C++17/C++20/C++23
#include <ranges>
#include <stdexcept> // For std::runtime_error
#include <string>
#include <vector>

// A small epsilon for floating-point comparisons to avoid issues with direct
// equality.
const double EPSILON = 1e-9;
const double PI      = 3.14159265358979323846; // Use a high-precision PI

// --- Mock Classes and Function under Test ---

/**
 * @brief Simple struct to represent a 3D point or vector.
 * Enhanced to include common vector operations.
 */
class Point
{
private:
  double x_, y_, z_;

public:
  // Default constructor (initializes to 0,0,0)
  Point()
    : x_(0.0)
    , y_(0.0)
    , z_(0.0)
  {
  }
  // Parameterized constructor
  Point(double x_val, double y_val, double z_val)
    : x_(x_val)
    , y_(y_val)
    , z_(z_val)
  {
  }

  // Methods to access components (for consistency with original code's x(),
  // y(), z())
  double
  x() const
  {
    return x_;
  }
  double
  y() const
  {
    return y_;
  }
  double
  z() const
  {
    return z_;
  }

  // Operator overloads for vector arithmetic
  Point
  operator+(const Point& other) const
  {
    return Point(x_ + other.x_, y_ + other.y_, z_ + other.z_);
  }
  Point
  operator-(const Point& other) const
  {
    return Point(x_ - other.x_, y_ - other.y_, z_ - other.z_);
  }
  Point
  operator*(double scalar) const
  {
    return Point(x_ * scalar, y_ * scalar, z_ * scalar);
  }

  // Division by scalar, handles division by zero by returning (0,0,0)
  Point
  operator/(double scalar) const
  {
    if (std::abs(scalar) < EPSILON) {
      // Handle division by zero by returning a zero vector
      return Point(0.0, 0.0, 0.0);
    }
    return Point(x_ / scalar, y_ / scalar, z_ / scalar);
  }

  // Dot product
  double
  dot(const Point& other) const
  {
    return x_ * other.x_ + y_ * other.y_ + z_ * other.z_;
  }

  // Cross product (used when Point acts as a Vector)
  Point
  cross(const Point& other) const
  {
    return Point(y_ * other.z_ - z_ * other.y_,
                 z_ * other.x_ - x_ * other.z_,
                 x_ * other.y_ - y_ * other.x_);
  }

  // Magnitude (length) of the vector
  double
  length() const
  {
    return std::sqrt(x_ * x_ + y_ * y_ + z_ * z_);
  }
  double
  lengthSq() const
  {
    return x_ * x_ + y_ * y_ + z_ * z_;
  }

  // Overload comparison operator for ease of testing in Google Tests
  // Uses EPSILON for floating-point comparison.
  bool
  operator==(const Point& other) const
  {
    return std::abs(x_ - other.x_) < EPSILON &&
           std::abs(y_ - other.y_) < EPSILON &&
           std::abs(z_ - other.z_) < EPSILON;
  }
};

// Define Vector as an alias for Point, as they behave similarly for 3D math.
// This preserves the original code's use of "Vector" while leveraging Point's
// functionality.
using Vector = Point;

/**
 * @brief Simple struct to represent an integer triplet for element
 * connectivity.
 */
struct IntVector
{
  int i, j, k;
  IntVector(int i_val, int j_val, int k_val)
    : i(i_val)
    , j(j_val)
    , k(k_val)
  {
  }

  // Overload comparison for testing
  bool
  operator==(const IntVector& other) const
  {
    return i == other.i && j == other.j && k == other.k;
  }
};

// --- Helper Functions for Crack Class ---

/**
 * @brief Calculates the direction cosine vector from point p1 to point p2.
 * @param p1 The starting point.
 * @param p2 The ending point.
 * @return A unit vector representing the direction from p1 to p2. Returns zero
 * vector if length is zero.
 */
Vector
TwoPtsDirCos(const Point& p1, const Point& p2)
{
  Vector v   = p2 - p1;
  double len = v.length();
  if (len < EPSILON) {
    return Vector(0.0, 0.0, 0.0); // Return zero vector if points are the same
  }
  return v / len;
}

/**
 * @brief Calculates the cross product of two vectors.
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return The cross product vector.
 */
Vector
Cross(const Vector& v1, const Vector& v2)
{
  return v1.cross(v2); // Use the cross method of the Point/Vector struct
}

/**
 * @brief Mock class to simulate the internal data structure of ProblemSpecP.
 * It holds configurations for "ellipse" blocks that the function will read.
 */
class MockProblemSpecData
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
  MockProblemSpecData()
    : ellipseBlocks()
    , currentBlockIdx(0)
  {
  }

  /**
   * @brief Constructor to initialize the mock data with predefined ellipse
   * blocks.
   * @param blocks A vector of EllipseBlock objects.
   */
  explicit MockProblemSpecData(std::vector<EllipseBlock> blocks)
    : ellipseBlocks(std::move(blocks))
    , currentBlockIdx(0)
  {
  }

  /**
   * @brief Simulates `geom_ps->findBlock("ellipse")`.
   * Returns a new shared_ptr pointing to a configured MockProblemSpecData
   * representing the first "ellipse" block. Returns `nullptr` if no such block
   * is found.
   * @param name The name of the block to find (e.g., "ellipse").
   * @return A shared_ptr to a MockProblemSpecData representing the block, or
   * nullptr.
   */
  std::shared_ptr<MockProblemSpecData>
  findBlock(const std::string& name) const
  {
    if (name == "ellipse" && !ellipseBlocks.empty()) {
      // Create a new shared_ptr that copies the current block data
      // but is configured to point to the first block.
      auto new_data = std::make_shared<MockProblemSpecData>(ellipseBlocks);
      new_data->currentBlockIdx = 0; // Set to the first block
      return new_data;
    }
    return nullptr; // Simulates `0` (null pointer)
  }

  /**
   * @brief Simulates `ellipse_ps->findNextBlock("ellipse")`.
   * Returns a new shared_ptr pointing to a configured MockProblemSpecData
   * representing the next "ellipse" block in the sequence. Returns `nullptr` if
   * no more blocks are found.
   * @param name The name of the block to find (e.g., "ellipse").
   * @return A shared_ptr to a MockProblemSpecData representing the next block,
   * or nullptr.
   */
  std::shared_ptr<MockProblemSpecData>
  findNextBlock(const std::string& name) const
  {
    // Ensure we are looking for "ellipse" and there is a next block
    if (name == "ellipse" && currentBlockIdx + 1 < ellipseBlocks.size()) {
      // Create a new shared_ptr for the next block
      auto new_data = std::make_shared<MockProblemSpecData>(ellipseBlocks);
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
using ProblemSpecP = std::shared_ptr<MockProblemSpecData>;

/**
 * @brief The Crack class containing the function to be fuzzed.
 * Members `d_ellipses`, `d_ellipseNCells`, and `d_ellipseCrkFrtSegID` are
 * implemented as maps to handle arbitrary `m` values.
 */
class Crack
{
public:
  // Existing members for ReadEllipticCracks
  std::map<int, std::vector<std::vector<Point>>> d_ellipses;
  std::map<int, std::vector<int>> d_ellipseNCells;
  std::map<int, std::vector<int>> d_ellipseCrkFrtSegID;

  // New members for DiscretizeEllipticCracks and
  // DiscretizePartialEllipticCracks
  std::map<int, std::vector<Point>> d_cx;       // Crack nodes
  std::map<int, std::vector<IntVector>> d_ce;   // Crack elements
  std::map<int, std::vector<int>> d_cfSegNodes; // Crack front segment nodes

  // New members for DiscretizePartialEllipticCracks
  std::map<int, std::vector<std::vector<Point>>>
    d_pellipses; // Partial ellipses, similar to d_ellipses
  std::map<int, std::vector<double>>
    d_pellipseExtent; // Extent for partial ellipses
  std::map<int, std::vector<int>>
    d_pellipseNCells; // Resolution for partial ellipses
  std::map<int, std::vector<int>>
    d_pellipseCrkFrtSegID; // Crack front segment ID for partial ellipses

  /**
   * @brief Reads elliptic crack data from the provided ProblemSpecP object.
   * This function is copied exactly from your prompt.
   * @param m An integer index.
   * @param geom_ps A ProblemSpecP object representing the geometric problem
   * specification.
   */
  void
  ReadEllipticCracks(const int& m, const ProblemSpecP& geom_ps)
  {
    // Loop through "ellipse" blocks
    // `ellipse_ps != 0` implicitly checks if the shared_ptr is not null.
    for (ProblemSpecP ellipse_ps = geom_ps->findBlock("ellipse");
         ellipse_ps != 0;
         ellipse_ps = ellipse_ps->findNextBlock("ellipse")) {

      // Three points on the arc (point1_axis1, point_axis2, point2_axis1)
      Point p;
      std::vector<Point> thisEllipsePts;
      // 'require' calls may throw std::runtime_error if the data is missing.
      ellipse_ps->require("point1_axis1", p);
      thisEllipsePts.push_back(p);
      ellipse_ps->require("point_axis2", p);
      thisEllipsePts.push_back(p);
      ellipse_ps->require("point2_axis1", p);
      thisEllipsePts.push_back(p);
      d_ellipses[m].push_back(thisEllipsePts);
      thisEllipsePts.clear(); // Important: clear points for the next ellipse

      // Resolution on circumference
      int n = 1; // Default value as per original code if not found by 'require'
      ellipse_ps->require("resolution_circumference", n);
      d_ellipseNCells[m].push_back(n);

      // Crack front segment ID, default to -1 if not specified
      int cfsID;
      if (!ellipse_ps->get("crack_front_segment_ID", cfsID)) {
        cfsID = -1; // Default value as per original code
      }
      d_ellipseCrkFrtSegID[m].push_back(cfsID);
    }
  }

  /**
   * @brief Discretizes elliptic cracks based on previously read data.
   * This function is copied exactly from your prompt.
   * @param m An integer index.
   * @param nnode0 Starting node ID, updated by the function.
   */
  void
  DiscretizeEllipticCracks(const int& m, int& nnode0)
  {
    // Ensure d_ellipses[m] exists to prevent accessing non-existent map entry
    if (d_ellipses.find(m) == d_ellipses.end()) {
      // Or handle as an error / log, depending on expected behavior.
      return;
    }

    for (int k = 0; k < (int)d_ellipses[m].size(); k++) {
      // Three points of the ellipse
      Point p1 = d_ellipses[m][k][0];
      Point p2 = d_ellipses[m][k][1];
      Point p3 = d_ellipses[m][k][2];

      // Center and half axial lengths of the ellipse
      double x0, y0, z0, a, b;
      Point origin = p3 + (p1 - p3) * 0.5;
      x0           = origin.x();
      y0           = origin.y();
      z0           = origin.z();
      a            = (p1 - origin).length();
      b            = (p2 - origin).length();

      // Local coordinates
      Vector v1, v2, v3;
      v1         = TwoPtsDirCos(origin, p1);
      v2         = TwoPtsDirCos(origin, p2);
      Vector v12 = Cross(v1, v2);

      // Handle case where v12 length is zero (e.g., v1 and v2 are collinear or
      // zero vectors)
      if (v12.length() < EPSILON) {
        v3 = Vector(0.0, 0.0, 0.0); // Or throw error/skip, depending on desired
                                    // behavior for invalid input.
      } else {
        v3 = v12 / v12.length();
      }

      double lx, mx, nx, ly, my, ny;
      lx = v1.x();
      mx = v1.y();
      nx = v1.z();
      ly = v2.x();
      my = v2.y();
      ny = v2.z();

      // Generate crack nodes
      d_cx[m].push_back(origin); // Center node
      for (int j = 0; j < d_ellipseNCells[m][k]; j++) {
        double thetai  = j * (2 * PI) / d_ellipseNCells[m][k];
        double xiprime = a * std::cos(thetai);
        double yiprime = b * std::sin(thetai);
        double xi      = lx * xiprime + ly * yiprime + x0;
        double yi      = mx * xiprime + my * yiprime + y0;
        double zi      = nx * xiprime + ny * yiprime + z0;
        d_cx[m].push_back(Point(xi, yi, zi));
      }

      // Generate crack elements
      for (int j = 1; j <= d_ellipseNCells[m][k]; j++) {
        int j1 = (j == d_ellipseNCells[m][k] ? 1 : j + 1);
        int n1 = nnode0;
        int n2 = nnode0 + j;
        int n3 = nnode0 + j1;
        d_ce[m].push_back(IntVector(n1, n2, n3));
        // Collect crack-front nodes
        if (d_ellipseCrkFrtSegID[m][k] == -1 ||
            d_ellipseCrkFrtSegID[m][k] == j) {
          d_cfSegNodes[m].push_back(n2);
          d_cfSegNodes[m].push_back(n3);
        }
      }
      nnode0 += d_ellipseNCells[m][k] + 1;
    }
  }

  /**
   * @brief Discretizes partial elliptic cracks based on previously read data.
   * This function is copied exactly from your prompt.
   * @param m An integer index.
   * @param nnode0 Starting node ID, updated by the function.
   */
  void
  DiscretizePartialEllipticCracks(const int& m, int& nnode0)
  {
    // Ensure d_pellipses[m] exists
    if (d_pellipses.find(m) == d_pellipses.end()) {
      return;
    }

    for (int k = 0; k < (int)d_pellipses[m].size(); k++) {
      double extent = d_pellipseExtent[m][k] / 360.;

      // Center, end points on major and minor axes
      Point origin  = d_pellipses[m][k][0];
      Point major_p = d_pellipses[m][k][1];
      Point minor_p = d_pellipses[m][k][2];
      double x0, y0, z0, a, b;
      x0 = origin.x();
      y0 = origin.y();
      z0 = origin.z();
      a  = (major_p - origin).length();
      b  = (minor_p - origin).length();

      // Local coordinates
      Vector v1, v2, v3;
      v1         = TwoPtsDirCos(origin, major_p);
      v2         = TwoPtsDirCos(origin, minor_p);
      Vector v12 = Cross(v1, v2);

      // Handle case where v12 length is zero
      if (v12.length() < EPSILON) {
        v3 = Vector(0.0, 0.0, 0.0); // Or throw error/skip
      } else {
        v3 = v12 / v12.length();
      }

      double lx, mx, nx, ly, my, ny;
      lx = v1.x();
      mx = v1.y();
      nx = v1.z();
      ly = v2.x();
      my = v2.y();
      ny = v2.z();

      // Generate crack nodes
      d_cx[m].push_back(origin); // Center node
      for (int j = 0; j <= d_pellipseNCells[m][k]; j++) {
        double thetai  = j * (2 * PI * extent) / d_pellipseNCells[m][k];
        double xiprime = a * std::cos(thetai);
        double yiprime = b * std::sin(thetai);
        double xi      = lx * xiprime + ly * yiprime + x0;
        double yi      = mx * xiprime + my * yiprime + y0;
        double zi      = nx * xiprime + ny * yiprime + z0;
        d_cx[m].push_back(Point(xi, yi, zi));
      }

      // Generate crack elements
      for (int j = 1; j <= d_pellipseNCells[m][k]; j++) {
        int n1 = nnode0;
        int n2 = nnode0 + j;
        int n3 = nnode0 + j + 1;
        d_ce[m].push_back(IntVector(n1, n2, n3));
        // Collect crack-front nodes
        if (d_pellipseCrkFrtSegID[m][k] == -1 ||
            d_pellipseCrkFrtSegID[m][k] == j) {
          d_cfSegNodes[m].push_back(n2);
          d_cfSegNodes[m].push_back(n3);
        }
      }
      nnode0 += d_pellipseNCells[m][k] + 2;
    }
  }
};

TEST(CrackTest, HandlesNoEllipses)
{
  Crack crack_obj;
  // Create a ProblemSpecP that contains no ellipse blocks.
  ProblemSpecP geom_ps = std::make_shared<MockProblemSpecData>();
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
  Crack crack_obj;
  int m = 10;

  std::vector<MockProblemSpecData::EllipseBlock> blocks;
  blocks.push_back({
    { { 1.0, 2.0, 3.0 },
      { 4.0, 5.0, 6.0 },
      { 7.0, 8.0, 9.0 } },   // Three required Point objects
    std::optional<int>(100), // Required resolution_circumference
    std::optional<int>(5)    // Optional crack_front_segment_ID
  });
  ProblemSpecP geom_ps = std::make_shared<MockProblemSpecData>(blocks);

  crack_obj.ReadEllipticCracks(m, geom_ps);

  // Verify that data was correctly parsed and stored.
  ASSERT_EQ(crack_obj.d_ellipses[m].size(), 1);
  ASSERT_EQ(crack_obj.d_ellipseNCells[m].size(), 1);
  ASSERT_EQ(crack_obj.d_ellipseCrkFrtSegID[m].size(), 1);

  ASSERT_EQ(crack_obj.d_ellipseNCells[m][0], 100);
  ASSERT_EQ(crack_obj.d_ellipseCrkFrtSegID[m][0], 5);
  ASSERT_EQ(crack_obj.d_ellipses[m][0][0].x(), 1.0); // Check a point coordinate
}

TEST(CrackTest, HandlesSingleEllipseNoCfsId)
{
  Crack crack_obj;
  int m = 20;

  std::vector<MockProblemSpecData::EllipseBlock> blocks;
  blocks.push_back({
    { { 10.0, 20.0, 30.0 },
      { 40.0, 50.0, 60.0 },
      { 70.0, 80.0, 90.0 } }, // Points
    std::optional<int>(200),  // Required resolution_circumference
    std::nullopt              // No optional crack_front_segment_ID
  });
  ProblemSpecP geom_ps = std::make_shared<MockProblemSpecData>(blocks);

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
  Crack crack_obj;
  int m = 0;

  std::vector<MockProblemSpecData::EllipseBlock> blocks;
  MockProblemSpecData::EllipseBlock block_missing_res = {
    { { 1.0, 1.0, 1.0 }, { 2.0, 2.0, 2.0 }, { 3.0, 3.0, 3.0 } }, // Valid points
    std::nullopt, // Missing resolution_circumference
    std::nullopt  // Optional cfsID also missing
  };
  blocks.push_back(block_missing_res);
  ProblemSpecP geom_ps = std::make_shared<MockProblemSpecData>(blocks);

  // Expect a std::runtime_error to be thrown because "resolution_circumference"
  // is required but marked as std::nullopt in our mock data.
  ASSERT_THROW(crack_obj.ReadEllipticCracks(m, geom_ps), std::runtime_error);
}

// --- Google Tests for DiscretizeEllipticCracks ---
namespace DiscretizeEllipticCracksTest {

TEST(DiscretizeEllipticCracksTest, CircleAroundOrigin)
{
  Crack crack_obj;
  int m      = 0;
  int nnode0 = 100; // Starting node ID

  // Define a full circle (p1, p2, p3 form an axis-aligned circle)
  // p1 = (a, 0, 0), p2 = (0, b, 0), p3 = (-a, 0, 0) for origin (0,0,0) and a=b
  double radius = 5.0;
  int resolution =
    4; // Generate 4 segments (4 nodes on circumference + 1 center)

  crack_obj.d_ellipses[m].push_back({ Point(radius, 0.0, 0.0),
                                      Point(0.0, radius, 0.0),
                                      Point(-radius, 0.0, 0.0) });
  crack_obj.d_ellipseNCells[m].push_back(resolution);
  crack_obj.d_ellipseCrkFrtSegID[m].push_back(
    -1); // All segments are crack front

  crack_obj.DiscretizeEllipticCracks(m, nnode0);

  // Expected nodes: center (0,0,0) + 4 points on circumference
  ASSERT_EQ(crack_obj.d_cx[m].size(),
            resolution + 1); // 1 center + 4 circumference points
  ASSERT_EQ(crack_obj.d_cx[m][0],
            Point(0.0, 0.0, 0.0)); // Center should be origin

  // Check some specific points for a circle with resolution 4
  // Nodes should be (0,0,0), (5,0,0), (0,5,0), (-5,0,0), (0,-5,0)
  ASSERT_EQ(crack_obj.d_cx[m][1], Point(radius, 0.0, 0.0));
  ASSERT_EQ(crack_obj.d_cx[m][2], Point(0.0, radius, 0.0));
  ASSERT_EQ(crack_obj.d_cx[m][3], Point(-radius, 0.0, 0.0));
  ASSERT_EQ(crack_obj.d_cx[m][4], Point(0.0, -radius, 0.0));

  // Expected elements: 4 triangles, all starting from center node (nnode0)
  ASSERT_EQ(crack_obj.d_ce[m].size(), resolution);
  ASSERT_EQ(crack_obj.d_ce[m][0], IntVector(100, 101, 102));
  ASSERT_EQ(crack_obj.d_ce[m][1], IntVector(100, 102, 103));
  ASSERT_EQ(crack_obj.d_ce[m][2], IntVector(100, 103, 104));
  ASSERT_EQ(crack_obj.d_ce[m][3],
            IntVector(100, 104, 101)); // Last element connects back to first

  // Expected crack front nodes (all segments)
  ASSERT_EQ(crack_obj.d_cfSegNodes[m].size(), resolution * 2);
  ASSERT_EQ(crack_obj.d_cfSegNodes[m][0], 101); // n2 for first element
  ASSERT_EQ(crack_obj.d_cfSegNodes[m][1], 102); // n3 for first element
  ASSERT_EQ(crack_obj.d_cfSegNodes[m][crack_obj.d_cfSegNodes[m].size() - 2],
            104);
  ASSERT_EQ(crack_obj.d_cfSegNodes[m][crack_obj.d_cfSegNodes[m].size() - 1],
            101);

  // Check final nnode0
  ASSERT_EQ(nnode0, 100 + resolution + 1);
}

TEST(DiscretizeEllipticCracksTest, DefaultFractureMPMTest)
{
  Crack crack_obj;
  int m          = 0;
  int resolution = 32;

  std::vector<MockProblemSpecData::EllipseBlock> blocks;
  blocks.push_back({
    { { 8.65e-3, 0.0, 4.98e-3 },
      { 0.0, 9.98e-3, 0.0 },
      { -8.65e-3, 0.0, -4.98e-3 } }, // Points
    std::optional<int>(resolution),  // Required resolution_circumference
    std::nullopt                     // No optional crack_front_segment_ID
  });
  ProblemSpecP geom_ps = std::make_shared<MockProblemSpecData>(blocks);

  crack_obj.ReadEllipticCracks(m, geom_ps);

  int nnode0 = 100; // Starting node ID

  crack_obj.DiscretizeEllipticCracks(m, nnode0);

  for (auto const& [idx, cx] :
       std::ranges::views::enumerate(crack_obj.d_cx[m])) {
    std::cout << std::format(
                   "idx: {} cx: {},{},{}", idx, cx.x(), cx.y(), cx.z())
              << std::endl;
  }

  for (auto const& [idx, ce] :
       std::ranges::views::enumerate(crack_obj.d_ce[m])) {
    std::cout << std::format("idx: {} ce: {},{},{}", idx, ce.i, ce.j, ce.k)
              << std::endl;
  }

  for (auto const& [idx, sn] :
       std::ranges::views::enumerate(crack_obj.d_cfSegNodes[m])) {
    std::cout << std::format("idx: {} cfSegNodes: {}", idx, sn)
              << std::endl;
  }

  // Expected nodes: center (0,0,0) + 32 points on circumference
  ASSERT_EQ(crack_obj.d_cx[m].size(),
            resolution + 1); // 1 center + circumference points
  ASSERT_EQ(crack_obj.d_cx[m][0],
            Point(0.0, 0.0, 0.0)); // Center should be origin

  // Expected elements: 32 triangles, all starting from center node (nnode0)
  ASSERT_EQ(crack_obj.d_ce[m].size(), resolution);
  ASSERT_EQ(crack_obj.d_ce[m][0], IntVector(100, 101, 102));
  ASSERT_EQ(crack_obj.d_ce[m][1], IntVector(100, 102, 103));
  ASSERT_EQ(crack_obj.d_ce[m][2], IntVector(100, 103, 104));
  ASSERT_EQ(crack_obj.d_ce[m][31],
            IntVector(100, 132, 101)); // Last element connects back to first

  // Expected crack front nodes (all segments)
  ASSERT_EQ(crack_obj.d_cfSegNodes[m].size(), resolution * 2);
  ASSERT_EQ(crack_obj.d_cfSegNodes[m][0], 101); // n2 for first element
  ASSERT_EQ(crack_obj.d_cfSegNodes[m][1], 102); // n3 for first element
  ASSERT_EQ(crack_obj.d_cfSegNodes[m][crack_obj.d_cfSegNodes[m].size() - 2],
            132);
  ASSERT_EQ(crack_obj.d_cfSegNodes[m][crack_obj.d_cfSegNodes[m].size() - 1],
            101);

  // Check final nnode0
  ASSERT_EQ(nnode0, 100 + resolution + 1);
}

TEST(DiscretizeEllipticCracksTest, EllipseWithOffset)
{
  Crack crack_obj;
  int m      = 1;
  int nnode0 = 0; // Starting node ID

  // An ellipse centered at (10, 10, 0) with a=4, b=2 on XY plane
  // p1 = (14, 10, 0), p2 = (10, 12, 0), p3 = (6, 10, 0)
  Point center(10.0, 10.0, 0.0);
  double a_len   = 4.0; // Half-axis length along X
  double b_len   = 2.0; // Half-axis length along Y
  int resolution = 8;

  crack_obj.d_ellipses[m].push_back(
    { Point(center.x() + a_len, center.y(), center.z()),
      Point(center.x(), center.y() + b_len, center.z()),
      Point(center.x() - a_len, center.y(), center.z()) });
  crack_obj.d_ellipseNCells[m].push_back(resolution);
  crack_obj.d_ellipseCrkFrtSegID[m].push_back(-1); // All crack front

  crack_obj.DiscretizeEllipticCracks(m, nnode0);

  // Verify center node
  ASSERT_EQ(crack_obj.d_cx[m][0], center);
  ASSERT_EQ(crack_obj.d_cx[m].size(), resolution + 1);

  // Check first node on circumference (at theta = 0)
  ASSERT_EQ(crack_obj.d_cx[m][1],
            Point(center.x() + a_len, center.y(), center.z()));

  // Check element count
  ASSERT_EQ(crack_obj.d_ce[m].size(), resolution);

  // Check final nnode0
  ASSERT_EQ(nnode0, resolution + 1);
}

TEST(DiscretizeEllipticCracksTest, CrackFrontSegments)
{
  Crack crack_obj;
  int m      = 2;
  int nnode0 = 50;

  // Ellipse 1: Specific segment as crack front (segment 2)
  crack_obj.d_ellipses[m].push_back(
    { Point(1, 0, 0), Point(0, 1, 0), Point(-1, 0, 0) });
  crack_obj.d_ellipseNCells[m].push_back(4);
  crack_obj.d_ellipseCrkFrtSegID[m].push_back(
    2); // Only segment 2 is crack front

  // Ellipse 2: All segments as crack front (-1)
  crack_obj.d_ellipses[m].push_back(
    { Point(10, 0, 0), Point(0, 10, 0), Point(-10, 0, 0) });
  crack_obj.d_ellipseNCells[m].push_back(2);
  crack_obj.d_ellipseCrkFrtSegID[m].push_back(-1); // All crack front

  crack_obj.DiscretizeEllipticCracks(m, nnode0);

  // After Ellipse 1 (resolution 4): nnode0 becomes 50 + 4 + 1 = 55
  // Nodes for Ellipse 1: 50(center), 51, 52, 53, 54
  // Elements: (50,51,52), (50,52,53), (50,53,54), (50,54,51)
  // Segment 2: Element 2 is (50,52,53) -> nodes 52, 53 are crack front
  ASSERT_EQ(crack_obj.d_cfSegNodes[m].size(),
            2 + (2 * 2)); // 2 from first ellipse, 4 from second
  ASSERT_EQ(crack_obj.d_cfSegNodes[m][0], 52);
  ASSERT_EQ(crack_obj.d_cfSegNodes[m][1], 53);

  // After Ellipse 2 (resolution 2): nnode0 becomes 55 + 2 + 1 = 58
  // Nodes for Ellipse 2: 55(center), 56, 57
  // Elements: (55,56,57), (55,57,56)
  // Both segments are crack front: -> nodes 56, 57, 57, 56 (or 56, 57)
  ASSERT_EQ(crack_obj.d_cfSegNodes[m][2], 56);
  ASSERT_EQ(crack_obj.d_cfSegNodes[m][3], 57);
  ASSERT_EQ(crack_obj.d_cfSegNodes[m][4], 57);
  ASSERT_EQ(crack_obj.d_cfSegNodes[m][5], 56);

  ASSERT_EQ(nnode0, 50 + (4 + 1) + (2 + 1)); // Total nnode0 update
}

TEST(DiscretizeEllipticCracksTest, ZeroLengthPoints)
{
  Crack crack_obj;
  int m      = 3;
  int nnode0 = 0;

  // Test case where points lead to zero length vectors (collinear points)
  crack_obj.d_ellipses[m].push_back(
    { Point(0, 0, 0), Point(0, 0, 0), Point(0, 0, 0) }); // All points same
  crack_obj.d_ellipseNCells[m].push_back(10);
  crack_obj.d_ellipseCrkFrtSegID[m].push_back(-1);

  // This should run without crashing, but the output nodes/elements will be
  // degenerate. The EPSILON checks in Point/Vector operations should prevent
  // division by zero.
  crack_obj.DiscretizeEllipticCracks(m, nnode0);

  // Expect the center to be (0,0,0) and all generated points to also be (0,0,0)
  ASSERT_EQ(crack_obj.d_cx[m][0], Point(0, 0, 0));
  // All subsequent points should also be (0,0,0) due to zero lengths 'a' and
  // 'b'
  for (size_t i = 1; i < crack_obj.d_cx[m].size(); ++i) {
    ASSERT_EQ(crack_obj.d_cx[m][i], Point(0, 0, 0));
  }
  ASSERT_EQ(crack_obj.d_cx[m].size(),
            10 + 1); // Still creates correct number of nodes, just at origin

  // Elements will be created but degenerate (all nodes point to the same few
  // IDs)
  ASSERT_EQ(crack_obj.d_ce[m].size(), 10);
  // nnode0 should still increment correctly
  ASSERT_EQ(nnode0, 10 + 1);
}

} // namespace DiscretizeEllipticCracksTest

// --- Google Tests for DiscretizePartialEllipticCracks ---
namespace DiscretizePartialEllipticCracksTest {

TEST(DiscretizePartialEllipticCracksTest, HalfCircle)
{
  Crack crack_obj;
  int m      = 0;
  int nnode0 = 1000;

  // A half-circle centered at (0,0,0), major axis along X, minor along Y
  // Points: origin(0,0,0), major_p(5,0,0), minor_p(0,5,0)
  double radius  = 5.0;
  double extent  = 180.0; // Half circle
  int resolution = 4;     // 4 segments for 180 degrees (5 points on arc)

  crack_obj.d_pellipses[m].push_back(
    { Point(0, 0, 0), Point(radius, 0, 0), Point(0, radius, 0) });
  crack_obj.d_pellipseExtent[m].push_back(extent);
  crack_obj.d_pellipseNCells[m].push_back(resolution);
  crack_obj.d_pellipseCrkFrtSegID[m].push_back(-1); // All crack front

  crack_obj.DiscretizePartialEllipticCracks(m, nnode0);

  // Expected nodes: center (0,0,0) + (resolution + 1) points on circumference
  ASSERT_EQ(crack_obj.d_cx[m].size(), 1 + (resolution + 1));
  ASSERT_EQ(crack_obj.d_cx[m][0],
            Point(0.0, 0.0, 0.0)); // Center should be origin

  // Check first and last generated points
  ASSERT_EQ(crack_obj.d_cx[m][1], Point(radius, 0.0, 0.0)); // j=0, thetai=0
  ASSERT_EQ(crack_obj.d_cx[m][crack_obj.d_cx[m].size() - 1],
            Point(-radius, 0.0, 0.0)); // j=resolution, thetai=PI

  // Expected elements: resolution triangles
  ASSERT_EQ(crack_obj.d_ce[m].size(), resolution);
  ASSERT_EQ(crack_obj.d_ce[m][0], IntVector(1000, 1001, 1002));
  ASSERT_EQ(crack_obj.d_ce[m][resolution - 1],
            IntVector(1000, 1000 + resolution, 1000 + resolution + 1));

  // Expected crack front nodes (all segments)
  ASSERT_EQ(crack_obj.d_cfSegNodes[m].size(), resolution * 2);

  // Check final nnode0
  ASSERT_EQ(nnode0, 1000 + resolution + 2);
}

TEST(DiscretizePartialEllipticCracksTest, QuarterEllipse)
{
  Crack crack_obj;
  int m      = 1;
  int nnode0 = 2000;

  // A quarter ellipse in first quadrant, centered at (0,0,0)
  // major_p = (4,0,0), minor_p = (0,2,0)
  double a_len   = 4.0;
  double b_len   = 2.0;
  double extent  = 90.0; // Quarter circle
  int resolution = 2;    // 2 segments for 90 degrees (3 points on arc)

  crack_obj.d_pellipses[m].push_back(
    { Point(0, 0, 0), Point(a_len, 0, 0), Point(0, b_len, 0) });
  crack_obj.d_pellipseExtent[m].push_back(extent);
  crack_obj.d_pellipseNCells[m].push_back(resolution);
  crack_obj.d_pellipseCrkFrtSegID[m].push_back(
    1); // Only segment 1 is crack front

  crack_obj.DiscretizePartialEllipticCracks(m, nnode0);

  // Expected nodes: center + resolution + 1 points on circumference
  ASSERT_EQ(crack_obj.d_cx[m].size(), 1 + (resolution + 1));
  ASSERT_EQ(crack_obj.d_cx[m][0], Point(0.0, 0.0, 0.0));

  // Check first and last generated points
  ASSERT_EQ(crack_obj.d_cx[m][1], Point(a_len, 0.0, 0.0)); // j=0, thetai=0
  ASSERT_TRUE(crack_obj.d_cx[m][3] ==
              Point(0.0, b_len, 0.0)); // j=resolution, thetai=PI/2

  // Elements
  ASSERT_EQ(crack_obj.d_ce[m].size(), resolution);
  ASSERT_EQ(crack_obj.d_ce[m][0], IntVector(2000, 2001, 2002));
  ASSERT_EQ(crack_obj.d_ce[m][1], IntVector(2000, 2002, 2003));

  // Crack front nodes (only segment 1, which corresponds to element with
  // n2=2001, n3=2002)
  ASSERT_EQ(crack_obj.d_cfSegNodes[m].size(), 2);
  ASSERT_EQ(crack_obj.d_cfSegNodes[m][0], 2001); // n2 for segment 1
  ASSERT_EQ(crack_obj.d_cfSegNodes[m][1], 2002); // n3 for segment 1

  // Check final nnode0
  ASSERT_EQ(nnode0, 2000 + resolution + 2);
}

TEST(DiscretizeEllipticCracksTest, EmptyEllipses)
{
  Crack crack_obj;
  int m      = 5;
  int nnode0 = 10;

  // Call DiscretizeEllipticCracks with empty d_ellipses[m]
  crack_obj.DiscretizeEllipticCracks(m, nnode0);
  ASSERT_TRUE(crack_obj.d_cx[m].empty());
  ASSERT_TRUE(crack_obj.d_ce[m].empty());
  ASSERT_TRUE(crack_obj.d_cfSegNodes[m].empty());
  ASSERT_EQ(nnode0, 10); // nnode0 should not change if no ellipses processed

  // Call DiscretizePartialEllipticCracks with empty d_pellipses[m]
  crack_obj.DiscretizePartialEllipticCracks(m, nnode0);
  ASSERT_TRUE(
    crack_obj.d_cx[m]
      .empty()); // Should still be empty or cleared if previous was not empty
  ASSERT_TRUE(crack_obj.d_ce[m].empty());
  ASSERT_TRUE(crack_obj.d_cfSegNodes[m].empty());
  ASSERT_EQ(nnode0, 10); // nnode0 should not change

  // Add a non-empty m for d_ellipses to ensure default map behavior is not an
  // issue
  crack_obj.d_ellipses[99].push_back(
    { Point(1, 0, 0), Point(0, 1, 0), Point(-1, 0, 0) });
  crack_obj.d_ellipseNCells[99].push_back(4);
  crack_obj.d_ellipseCrkFrtSegID[99].push_back(-1);

  int another_nnode0 = 0;
  crack_obj.DiscretizeEllipticCracks(99, another_nnode0);
  ASSERT_FALSE(crack_obj.d_cx[99].empty());
  ASSERT_FALSE(crack_obj.d_ce[99].empty());
  ASSERT_FALSE(crack_obj.d_cfSegNodes[99].empty());
  ASSERT_EQ(another_nnode0, 4 + 1); // Should be updated
}

} // namespace DiscretizePartialEllipticCracksTest
