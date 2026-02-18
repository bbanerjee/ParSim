#include <Core/Containers/SuperBox.h>
#include <Core/Geometry/IntVector.h>
#include <gtest/gtest.h>
#include <vector>

using namespace Uintah;

// Mock Box class for testing SuperBox
class MockBox {
public:
  MockBox(IntVector low, IntVector high, int id) 
    : low_(low), high_(high), id_(id) {}
  
  IntVector getLow() const { return low_; }
  IntVector getHigh() const { return high_; }
  int getVolume() const {
    IntVector diff = high_ - low_;
    return diff.x() * diff.y() * diff.z();
  }
  int getVolume(IntVector low, IntVector high) const {
    IntVector diff = high - low;
    return diff.x() * diff.y() * diff.z();
  }
  int getID() const { return id_; }
  int getArea(int side) const {
    IntVector diff = high_ - low_;
    if (side == 0) return diff.y() * diff.z();
    if (side == 1) return diff.x() * diff.z();
    return diff.x() * diff.y();
  }

private:
  IntVector low_, high_;
  int id_;
};

using MockBoxP = MockBox*;

// Simple evaluator for testing
struct MockEvaluator {
  template<class BoxPIterator>
  int operator()(BoxPIterator begin, BoxPIterator end, 
    [[maybe_unused]] IntVector low, 
    [[maybe_unused]] IntVector high) {
    int count = 0;
    for (auto it = begin; it != end; ++it) count++;
    return count;
  }
};

using TestSuperBox = SuperBox<MockBoxP, IntVector, int, int, MockEvaluator>;
using TestBasicBox = BasicBox<MockBoxP, IntVector, int, int, MockEvaluator>;

TEST(SuperBoxTest, BasicBoxConstruction) {
  MockBox box(IntVector(0,0,0), IntVector(1,1,1), 1);
  TestBasicBox bb(&box);
  
  EXPECT_EQ(bb.getLow(), IntVector(0,0,0));
  EXPECT_EQ(bb.getHigh(), IntVector(1,1,1));
  EXPECT_EQ(bb.getVolume(), 1);
  EXPECT_EQ(bb.getBoxes().size(), 1u);
  EXPECT_EQ(bb.getBoxes()[0]->getID(), 1);
}

TEST(SuperBoxTest, RegionOverlaps) {
  TestSuperBox::Region r1(IntVector(0,0,0), IntVector(2,2,2));
  TestSuperBox::Region r2(IntVector(1,1,1), IntVector(3,3,3));
  TestSuperBox::Region r3(IntVector(3,3,3), IntVector(4,4,4));
  
  EXPECT_TRUE(r1.overlaps(r2));
  EXPECT_TRUE(r2.overlaps(r1));
  EXPECT_FALSE(r1.overlaps(r3));
}

TEST(SuperBoxTest, RegionContains) {
  TestSuperBox::Region r1(IntVector(0,0,0), IntVector(5,5,5));
  TestSuperBox::Region r2(IntVector(1,1,1), IntVector(2,2,2));
  
  EXPECT_TRUE(r1.contains(r2));
  EXPECT_FALSE(r2.contains(r1));
}
