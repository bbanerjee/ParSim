#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/IntVector.h>
#include <gtest/gtest.h>

using namespace Uintah;

TEST(GridTest, Construction) {
  GridP grid = scinew Grid();
  EXPECT_EQ(grid->numLevels(), 0);
  
  Point anchor(0, 0, 0);
  Vector dcell(1, 1, 1);
  LevelP level = grid->addLevel(anchor, dcell);
  
  EXPECT_EQ(grid->numLevels(), 1);
  EXPECT_EQ(grid->getLevel(0), level);
}

TEST(LevelTest, Construction) {
  GridP grid = scinew Grid();
  Point anchor(0, 0, 0);
  Vector dcell(0.1, 0.1, 0.1);
  LevelP level = grid->addLevel(anchor, dcell);
  
  EXPECT_EQ(level->getIndex(), 0);
  EXPECT_EQ(level->getAnchor(), anchor);
  EXPECT_EQ(level->dCell(), dcell);
  EXPECT_DOUBLE_EQ(level->cellVolume(), 0.001);
}

TEST(LevelTest, AddPatch) {
  GridP grid = scinew Grid();
  LevelP level = grid->addLevel(Point(0,0,0), Vector(1,1,1));
  
  IntVector low(0, 0, 0);
  IntVector high(10, 10, 10);
  Patch* patch = level->addPatch(low, high, low, high, grid.get_rep());
  
  EXPECT_EQ(level->numPatches(), 1);
  EXPECT_EQ(level->getPatch(0), patch);
}
