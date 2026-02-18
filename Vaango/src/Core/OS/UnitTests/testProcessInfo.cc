#include <Core/OS/ProcessInfo.h>
#include <gtest/gtest.h>
#include <string>

using namespace Uintah;

TEST(ProcessInfoTest, BasicFunctions) {
  if (ProcessInfo::isSupported(ProcessInfo::MEM_SIZE)) {
    unsigned long mem = ProcessInfo::getMemoryUsed();
    EXPECT_GT(mem, 0u);
    
    std::string human = ProcessInfo::toHumanUnits(mem);
    EXPECT_TRUE(human.find("MBs") != std::string::npos);
  }
}

TEST(ProcessInfoTest, RSS) {
  if (ProcessInfo::isSupported(ProcessInfo::MEM_RSS)) {
    unsigned long rss = ProcessInfo::getMemoryResident();
    EXPECT_GT(rss, 0u);
  }
}
