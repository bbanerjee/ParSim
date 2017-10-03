#ifndef ELLIP3D_BOUNDARY_READER_H_H
#define ELLIP3D_BOUNDARY_READER_H_H

#include <Boundary/BoundaryContainers.h>
#include <Core/Geometry/Box.h>
#include <Core/Types/realtypes.h>

#include <InputOutput/zenxml/xml.h>

namespace dem {

class BoundaryFileWriter
{

public:
  BoundaryFileWriter() = default;
  ~BoundaryFileWriter() = default;

  void writeCSV(std::size_t boundaryNum,
                const std::string& outputFileName, 
                Box& allContainer) const;

  void writeXML(std::size_t boundaryNum,
                const std::string& outputFileName, 
                Box& allContainer) const;

private:
  enum BoundaryFlag
  {
    ONLY_BOTTOM_BOUNDARY,
    NO_TOP_BOUNDARY,
    ALL_BOUNDARIES
  };

  enum BoundaryID
  {
    NONE,
    XMINUS,
    XPLUS,
    YMINUS,
    YPLUS,
    ZMINUS,
    ZPLUS
  };

  BoundaryFlag getBoundaryFlag(std::size_t boundaryNum) const
  {
    if (boundaryNum == 1) 
      return BoundaryFlag::ONLY_BOTTOM_BOUNDARY;
    else if (boundaryNum == 5)
      return BoundaryFlag::NO_TOP_BOUNDARY;
    else
      return BoundaryFlag::ALL_BOUNDARIES;
  }

  inline void writeHeader(std::size_t extraBoundaries, std::ofstream& ofs) const;
  inline void writeBoundaryID(BoundaryID id, std::ofstream& ofs) const;
  inline void writeEmptyLine(std::ofstream& ofs) const;
  inline void writeXMinusCSV(const Vec& domainMin, const Vec& domainMax,
                             const Vec& domainCen, std::ofstream& ofs) const;
  inline void writeXPlusCSV(const Vec& domainMin, const Vec& domainMax,
                            const Vec& domainCen, std::ofstream& ofs) const;
  inline void writeYMinusCSV(const Vec& domainMin, const Vec& domainMax,
                             const Vec& domainCen, std::ofstream& ofs) const;
  inline void writeYPlusCSV(const Vec& domainMin, const Vec& domainMax,
                            const Vec& domainCen, std::ofstream& ofs) const;
  inline void writeZMinusCSV(const Vec& domainMin, const Vec& domainMax,
                             const Vec& domainCen, std::ofstream& ofs) const;
  inline void writeZPlusCSV(const Vec& domainMin, const Vec& domainMax,
                            const Vec& domainCen, std::ofstream& ofs) const;
  inline void writeExtraZPlusCSV(const Vec& domainMin, const Vec& domainMax,
                                 const Vec& domainCen, std::ofstream& ofs) const;

  inline void writeHeader(std::size_t extraBoundaries, zen::XmlOut& xml) const;
  inline void writeBoundaryID(BoundaryID id, zen::XmlOut& xml) const;
  inline void writeXMinusXML(const Vec& domainMin, const Vec& domainMax,
                             const Vec& domainCen, zen::XmlOut& xml) const;
  inline void writeXPlusXML(const Vec& domainMin, const Vec& domainMax,
                            const Vec& domainCen, zen::XmlOut& xml) const;
  inline void writeYMinusXML(const Vec& domainMin, const Vec& domainMax,
                             const Vec& domainCen, zen::XmlOut& xml) const;
  inline void writeYPlusXML(const Vec& domainMin, const Vec& domainMax,
                            const Vec& domainCen, zen::XmlOut& xml) const;
  inline void writeZMinusXML(const Vec& domainMin, const Vec& domainMax,
                             const Vec& domainCen, zen::XmlOut& xml) const;
  inline void writeZPlusXML(const Vec& domainMin, const Vec& domainMax,
                            const Vec& domainCen, zen::XmlOut& xml) const;
  inline void writeExtraZPlusXML(const Vec& domainMin, const Vec& domainMax,
                                 const Vec& domainCen, zen::XmlOut& xml) const;

  BoundaryFileWriter(BoundaryFileWriter const&) = delete;
  void operator=(BoundaryFileWriter const&) = delete;
};
}
#endif
