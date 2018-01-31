#include <Core/Const/Constants.h>
#include <Core/Math/IntVec.h>
#include <iostream>

using namespace dem;

bool
IntVec::operator==(const IntVec v)
{
  return d_data[0] == v.d_data[0] && d_data[1] == v.d_data[1] && d_data[2] == v.d_data[2];
}

bool
IntVec::operator==(const int d)
{
  return d_data[0] == d && d_data[1] == d && d_data[2] == d;
}

bool
IntVec::operator!=(const IntVec v)
{
  return d_data[0] != v.d_data[0] || d_data[1] != v.d_data[1] || d_data[2] != v.d_data[2];
}

void
IntVec::operator+=(IntVec v)
{
  d_data[0] += v.d_data[0];
  d_data[1] += v.d_data[1];
  d_data[2] += v.d_data[2];
}

void
IntVec::operator-=(IntVec v)
{
  d_data[0] -= v.d_data[0];
  d_data[1] -= v.d_data[1];
  d_data[2] -= v.d_data[2];
}

void
IntVec::operator*=(int d)
{
  d_data[0] *= d;
  d_data[1] *= d;
  d_data[2] *= d;
}

IntVec
IntVec::operator+(IntVec v) const
{
  return IntVec(d_data[0] + v.d_data[0], d_data[1] + v.d_data[1], d_data[2] + v.d_data[2]);
}

IntVec
IntVec::operator-(IntVec v) const
{
  return IntVec(d_data[0] - v.d_data[0], d_data[1] - v.d_data[1], d_data[2] - v.d_data[2]);
}

IntVec IntVec::operator*(int d) const
{
  return IntVec(d_data[0] * d, d_data[1] * d, d_data[2] * d);
}

int IntVec::operator*(IntVec p) const
{
  return (d_data[0] * p.d_data[0] + d_data[1] * p.d_data[1] + d_data[2] * p.d_data[2]);
}

void
IntVec::print(std::ostream& ofs) const
{
  ofs << std::setw(OWID) << d_data[0] << std::setw(OWID) << d_data[1] << std::setw(OWID)
      << d_data[2];
}

IntVec
IntVec::fromString(const std::string& str)
{

  // Parse out the [num,num,num]
  std::string::size_type i1 = str.find("[");
  std::string::size_type i2 = str.find_first_of(",");
  std::string::size_type i3 = str.find_last_of(",");
  std::string::size_type i4 = str.find("]");

  std::string x_val(str, i1 + 1, i2 - i1 - 1);
  std::string y_val(str, i2 + 1, i3 - i2 - 1);
  std::string z_val(str, i3 + 1, i4 - i3 - 1);

  // Validate the values
  std::string validChars(" -0123456789");
  std::string::size_type xpos = x_val.find_first_not_of(validChars);
  std::string::size_type ypos = y_val.find_first_not_of(validChars);
  std::string::size_type zpos = z_val.find_first_not_of(validChars);

  if (xpos != std::string::npos || ypos != std::string::npos ||
      zpos != std::string::npos) {
    std::cerr << "Unexpected non-integer value in IntVector."
              << " Converting to int if possible\n";
  }

  // Create IntVec from string
  IntVec result;
  result.set(std::stoi(x_val), std::stoi(y_val), std::stoi(z_val));
  return result;
}

IntVec operator*(int d, IntVec v)
{
  return IntVec(v.x() * d, v.y() * d, v.z() * d);
}

IntVec
operator-(IntVec v)
{
  return IntVec(-v.x(), -v.y(), -v.z());
}

namespace dem {

std::ostream&
operator<<(std::ostream& os, const IntVec& v)
{
  os << ' ' << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
  return os;
}

} // end namespace dem