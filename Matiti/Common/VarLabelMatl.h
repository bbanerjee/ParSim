#ifndef MATITI_VARLABELMATL_H
#define MATITI_VARLABELMATL_H

#include <Common/VarLabel.h>

namespace Matiti {

struct VarLabelMatl {
  VarLabelMatl(const VarLabel* label, int matlIndex)
    : label_(label), matlIndex_(matlIndex) {}
  VarLabelMatl(const VarLabelMatl& copy)
    : label_(copy.label_), matlIndex_(copy.matlIndex_)
  {}
  VarLabelMatl& operator=(const VarLabelMatl& copy)
  {
    label_=copy.label_; matlIndex_=copy.matlIndex_; 
    return *this;
  }
  
  bool operator<(const VarLabelMatl& other) const
  {
    if (label_->equals(other.label_)) {
      if (matlIndex_ == other.matlIndex_)
        return false;
      else
        return matlIndex_ < other.matlIndex_;
    }
    else {
      VarLabel::Compare comp;
      return comp(label_, other.label_);
    }
  };
  
  bool operator==(const VarLabelMatl& other) const
  {
    return ((label_->equals(other.label_)) && (matlIndex_ == other.matlIndex_));
  };
 
  const VarLabel* label_;
  int matlIndex_;
};  

} // End namespace Matiti

#endif
