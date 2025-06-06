// Copyright 2016, Tobias Hermann.
// https://github.com/Dobiasd/frugally-deep
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "submodules/fdeep/layers/layer.hpp"

namespace fdeep { namespace internal
{

class concatenate_layer : public layer
{
public:
    explicit concatenate_layer(const std::string& name)
        : layer(name)
    {
    }
protected:
    tensor3s apply_impl(const tensor3s& input) const override
    {
        return {concatenate_tensor3s(input)};
    }
};

} } // namespace fdeep, namespace internal
