// Copyright 2016, Tobias Hermann.
// https://github.com/Dobiasd/frugally-deep
// Distributed under the MIT License.
// (See accompanying LICENSE file or at
//  https://opensource.org/licenses/MIT)

#pragma once

#include "submodules/fdeep/common.hpp"

#include "submodules/fdeep/tensor3.hpp"

#include "submodules/fdeep/node.hpp"

#include <cstddef>
#include <memory>

namespace fdeep { namespace internal
{

class layer;
typedef std::shared_ptr<layer> layer_ptr;
typedef std::vector<layer_ptr> layer_ptrs;

class activation_layer;
typedef std::shared_ptr<activation_layer> activation_layer_ptr;
tensor3s apply_activation_layer(const activation_layer_ptr& ptr,
    const tensor3s& input);

class layer
{
public:
    explicit layer(const std::string& name)
        : name_(name), nodes_(), activation_(nullptr)
    {
    }
    virtual ~layer()
    {
    }

    void set_activation(const activation_layer_ptr& activation)
    {
        activation_ = activation;
    }

    void set_nodes(const nodes& layer_nodes)
    {
        nodes_ = layer_nodes;
    }

    virtual tensor3s apply(const tensor3s& input) const final
    {
        const auto result = apply_impl(input);
        if (activation_ == nullptr)
            return result;
        else
            return apply_activation_layer(activation_, result);
    }

    virtual tensor3 get_output(const layer_ptrs& layers,
        output_dict& output_cache,
        std::size_t node_idx, std::size_t tensor_idx) const
    {
        const node_connection conn(name_, node_idx, tensor_idx);

        if (!fplus::map_contains(output_cache, conn.without_tensor_idx()))
        {
            assertion(node_idx < nodes_.size(), "invalid node index");
            output_cache[conn.without_tensor_idx()] =
                nodes_[node_idx].get_output(layers, output_cache, *this);
        }

        const auto& outputs = fplus::get_from_map_unsafe(
            output_cache, conn.without_tensor_idx());

        assertion(tensor_idx < outputs.size(),
            "invalid tensor index");
        return outputs[tensor_idx];
    }

    std::string name_;
    nodes nodes_;

protected:
    virtual tensor3s apply_impl(const tensor3s& input) const = 0;
    activation_layer_ptr activation_;
};

inline tensor3 get_layer_output(const layer_ptrs& layers,
    output_dict& output_cache,
    const layer_ptr& layer,
    std::size_t node_idx, std::size_t tensor_idx)
{
    return layer->get_output(layers, output_cache, node_idx, tensor_idx);
}

inline tensor3s apply_layer(const layer& layer, const tensor3s& inputs)
{
    return layer.apply(inputs);
}

inline layer_ptr get_layer(const layer_ptrs& layers,
    const std::string& layer_id)
{
    const auto is_matching_layer = [layer_id](const layer_ptr& ptr) -> bool
    {
        return ptr->name_ == layer_id;
    };
    return fplus::throw_on_nothing(
        error("dangling layer reference: " + layer_id),
        fplus::find_first_by(is_matching_layer, layers));
}

} } // namespace fdeep, namespace internal
