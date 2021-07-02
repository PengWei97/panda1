//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GrainTracker.h"

/**
 * GrainTrackerCP derived class template to base objects on which maintain physical
 * parameters for individual grains.
 */
template <typename T>
class GrainDataTrackerCP : public GrainTracker
{
public:
  GrainDataTrackerCP(const InputParameters & parameters);

  /// return data for selected grain
  const T & getData(unsigned int grain_id) const;

protected:
  /// implement this method to initialize the data for the new grain
  virtual T newGrain(unsigned int new_grain_id) = 0;

  virtual void newGrainCreated(unsigned int new_grain_id);

  /// per grain data
  std::vector<T> _grain_data;
  // 如四阶张量-弹性模量
};

template <typename T>
GrainDataTrackerCP<T>::GrainDataTrackerCP(const InputParameters & parameters) : GrainTracker(parameters)
{
}
// 构造函数

template <typename T>
const T &
GrainDataTrackerCP<T>::getData(unsigned int grain_id) const
{
  mooseAssert(grain_id < _grain_data.size(), "Requested data for invalid grain index.");
  // 若grain_id >= _grain_data.size(),这断言
  // _grain_data.size()：数据结构的大小，即向量的分量
  // mooseAssert：如果它的条件返回错误，则终止程序执行
  // 如bicrystal：grain_id = 0,1 _grain_data.size() = 2
  return _grain_data[grain_id];
}
// getData在ComputePolycrystalElasticityTensorCP中被调用

template <typename T>
void
GrainDataTrackerCP<T>::newGrainCreated(unsigned int new_grain_id)
{
  if (_grain_data.size() <= new_grain_id)
    _grain_data.resize(new_grain_id + 1);

  _grain_data[new_grain_id] = newGrain(new_grain_id);
  // 在ComputePolycrystalElasticityTensorCP中重新定义
}
