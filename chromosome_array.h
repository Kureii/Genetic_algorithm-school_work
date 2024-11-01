/*
 * Created by kureii on 10/18/24.
 */

#pragma once
#include <cstdint>
#include <concepts>
#include "real_structures.h"

template <typename T>
requires std::totally_ordered<T>
class ChromosomeArray {
 public:
  ChromosomeArray();
  explicit ChromosomeArray(uint64_t size, uint64_t chromosome_length);
  explicit ChromosomeArray(uint64_t size, uint64_t chromosome_length,
      const mapping_structure_t& mapping_structure);
  ChromosomeArray(const ChromosomeArray &other);
  ~ChromosomeArray();

  [[nodiscard]] uint64_t size() const;

  void fill(const T& value);

  void reset();
  T* data();

  T* begin();
  T* end();
  const T* begin() const;
  const T* end() const;
  std::reverse_iterator<T*> rbegin();
  std::reverse_iterator<T*> rend();
  std::reverse_iterator<const T*> rbegin() const;
  std::reverse_iterator<const T*> rend() const;

  const T& operator[](std::size_t index) const;
  T& operator[](std::size_t index);
  ChromosomeArray& operator=(const ChromosomeArray& other);

  [[nodiscard]] uint64_t chromosome_length() const;
  [[nodiscard]] mapping_structure_t mapping_structure() const;

  void SetMappingStructure(const mapping_structure_t& mapping_structure);

 private:
  uint64_t size_;
  uint64_t chromosome_length_;
  mapping_structure_t mapping_structure_;
  T* data_;
};

#include "chromosome_array.tpp"
