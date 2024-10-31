/*
 * Created by kureii on 10/18/24.
 */

#pragma once
#include <cstdint>
#include <concepts>

template <typename T>
requires std::totally_ordered<T>
class ChromosomeArray {
 public:
  ChromosomeArray();
  explicit ChromosomeArray(uint64_t size, uint64_t chromosome_length);
  ChromosomeArray(const ChromosomeArray<T> &other);
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
  ChromosomeArray<T>& operator=(const ChromosomeArray<T>& other);

  [[nodiscard]] uint64_t chromosome_length() const;
 private:
  uint64_t size_;
  uint64_t chromosome_length_;

  T* data_;
};

#include "chromosome_array.tpp"
