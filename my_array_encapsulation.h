/*
 * Created by kureii on 10/18/24.
 */

#pragma once
#include <cstdint>
#include <concepts>

template <typename T>
requires std::totally_ordered<T>
class myArrayEncapsulation {
 public:
  myArrayEncapsulation();
  explicit myArrayEncapsulation(uint64_t size);
  myArrayEncapsulation(const myArrayEncapsulation<T> &other);
  ~myArrayEncapsulation();

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
  myArrayEncapsulation<T>& operator=(const myArrayEncapsulation<T> &other);

 private:
  uint64_t size_;
  T* data_;
};

#include "my_array_encapsulation.tpp"
