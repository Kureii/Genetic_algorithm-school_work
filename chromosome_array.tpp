#pragma once

#include <stdexcept>

template <typename T>
  requires std::totally_ordered<T>
ChromosomeArray<T>::ChromosomeArray() : size_(0), data_(nullptr) {}

template <typename T>
  requires std::totally_ordered<T>
ChromosomeArray<T>::ChromosomeArray(uint64_t size, uint64_t chromosome_length) : size_(size), chromosome_length_(chromosome_length) {
  data_ = new T[size_];
}

template <typename T>
  requires std::totally_ordered<T>
ChromosomeArray<T>::~ChromosomeArray() {
  delete[] data_;
}

template <typename T>
  requires std::totally_ordered<T>
ChromosomeArray<T>::ChromosomeArray(
    const ChromosomeArray<T>& other) {
  if (this != &other) {
    if (other.size_ > 0) {
      size_ = other.size_;
      data_ = new T[size_];
      for (uint64_t i = 0; i < size_; i++) {
        data_[i] = other.data_[i];
      }
    } else {
      delete[] data_;
      data_ = nullptr;
    }
  }
  size_ = other.size_;
  chromosome_length_ = other.chromosome_length_;
}

template <typename T>
  requires std::totally_ordered<T>
ChromosomeArray<T>& ChromosomeArray<T>::operator=(
    const ChromosomeArray<T>& other) {
  if (this != &other) {
    if (other.size_ > 0) {
      if (size_ != other.size_) {
        delete[] data_;
        size_ = other.size_;
        data_ = new T[size_];
      }
      for (uint64_t i = 0; i < size_; i++) {
        data_[i] = other.data_[i];
      }
    } else {
      delete[] data_;
      data_ = nullptr;
    }
  }
  chromosome_length_ = other.chromosome_length_;
  return *this;
}

template <typename T>
  requires std::totally_ordered<T>
uint64_t ChromosomeArray<T>::size() const {
  return size_;
}

template <typename T>
  requires std::totally_ordered<T>
void ChromosomeArray<T>::reset() {
  fill(T{});
}

template <typename T>
  requires std::totally_ordered<T>
void ChromosomeArray<T>::fill(const T& value) {
  for (uint64_t i = 0; i < size_; i++) {
    data_[i] = value;
  }
}

template <typename T>
  requires std::totally_ordered<T>
T* ChromosomeArray<T>::data() {
  return data_;
}

template <typename T>
  requires std::totally_ordered<T>
T* ChromosomeArray<T>::begin() {
  return data_;
}
template <typename T>
  requires std::totally_ordered<T>
const T* ::ChromosomeArray<T>::begin() const {
  return data_;
}

template <typename T>
  requires std::totally_ordered<T>
T* ::ChromosomeArray<T>::end() {
  return data_ + size_;
}

template <typename T>
  requires std::totally_ordered<T>
const T* ChromosomeArray<T>::end() const {
  return data_ + size_;
}

template <typename T>
  requires std::totally_ordered<T>
std::reverse_iterator<T*> ChromosomeArray<T>::rbegin() {
  return std::reverse_iterator<T*>(end());
}

template <typename T>
  requires std::totally_ordered<T>
std::reverse_iterator<T*> ChromosomeArray<T>::rend() {
  return std::reverse_iterator<T*>(begin());
}

template <typename T>
  requires std::totally_ordered<T>
std::reverse_iterator<const T*> ChromosomeArray<T>::rbegin() const {
  return std::reverse_iterator<const T*>(end());
}

template <typename T>
  requires std::totally_ordered<T>
std::reverse_iterator<const T*> ChromosomeArray<T>::rend() const {
  return std::reverse_iterator<const T*>(begin());
}

template <typename T>
  requires std::totally_ordered<T>
const T& ChromosomeArray<T>::operator[](std::size_t index) const {
  if (index >= size_) {
    throw std::out_of_range("Index out of bounds");
  }
  return data_[index];
}

template <typename T>
  requires std::totally_ordered<T>
T& ChromosomeArray<T>::operator[](std::size_t index) {
  if (index >= size_) {
    throw std::out_of_range("Index out of bounds");
  }
  return data_[index];
}

template <typename T>
  requires std::totally_ordered<T>
uint64_t ChromosomeArray<T>::chromosome_length() const {
  return chromosome_length_;
}
