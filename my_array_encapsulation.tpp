#pragma once

#include <stdexcept>

template <typename T>
requires std::totally_ordered<T>
myArrayEncapsulation<T>::myArrayEncapsulation() : size_(0), data_(nullptr) {
}

template <typename T>
requires std::totally_ordered<T>
myArrayEncapsulation<T>::myArrayEncapsulation(uint64_t size) : size_(size) {
  data_ = new T[size_];
}

template <typename T>
requires std::totally_ordered<T>
myArrayEncapsulation<T>::~myArrayEncapsulation() {
  delete[] data_;
}

template <typename T>
requires std::totally_ordered<T>
myArrayEncapsulation<T>::myArrayEncapsulation(const myArrayEncapsulation<T>& other) {
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
}

template <typename T>
requires std::totally_ordered<T>
myArrayEncapsulation<T>& myArrayEncapsulation<T>::operator=(const myArrayEncapsulation<T>& other) {
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
  return *this;
}

template <typename T>
requires std::totally_ordered<T>
uint64_t myArrayEncapsulation<T>::size() const {
  return size_;
}

template <typename T>
requires std::totally_ordered<T>
void myArrayEncapsulation<T>::reset() {
  fill(T{});
}

template <typename T>
requires std::totally_ordered<T>
void myArrayEncapsulation<T>::fill(const T& value) {
  for (uint64_t i = 0; i < size_; i++) {
    data_[i] = value;
  }
}

template <typename T>
  requires std::totally_ordered<T>
T* myArrayEncapsulation<T>::data() {
  return data_;
}


template <typename T>
requires std::totally_ordered<T>
T* myArrayEncapsulation<T>::begin() {
  return data_;
}
template <typename T>
requires std::totally_ordered<T>
const T* ::myArrayEncapsulation<T>::begin() const {
  return data_;
}


template <typename T>
requires std::totally_ordered<T>
T* ::myArrayEncapsulation<T>::end() {
  return data_ + size_;
}


template <typename T>
requires std::totally_ordered<T>
const T* myArrayEncapsulation<T>::end() const {
  return data_ + size_;
}


template <typename T>
requires std::totally_ordered<T>
const T& myArrayEncapsulation<T>::operator[](std::size_t index) const {
  if (index >= size_) {
    throw std::out_of_range("Index out of bounds");
  }
  return data_[index];
}

template <typename T>
requires std::totally_ordered<T>
T& myArrayEncapsulation<T>::operator[](std::size_t index) {
  if (index >= size_) {
    throw std::out_of_range("Index out of bounds");
  }
  return data_[index];
}
