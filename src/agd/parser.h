#pragma once

#include "format.h"
#include "errors.h"
#include "buffer.h"
#include <vector>
#include <array>
#include <string>
#include <cstdint>

namespace agd {
  extern unsigned char nst_nt4_table[256];

  template <size_t N>
  class BaseMapping {

  private:
    std::array<char, N> characters_;
    std::size_t effective_characters_;

  public:

    BaseMapping(std::array<char, N> chars, std::size_t effective_characters) : characters_(chars),
      effective_characters_(effective_characters) {}

    BaseMapping() {
      characters_.fill('\0');
      characters_[0] = 'Z'; // TODO hack: an arbitrary bad value, used to indicate an impossible issue
    }

    const std::array<char, N>& get() const {
      return characters_;
    }

    const std::size_t effective_characters() const {
      return effective_characters_;
    }
  };

  const BaseMapping<3>*
  lookup_triple(const std::size_t bases);

  class RecordParser
  {
  public:
    explicit RecordParser();

    Status ParseNew(const char* data, const std::size_t length, const bool verify, Buffer *result_buffer, 
        uint64_t *first_ordinal, uint32_t *num_records, std::string &record_id, bool unpack=true);

  private:

    void reset();

    Buffer conversion_scratch_, index_scratch_;
    const format::RelativeIndex *records = nullptr;
  };

}  //  namespace agd
