
#pragma once
#include "agd_record_reader.h"
#include "json.hpp"

namespace agd {
using json = nlohmann::json;

// AGD Dataset loads an entire dataset into memory.
// This may not be suitable for very large datasets.
// Consider adding AGDBufferedDataset, which will buffer one or two chunks
// from each requested column for sequential access.
class AGDDataset {
 public:
  // create dataset
  // load only specific columns if given
  static Status Create(const std::string& agd_json_path,
                       std::unique_ptr<AGDDataset>& dataset,
                       std::vector<std::string> columns = {});

  // allows iteration over a complete agd column
  // associated AGDDataset must outlive the ColumnIterator
  class ColumnIterator {
    friend class AGDDataset;

   public:
    ColumnIterator() = default;
    Status GetNextRecord(const char** data, size_t* size);
    Status GetNextAt(size_t index, const char** data, size_t* size);
    void Reset();

   private:
    ColumnIterator(std::vector<AGDRecordReader>* readers, uint32_t total_records)
        : column_readers_(readers), total_records_(total_records) {}

    std::vector<AGDRecordReader>* column_readers_;
    uint32_t current_reader_ = 0;
    uint32_t total_records_;
  };

  Status Column(const std::string& column, ColumnIterator* iter) {
    if (column_map_.find(column) != column_map_.end()) {
      *iter = ColumnIterator(&column_map_[column], total_records_);
      return Status::OK();
    } else {
      return NotFound("column ", column, " not found in column map.");
    }
  }

  uint32_t Size() const { return total_records_; }

  const std::string& Name() const { return name_; }

 private:
  AGDDataset() = default;
  Status Initialize(const std::string& agd_json_path,
                    const std::vector<std::string>& columns);
  std::string name_;

  json agd_metadata_;

  std::vector<Buffer> chunks_;
  // map column name to its set of readers
  std::unordered_map<std::string, std::vector<AGDRecordReader>> column_map_;
  std::vector<size_t> chunk_sizes_;
  uint32_t total_records_;
};

}  // namespace agd
